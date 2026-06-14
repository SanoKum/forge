# 非平衡凝縮 (4 モーメント方程式) の forge 実装

## メタ

- **area**: `condensation`
- **status**: `in_progress`  <!-- Phase 1 実装中。Phase 2/3 は planned -->
- **related_docs**:
  - `docs/condensation/theory.md`
  - `docs/condensation/implementation.md`
- **related_plans**:
  - `.github/plans/condensation-nonequilibrium-session-prompt.md` (着手用プロンプト・コンテキスト)
  - `.github/plans/thermophysics-multicomponent-tpgas.md` (多成分 TP gas、気相多成分化の前提)
  - `.github/plans/architecture-rans-sst.md` (追加保存スカラー+ソースの骨格の手本)
- **created**: `2026-06-14`
- **owner**: `CFD Dev`

## 1. 目的

極超音速ノズルの強膨張で試験ガス (N2 等) が非平衡凝縮し潜熱放出で静圧・静温が上がる現象を、
気相 Euler/NS に**液相のモーメント輸送方程式 4 本** ($\rho g,\rho Q_2,\rho Q_1,\rho Q_0$) を結合して
計算できるようにする。最終的に case/34 Arthur ノズルで dry に対し凝縮による静圧上振れ (論文
Lin 2014 Fig.11) を再現する。

## 2. スコープ

- **やる (本 plan 全体)**: CNT+Iland 核生成 / 修正 Gyarmathy(Goodheart) 成長 / Hill $T_d$ /
  method of moments / 二相 EOS / fractional-step + point-implicit ソース。凝縮種ごとに 4 モーメント +
  物性を持ち、核生成/成長/表面張力モデルを種ごとに config 切替 (N2=CNT_Iland+Goodheart /
  H2O=CNT_Kantrowitz+Hertz–Knudsen)。
- **Phase 分け**:
  - **Phase 1 (本セッション)**: 4 モーメントを受動スカラー (ソース=0) として輸送する骨格 +
    case/34 dry 回帰一致。
  - **Phase 2 (planned)**: N2 凝縮物理 + 二相 EOS + point-implicit ソース + $\mu_n$ 無次元化 →
    case/34 で静圧上振れ検証。
  - **Phase 3 (planned)**: H2O モデル + case/16 Wyslouzil 検証。
- **やらない (当面)**: Roe の一般 EOS 固有構造対応 (後段)、気相+モーメント密結合 block (収束不良時のみ)、
  double 化 (float の桁落ちが顕在化したらフォールバック)。

## 3. 関連 docs と前提

理論は [`docs/condensation/theory.md`](../../docs/condensation/theory.md)、実装対応は
[`docs/condensation/implementation.md`](../../docs/condensation/implementation.md)。物理係数の出典は
[`papers/.../_summary.md`](../../papers/on%20nitrogen%20condensation%20in%20hypersonic%20nozzle%20flows_summary.md)。
前提: case/34 で dry 膨張再現済み ([case/34 README](../../case/34.arthur_n2_nozzle/README.md))。

## 4. 設計方針 (論点 A–D)

- **(A) 多成分一般化**: モーメントは種インデックス付きで確保し `nCondSpecies` でループ。当面 N2 1 種。
  核生成/成長/表面張力/飽和圧は enum+switch (device) + 種ごと係数構造体で切替。
- **(B) 対流**: 受動スカラー、SLAU 優先、既存 `ScalarTransportDesc` 流用。二相 EOS の圧力 g 依存は
  Phase 2。Roe 後段。
- **(C) 連成・陰解法**: NS と**分離 (segregated)**。fractional-step + point-implicit ソース。
  二相逆結合は loose coupling (T/P 再計算)。密結合は収束不良時のフォールバック。
- **(D) 精度**: まず全 float + $\mu_n$ 無次元化で O(1) 化。double 化はフォールバック。

二相 EOS の圧力微分 $\kappa=(\rho-\rho g)R/C_v$ (=$\gamma-1$ 相当)、
$\xi_g=\partial p/\partial(\rho g)=-RT+\kappa(L-RT)$ (新規列)、flux Jacobian への入り方は
[theory.md](../../docs/condensation/theory.md) 5 節 / [implementation.md](../../docs/condensation/implementation.md) 3 節を参照。

## 5. 実装ステップ

### Phase 1 (本セッション)

1. **config**: `solverConfig.{hpp,cpp}` に `condensation` (0/1)・`nCondSpecies` (既定 0) を追加。
2. **変数登録**: `variables.{hpp,cpp}` に `registerCondensation(nCondSpecies)` を新設。種ごとに
   `rog/roQ2/roQ1/roQ0` と N/M・残差・point-implicit 対角・primitive を追加 (roK/roOmega 構成踏襲)。
3. **移流 (新規)**: `cuda_forge/condensationTransport_d.{cu,cuh}` — `ScalarTransportDesc` を種×4 本
   構築し受動スカラー移流 (diffusion=0, src_jac=0)。residual/timeIntegration wrapper。
4. **境界 (新規)**: `cuda_forge/condensationBoundary_d.{cu,cuh}` — inlet=0 / 他 zero-grad。
5. **更新・残差・初期**: `update_d.cu` に N/M ステージ + `applyCondensationPointImplicit`、`main.cpp`
   に wrapper 呼出・残差列追加・register 呼出、`setInitial.hpp` に H2D 既定 0・restart フォールバック。
6. **検証**: case/34 run_0003 複製で dry 回帰一致。

### Phase 2 (N2 凝縮物理) — 分解と確定事項

1. **N2 物性 device モジュール** (`condensationProperties_d.cuh`) — 実装済。$p_{sat},\rho_l,L,\sigma$ を
   double で評価。**相は過冷却液 (supercooled liquid)、固体ではない** (論文確認: "liquid phase",
   "liquid droplets"、$r_*$ は $\rho_l$ 使用。固相フィットは図比較用で不使用)。液フィットは ~40K 以下で
   外挿破綻するため**物性評価温度を $[45\text{K},T_c)$ にクランプ** (`COND_T_PROP_FLOOR`)。anchors 検証済
   (NBP/三重点/液密度/潜熱/σ が ~1%)。
2. **一温度 二相 EOS** — 実装済。$T_v=T_d=T$、$e=c_v T-g L(T)$ から Newton で $T$ 逆算、
   $p=(1-g)\rho R T$、$L(T)$ 温度依存あり。$T_d$ は輸送変数にしない。
   - **初期実装 = CPG** (`thermalMethod 0` 分岐, `cond_T_from_e_cpg`)。論文が完全理想気体なのに合わせ、
     case/34 を `thermalMethod 0` のまま使える (TP 初期化不要)。**これを既定で進める**。
   - **後続 B = thermally perfect** (`thermalMethod 2` 分岐, `cond_T_from_e_onetemp`, NASA-9 $e_v$)。実装済。
   - host unit test: CPG/TP とも (b) g=0 厳密縮約 (diff 0)、(c) g>0 安定 (潜熱で T 50→72K, 残差 1e-9)。
3. **核生成 $J$** (CNT × Iland 補正)、**成長 $dr/dt$** (Goodheart)、**Hill $T_d$** 反復 — device 関数。
4. **ソース kernel** $S_{Q_n},S_g$ + point-implicit `src_jac` + $\mu_n$ 無次元化 (float)。
5. **検証**: case/34 で貯気を振り中心線静圧が dry 等エントロピー線より上振れ (論文 Fig.11)。

### Phase 3 (planned)

H2O モデル (CNT_Kantrowitz + Hertz–Knudsen) + Wyslouzil 検証。詳細は implementation.md 4–8 節。

## 6. 検証

- **ビルド**: `cmake --build solver_density_cuda/build-native -j`。
- **検証ケース**: `case/34.arthur_n2_nozzle/` (dry 回帰)。Phase 2 で N2 凝縮、Phase 3 で
  `case/16.nozzle_wys` (H2O)。
- **判定基準 (Phase 1)**: ① NaN/Inf 無し (序盤+最終)、② `tools/check_convergence.py` VERDICT=PASS、
  ③ 中心線 P/P0・出口 Mach が run_0003 と一致。

## 7. 影響範囲

- 触るファイル: `input/solverConfig.{hpp,cpp}`, `variables.{hpp,cpp}`,
  `cuda_forge/condensationTransport_d.*` (新規), `cuda_forge/condensationBoundary_d.*` (新規),
  `cuda_forge/update_d.cu`, `main.cpp`, `input/setInitial.hpp`。
- 既存ケースへの影響: `condensation` 既定 0 で既存経路はビット不変。
- docs: `docs/condensation/{theory,implementation}.md`, `docs/index.md`。

## 8. 完了条件

- [x] `docs/condensation/theory.md` 作成
- [x] `docs/condensation/implementation.md` 作成
- [x] Phase 1 実装・検証完了 (§6)
- [x] `.github/plans/README.md` 状態同期
- [ ] Phase 2/3 着手時に本 plan を更新

## 9. 変更ログ

- `2026-06-14` — 初稿。docs (theory/implementation) 作成、Phase 1 (受動スカラー骨格) 着手。
- `2026-06-14` — **Phase 2 着手**: 方針修正 (相=過冷却液 not 固体 / 一温度 $T_v=T_d=T$ /
  $T_d$ は輸送変数にしない / 2 温度・thermally perfect は後続)。N2 物性 device モジュール
  (`condensationProperties_d.cuh`, 液フィット 45K クランプ, anchors 検証) と一温度 二相 EOS
  (`condensationEOS_d.cuh`) を実装。**初期実装は CPG (`thermalMethod 0`, `cond_T_from_e_cpg`)** に決定
  (論文が完全理想気体、case/34 を thermalMethod 0 のまま使え TP 初期化不要)。TP 版 (`thermalMethod 2`,
  `cond_T_from_e_onetemp`) も実装済 (後続 B)。両分岐とも host unit test で g=0 厳密縮約・g>0 安定を確認。
  `L(T)` 温度依存あり。次: 核生成/成長ソース → case/34 (CPG) で静圧上振れ検証。
- `2026-06-14` — **Phase 1 実装・検証完了**。`condensation`/`nCondSpecies` config、`registerCondensation`
  (種ごと ρg,ρQ2,ρQ1,ρQ0 + N/M/残差/point-implicit 対角/primitive)、`condensationTransport_d.{cu,cuh}`
  (ScalarTransportDesc 流用の受動スカラー移流・境界 inlet=0/他 zero-grad・更新・時間積分) を実装。main.cpp の
  assembleResidual/explicit RK/implicit/dual-time に wrapper を配線、残差列 (rms_rog_0 等) と restart 0
  フォールバックを追加。native ビルド成功。**検証 (case/34 run_0004 vs run_0005 cond OFF, 12000 step)**:
  液相モーメント全セル厳密 0・NaN/Inf なし・出口 M=6.874 で run_0003 と一致。cond ON/OFF の場差
  (ro maxrel 9.5e-4) は cond OFF 2 回の run-to-run ノイズ (1.0e-3, Venkata リミットサイクル+GPU atomic
  非決定) 以下=凝縮スカラーは気相に不干渉。Phase 2 (N2 凝縮物理) へ。

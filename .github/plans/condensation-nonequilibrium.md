# 非平衡凝縮 (4 モーメント方程式) の forge 実装

## メタ

- **area**: `condensation`
- **status**: `in_progress`  <!-- Phase 1/2/3 実装・検証済 (N2 Fig.2, H2O Wyslouzil Fig.3)。残: 二温度モデル / TP carrier 低温対応 / μ_n 無次元化は future -->
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

### Phase 3 (H2O / Wyslouzil, carrier+condensible) — 着手

**系**: N2 キャリア(非凝縮)+ 水蒸気(凝縮、希薄)。case/16 ノズル(M≈2, 整備済)。検証=Wyslouzil
JCP 113,7317 (2000) Fig.3 の 1kPa ケース(p0=59.07kPa, T0=286.65K, pv0=1.0kPa, y=1.7%, 水~1.1%質量)。
digitize: `case/16.nozzle_wys/wyslouzil_fig3_1kPa.csv`。

**設計 (Option A 完全多成分, ユーザ確定)**: 化学種は **N2 と「総水(蒸気+液)」** を保存
(thermalMethod 2, ΣY_s=1 不変、既存化学種輸送そのまま)。**液は moment `rog`、蒸気 = ρY_H2O − ρg**。
凝縮は総水保存のまま蒸気→液を移す → **化学種側に sink 不要**(最小結合)。EOS/ソースで蒸気=Y_H2O−g
を使い、過飽和 S=pv/psat(T)、pv=(水蒸気モル分率)·p。realizability: g≤Y_H2O。

**実装ステップ**:
1. **H2O thermo (NASA-9) を thermo DB に追加** — 済。
2. **水凝縮物性** (`h2o_psat`=Murphy-Koop 過冷却水, `h2o_rho_cond/latent/sigma`, model 切替) — 済 (anchors 検証)。
3. **carrier+condensible EOS** (dependentVariables thermalMethod 2): 蒸気=Y_H2O−g で気相混合圧力、
   潜熱 g·L_H2O を含む T 逆算。
4. **ソース結合**: 蒸気分圧から S、CNT+Kantrowitz 核生成・Hertz–Knudsen 成長 (enum で N2 と切替)、g≤Y_H2O。
5. **エネルギー流束** 二相 (水潜熱・気相混合) — SLAU の補正を carrier 系へ一般化。
6. case/16 に H2O run、Fig.3 と比較。

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

- `2026-06-15` — **核生成/成長モデルの感度スイッチ追加 (Kantrowitz / Gyarmathy)**。ユーザ要望で、
  核生成の **Kantrowitz 非等温補正**と **Gyarmathy 熱伝導律速成長**を config フラグ化:
  - `condKantrowitz` (0/1): $J\to J/(1+\theta)$, $\theta=\frac{2(\gamma-1)}{\gamma+1}b(b-\tfrac12)$, $b=L/(R_vT)$。
  - `condGrowthModel` (0=Hertz-Knudsen/Goodheart, 1=Gyarmathy): Gyarmathy は N2 Goodheart と前因子を共有し
    Knudsen 内挿を $1/(1+3.18Kn)$ にした熱伝導律速形。carrier では Kn 用平均自由行程に全圧を使用。
  - 実装: `condensationSource_d.{cuh,cu}` の `cond_nucleation`/`cond_growth`/`cond_source_vector` に
    引数追加 (既定 off でビット不変)、`solverConfig.{hpp,cpp}` にフラグ。docs/condensation/implementation.md 更新。
  - 検証: case/16 SST で 2×2 (Kantrowitz off/on × HK/Gyarmathy) = run_0010〜0013 を比較。
- `2026-06-15` — **Phase 3 H2O / Wyslouzil Fig.3 検証成功 + carrier=CPG へ方針変更**。
  - **発散原因究明**: 当初の Option A (thermalMethod 2 / NASA-9 TP で N2+H2O) は、ノズル膨張で気相が
    **<200K (NASA-9 フィット下限) まで冷却**し外挿が不安定→初手から発散 (T が 50K フロアに張り付き P が
    入口の 8 倍へ発散)。dry CPG (run_0001) は同条件で T 163–293K に安定収束することを確認し、**低温域では
    N2 cp がほぼ一定で CPG (γ=1.4) が正確かつ頑健**と判断。**carrier を CPG で解く方針に変更** (H2O は
    `species` として移流し vapor budget に使用、潜熱・物性は CPG エネルギー経路へ組込み)。
  - **CPG carrier path の H2O 対応**: `cond_T_from_e_cpg` / dependentVariables CPG 分岐 / 二相エネルギー流束
    (SLAU・HLLE・ROE) の潜熱を `n2_latent` 固定から **`cond_latent(condModel)` dispatch** へ一般化
    (N2/H2O 切替)。CPG では per-cell `cp`/`Rmix` 配列が未充填のため、ソースへは `nullptr` を渡し
    `cfg.cp`/γ フォールバック (=carrier N2 の cp/R) を使用。
  - **実現可能性クランプ** (`cond_realizability_clamp_d`): 陰解法の実効 Δt と θ 律速の dt_local 差で僅かに
    過凝縮 (g/Y1≈1.07) しうるため、毎ステップ `0≤rog≤roY_w` (carrier) と `roQn≥0` を硬クランプ。
    クランプ後 g/Y1∈[0,1] で physics 不変。
  - **検証 (case/16, Wyslouzil Fig.3 pv0=1kPa)**: `run_0007`(dry)/`run_0008`(cond) 非粘性。凝縮潜熱で
    中心線 T 143→185K, p/p0 が dry 比 **1.19〜1.28 倍** (2.5–8cm) に上昇し、Wyslouzil 実験の
    cond/isentrope 比 1.18〜1.23 と **~5% 以内で一致**。`run_0009`(dry)/`run_0010`(cond) で粘性+SST 乱流版も
    実施 (ユーザ要望)。Y_H2O=0.01095 は全域保存、全エンタルピー H=cpT−gL+ek 保存。
- `2026-06-14` — **後処理 3 件 (1 各スキーム展開 / 2 dry baseline / 3 低優先トリオ)**。
  - **1**: 二相エネルギー流束補正を HLLE/ROE へ展開し各スキームで case/34 を実行。全て全エンタルピー保存
    (dev ±0.0%)・過剰凝縮解消 (g mean 0.8–1.1%)・壁圧上昇あり (SLAU 1.33/1.43x = 論文最良、ROE 1.28/1.28x、
    HLLE 1.20/1.25x = 最散逸)。AUSM/KEEP (dispatch 外) と TP 分岐は未対応 (CPG で完了)。
  - **2 (dry baseline ① 解明)**: P/P0 を A/A* で比較すると **forge dry = 等エントロピー (非粘性) に厳密一致**、
    **論文 dry は等エントロピーより ~1.35x 上** (実験はさらに上)。つまり forge は正しく、論文 dry/実験は粘性
    境界層変位ぶん上振れ。**cond/dry 比が一致**するので凝縮物理は正。絶対値一致には粘性計算が要る (別途)。
  - **3 (低優先トリオ評価)**: (a) チェッカーボードは**気相 odd-even** (dry でも sign-change 0.79、凝縮セルの
    モーメントもこれを継承) で**モーメント精度問題ではない** → 無次元化は効かず不要 (float で g~1%・r̄ 物理)。
    (b) 二相 BC は flux 修正後 H が境界でも保存 (inlet +0.01%/wall −0.00%/outlet −0.01%) で実害なし。
    (c) src_jac 自己抑制線形化は実装済 (commit 7cfdb84)。**残る実問題は気相 odd-even (pre-existing, 別件)**。
- `2026-06-14` — **N2 検証成功 (Fig.2 一致)**。全エンタルピー診断 (ユーザ提案) で **SLAU エネルギー流束が
  単相エンタルピー $h=c_p(1-g)T+ek$ を運び潜熱 $-gL$ を落としていた**ことが判明 (dry は H 保存 0.0%、
  凝縮で −8.7%/局所 −16% の非保存)。このエネルギー漏れが過冷却→過飽和維持→過剰凝縮 (g~0.34) の根本。
  **修正**: SLAU_d CPG 分岐で $h$ を二相全エンタルピー $h=c_pT-gL+ek$ に補正 (差 $g(c_pT-L)$、セル値 1 次、
  $g$=nullptr でビット不変)。結果: 全エンタルピー保存回復 (dev −0.0%)、**過剰凝縮解消 (g max 0.34→0.08,
  mean 1.0%)**、**壁面静圧の凝縮上昇が論文 Fig.2 と 1–2% 一致** (cond/dry 比 3in 1.21x/5in 1.43x vs 論文
  1.20/1.45x、onset バンプ再現)。Fig.2 は `case/34/arthur_fig2_digitized.csv`。**残**: 同バグを
  HLLE/ROE/AUSM/TP 分岐へ展開 (case/34 は SLAU+CPG で完了)、dry baseline ズレ (①) 切り分け。


- `2026-06-14` — 初稿。docs (theory/implementation) 作成、Phase 1 (受動スカラー骨格) 着手。
- `2026-06-14` — **Phase 2 ソース項 初期実装 (bounded 確認)**。核生成 (CNT+Iland)・成長 (Goodheart) を
  `condensationSource_d.{cu,cuh}` に実装、一温度・rates freeze・明示的ソース+安定化 clamp/limiter
  ($J$ 上限, $dr/dt\!<\!0\to0$, $\bar r\le r_*$ 成長停止, 非負, $g\le1$, 1step $\Delta g/\Delta T/$蒸気枯渇を
  $\theta$ 律速)。src_jac はソース由来 0 (時間+移流対角のみ; 完全自己抑制線形化は後続)。**psat を C-C 外挿**
  (45K クランプだと過飽和が潰れ核生成せず → C-C で物理外挿し onset 復活)。**case/34 run_0006 (CPG, dry
  restart)**: NaN なし・bounded ($g\in[0,1]$)・最小 T 27.8→34.6K (潜熱)・g=0 域は dry 一致。**ユーザ提示の
  初期目標 (NaN なし/g=0 単相復帰/bounded) を達成**。g 過大 (22%) は run_0007 (一様初期から develop) でも
  同様 (g~21%) → restart アーティファクトでなく定常解の過大予測 (定常で S~100 のまま、レート過大)。
  **定量一致は後続課題** (paper 壁圧トレース照合・レート較正・完全 src_jac 線形化・膨張冷却とのバランス)。
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

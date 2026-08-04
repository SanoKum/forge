# 軸対称の SU2 流定式化 (planar 幾何 + 1/y ソース項) — `axisymMethod: 1`

## メタ

- **area**: `axisymmetric / discretization`
- **status**: `in-progress`
- **related_docs**:
  - `methods/axisymmetric/implementation.md` (現行 B 流儀 = r 重み幾何の正本)
- **related_plans**:
  - [`boundary-node-nozzle-wall-outlet-stability.md`](boundary-node-nozzle-wall-outlet-stability.md)
    (§2.6: r 重み幾何の軸半 CV 悪条件化 → `nodeAxisDirichlet` 対症の経緯)
- **created**: `2026-08-04`
- **owner**: `sano` (指示)・Claude (実装)

## 1. 目的

forge の軸対称は r 重み幾何 (B 流儀: 面積・体積に r を乗じる) で、node (median-dual) では
軸上ノードの半 CV 体積・面積が r̄ ∝ Δr で消え、軸行が radial 圧力平衡からデカップルして
真空まで過膨張する基底欠陥がある (対症 = `nodeAxisDirichlet`)。ユーザ指示
「軸の処理は全部 SU2 を真似する」に従い、**SU2 と同じ「planar 幾何 + 1/y ソース項」方式**を
`axisymMethod: 1` として実装する。軸ノードは通常 DOF として解く (Dirichlet 置換なし)。

## 2. SU2 の定式 (.external/su2-src で確認済み)

- **幾何**: 2D planar のまま (体積・面積に r を乗じない)。
- **非粘性ソース** (`CSourceAxisymmetric_Flow::ComputeResidual`): y>EPS で
  LinSysRes += (V/y)·[ρv, ρuv, ρv², ρvH] (flux-out 側に加算 = 物理ソースは負号)。
  y≤EPS (軸上) は**ソース 0**。解析 Jacobian (4×4) を対角に加算。
- **粘性ソース** (`ResidualDiffusion`): AuxVar A0=μ_tot·v/y, A1=A0·v, A2=A0·u の GG/LSQ 勾配を
  使い、運動量 2 成分とエネルギーに 1/y 粘性項 + (2/3)∇A 項。κ_tot = κ + cp·μt/Pr_t (0.9)。
  ※ SU2 のエネルギー項中 `+ρk` (turb_ke) は次元不整合が疑われるため移植しない (!rans 相当)。
- **SST ソース** (`ResidualAxisymmetricConvectionDiffusion`): res_φ −= (V/y)·(ρv·φ −
  (μ+σ_φ μt)·∂φ/∂y)、対角 Jacobian −= (V/y)·ρv (φ = k, ω)。
- **軸 BC**: `BC_Sym_Plane` = 鏡映 flux + 法線運動量残差の射影 — forge の
  slip ghost + `zeroAxisRadialResidual` が既に等価。

## 3. forge への移植設計

- **config**: `physProp.axisymMethod` (int, 既定 0 = 現行 r 重み・ビット不変 / 1 = SU2 流)。
- **幾何** (`variables.cpp`): method 1 では r 重み付けブロックをスキップ (planar のまま)。
  `A_planar = volume` を入れておく (r_eff 系の既存式は method 1 では使わない)。
- **幾何項** (`axisymmetricGeomTerms_d`): method 1 の `axisym_uy_over_r` は
  `Uy/ccy` (軸ノード = `axis_flag` は 0)。S_θθ・divU の下流利用 (SST 生産・τθθ) は不変。
- **ソース** (`axisymmetricSource_d.cu` に新カーネル `axisymmetricSourceSU2_d`):
  §2 の非粘性+粘性を forge 符号 (res += 物理ソース) で実装。AuxVar 3 種を毎ステップ計算し
  専用 GG 勾配 (planar, 境界は owner 値) で ∇A を得る。軸ノード (`axis_flag`) はソース 0。
- **SST** (`rans_sst_source_d` 内, method 1 ゲート): §2 の k/ω ソース + `src_jac_*` へ
  max(v/y, 0) を加算 (対角正値性維持)。F1 ブレンド σ_k/σ_ω は既存計算を利用。
- **陰解法** (`implicit_defect_correction_block_d`): method 0 の hoop 対角はそのまま、
  method 1 では SU2 4×4 Jacobian を 5×5 の行列 (行/列 0,1,2,4) に加算。軸ノードはスキップ
  (`axis_flag_src` を新引数で渡す — 全行 decouple 用の `axis_flag` とは別ゲート)。
- **不要になるもの (method 1)**: `nodeAxisDirichlet` (軸 CV は有限 planar 体積で普通に解ける)、
  `axisTimestepBeta`、hoop ソース (`axisymmetricSource_d`)。`zeroAxisRadialResidual` は
  SU2 の対称面射影に相当するため**残す**。

## 4. 検証

1. **ビット不変**: `axisymMethod: 0` (既定) は全経路不変 (cell/node とも)。
2. **node SST 生産構成** (case/40, 2 次+陰解法 cfl4, `nodeAxisDirichlet` **なし**):
   軸行が健全 (床ピン 0・真空化なし) で 12000 step 完走、η_CF が method 0 (0.9896) と整合。
3. **laminar node** smoke。cell + method 1 は当面スコープ外 (cell は method 0 で健全)。
4. 併せて壁 T 市松 (別欠陥, boundary plan §2.8) が method で変わらないことを記録。

## 5. 変更ログ

- `2026-08-04` — 起票 (ユーザ指示: 軸処理は SU2 を全面的に踏襲)。SU2 コード精読・移植設計。
- `2026-08-04` — **実装・一次検証完了**。§3 の全要素を実装 (config/幾何スキップ/
  `axisymmetricSourceSU2_d` 非粘性+粘性+AuxVar GG 勾配/SST 1/y ソース/block-DPLUR に
  SU2 4×4 Jacobian + 粘性 stiff 主対角 2μ/(ρy²))。full rebuild。
  - **検証 (case/40 node, `run_0026_node_su2axis`)**: 軸行は**真に健全** (床ピン 0、
    k シート消滅 (row1 k max 1.35)、軸線 T 滑らか)。`nodeAxisDirichlet`/`bndFirstOrder`
    なしの全域 2 次+陰解法 cfl4 で 12000 step 完走、η_CF=0.9904 (method 0 の 0.9896 と
    0.08% 差)・quasisteady ALL STEADY。
  - **注意 (warm start)**: method 0 世代の warm 場は軸行が真空のままなので、method 1 へ
    引き継ぐ際は**軸が健全な場** (nodeAxisDirichlet 済み run 等) から warm すること
    (真空軸 IC からは初手発散)。
  - **残課題 = 陰解法の収束深さ**: implicit は喉部近軸のリミットサイクルで rms_ro ~6-8e-5
    プラトー (cfl 2/4 とも、implicitRelaxSST・粘性主対角追加でも不変)。**explicit は
    rms_ro 3.0e-7 に到達** = 空間スキームは健全で、stiff 1/y ソースの陰的線形化の不足が
    原因。1/y ソースの近傍結合 (off-diagonal) を含む線形化、または近軸限定の
    under-relaxation が次の打ち手。
  - **非退行**: method 0 (既定) は新旧バイナリ差が過渡ノイズ増幅の床内 (同一バイナリ再走
    17.5% vs 新旧 11%) で不変。
  - **当面の運用**: 生産既定 (runner) は method 0 + `nodeAxisDirichlet` を維持
    (PASS 品質の収束)。method 1 は `mesh.axisymMethod: 1` のオプトイン。陰解法の
    線形化改善後に既定切替を再評価。
- `2026-08-04` — **陰解法プラトーの Jacobian 作戦を切り分け** (env ゲート診断
  `FORGE_DIAG_SU2JAC_OFF` / `SU2VISC_OFF` / `SU2SST_OFF` を追加, 同一 warm 場から各 12000 step):
  - T1 `implicitSolvePrecision: 1` (5×5 solve double 化): **不変** (8.3e-5) → 条件数/精度説は棄却。
  - T2 非粘性ソース Jacobian OFF (完全 lag): rms_ro **改善** (5.9e-5) だが rms_roK/roOmega が
    rising に転落 → **解析 Jacobian は SST 安定に寄与する一方、平均流サイクルを僅かに悪化**させる。
  - T3 粘性軸対称ソース OFF: 不変 (8.4e-5) → 粘性ソース lag は無罪。
  - T4 SST 軸対称ソース OFF: 不変 (8.3e-5) → SST ソースも無罪。
  - ソース CFL の見積り: dt·|v|/y ~ 0.1 @cfl4 で形式的に stiff でない。以上より**対角 Jacobian の
    工夫 (完全化・精度・選択的無効化) だけではプラトーは解けない**と結論。残る打ち手は
    ①ソース込み完全線形化のオフ対角 (DPLUR sweep の近傍演算子にソース結合を追加)、
    ②ステージ内セミ陰的ソース更新 (Rosenbrock 風)、いずれも DPLUR の構造拡張を要する中規模作業。
    explicit が 3e-7 に達する事実から定常解は健全で、必要なら「implicit で場を作り explicit で
    磨く」二段運用も実用解。
- `2026-08-04` — **対案 `axisRFloor` (ユーザ提案: 規定 r 以下は r 重み・ソース・Jacobian とも
  課さない) を実装・検証** (methods/axisymmetric/implementation.md の節参照)。素朴な面床は
  混在 CV の圧力閉性を壊し step 19 で発散 → **離散閉性面積 (A_closure_x/y) への置換**で解決。
  case/40 run_0027 + soak: 軸健全・η=0.9906・rms_ro ~1e-7 プラトー・床帯縁の k 帯 (~230) は
  24000 step でビット不変の有界平衡。3 方式の比較は case/40 README を参照。

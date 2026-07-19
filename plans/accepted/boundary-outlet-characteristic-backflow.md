# 静圧固定流出 BC の特性ベース化と逆流統一

## メタ

- **area**: `boundary`
- **status**: `done`
- **related_docs**:
  - `methods/boundary.md`
  - `methods/boundary.md`
- **related_plans**:
  - `.github/plans/discretization-node-boundary-ghostless.md` (node 境界の弱形式化、壁∩出口コーナー退化幾何)
  - `.github/plans/diffusion-node-wall-viscous-distance.md` (AddTauWall, node 壁関数)
- **created**: `2026-06-25`
- **owner**: `CFD Dev`

## 1. 目的

亜音速静圧固定流出 (`outlet_statPress`) を SU2 `CEulerSolver::BC_Outlet` 同型の特性ベース構成に統一し、
特に**局所逆流 (出口で $U_n<0$) を静圧アンカーのまま扱う**ことで、剥離 BL が出口に達するケース
(擬似衝撃ダクト, 壁∩出口コーナー) の発散を解消する。

## 2. スコープ

- **やる**:
  - forward (流出) を `P=P_exit` 規定 + 内部エントロピー + 外向き Riemann で構成 (自己整合)。
  - backflow ($U_n<0$) を同じ静圧アンカーで統一 ($V_n^{\text{exit}}<0$ 許容、クランプ無)。
  - 旧 backflow 分岐 (全圧 $P_t,T_t$ stagnation 流入 + 歪速度 $-U_x[ic]\,n_x$) を撤去。
- **やらない**:
  - 乱流スカラー $k,\omega$ の逆流処理変更 (現状ゼロ勾配=SU2 整合・無害、別途精度改善の余地のみ)。
  - case36 SST の残課題 (node 171mm vs cell 142/SU2 132 の差・残差プラトー) — 出口とは無関係と判明、別 plan。
  - TP (`thermalMethod==2`) の特性構成厳密化 (現状 CPG 同型近似)。

## 3. 関連 docs と前提

- `methods/boundary.md` / `implementation.md` の `outlet_statPress` 節。
- SU2 クロスチェック: `CEulerSolver::BC_Outlet` は常に静圧 `P_Exit` をアンカーし、方向分岐せず、
  $V_n$ のクランプもしない (全圧 stagnation 構成は `BC_Inlet` TOTAL_CONDITIONS 専用)。

## 4. 設計

実装: `solver_density_cuda/cuda_forge/boundaryCond_d.cu` `outlet_statPress_d`。

- 超音速流出 ($M_L\ge1$ かつ $U_n>0$): 全量外挿 ($P=P_L,\rho=\rho_L$)。
- それ以外 (亜音速流出 / 逆流): $s=P_L/\rho_L^\gamma$, $R^+=U_n+2c_L/(\gamma-1)$,
  $\rho_R=(P_{\text{exit}}/s)^{1/\gamma}$, $c_R=\sqrt{\gamma P_{\text{exit}}/\rho_R}$,
  $V_n=R^+-2c_R/(\gamma-1)$, $\mathbf u_R=\mathbf u_L+(V_n-U_n)\hat{\mathbf n}$, $P_R=P_{\text{exit}}$。
- 値の取得元: 静圧のみ境界 `Psb[ib]`、エントロピー/Riemann/接線速度は内部 `ic`、面法線は plane。

## 5. 検証

case: `case/36.passive_pseudoshock_control` (擬似衝撃ダクト, 壁∩出口コーナーで逆流)。

| run | 出口処理 | 結果 |
| --- | --- | --- |
| `run_node_lam1st_bp1p70` (main, naive) | naive forward + 全圧 backflow | step7185 で壁∩出口コーナー発散 |
| `run_0071_node_lam1st_bffix_bp1p70` | 特性 forward + 特性 backflow | step9999 まで生存・NaN 無 (発散解消) |
| `run_0067` (naive) vs `run_0072` (bffix) — case36 SST Ps1.90 | — | 衝撃軌跡ビット一致 (171mm/Mmax1.92)、回帰無害 |

- 逆流が効く層流コーナーで発散解消、効かない SST Ps1.90 で完全不変、という理想挙動。
- VERDICT (`check_convergence.py`): run_0071=NOT CONVERGED (still converging, 発散せず生存),
  run_0072=NOT CONVERGED (plateau, run_0067 と同一)。いずれも「収束」とは報告しない。

## 変更ログ

- 2026-06-25: 実装・検証完了。`outlet_statPress_d` を特性ベース forward + 静圧アンカー backflow に統一。
  旧全圧 stagnation backflow 分岐 (剥離域過加圧・歪速度・forward↔backflow ばたつきで発散) を撤去。
  層流コーナー発散解消・SST 回帰無害を確認。診断 run: `run_0068` (forward-only 対照, case36 で inert),
  `run_0073/0074` (node vs cell 層流比較, 別調査)。

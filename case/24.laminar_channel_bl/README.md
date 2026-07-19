# case/24.laminar_channel_bl — 層流チャネル入口境界層

平行平板チャネルの層流発達流れ。小規模 (~13k cells) で回るため、壁粘性・境界処理・離散化
(cell/node) の切り分け用テストベッドとして使われてきた (SU2 クロスチェック、node 境界
ghostless 化の A/B、median-dual 検証、WMLES 壁モデルのスモーク)。

- メッシュ: `mesh/` (構造化四角形、壁 y+ は層流想定)
- 比較スクリプト: `compare_ghostfix.py` (ghost 撤廃 A/B)、`comparison_*.png` / `m1_cell_vs_node_residuals.png`

## 計算 run 一覧

註: `run_0001`〜`run_regr_*` は本 README 整備 (2026-07-20) 以前の run。目的は名称・config・
関連 plan から要約したもので、詳細は各 run 内の config / 関連 plan を参照。

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_laminar` | SLAU 層流ベースライン | — | ref |
| `run_0002_su2_laminar` / `run_0006_su2_laminar_converged` | 同一メッシュ・同一 BC の SU2 クロスチェック ([procedures/su2-cross-check.md](../../procedures/su2-cross-check.md)) | `comparison_outlet.png` 等 | ref |
| `run_0003_slau_laminar_viscfix` / `run_0004_..._long` | 壁粘性 flux 修正の A/B (+長時間版) | `comparison_residuals.png` | ref |
| `run_0005_roe_laminar_long` | ROE 版長時間 | — | ref |
| `run_0006_slau_ghost_baseline` / `run_0007_slau_ghost_fix` | node 境界 ghostless 化の A/B ([discretization-node-boundary-ghostless plan](../../plans/active/discretization-node-boundary-ghostless.md)) | `compare_ghostfix.py` | ref |
| `run_dual_*` (cell/node × m1/m2chk/m3diag) | median-dual M1-M3 の cell vs node 検証 ([discretization-median-dual plan](../../plans/archived/discretization-median-dual.md)) | `m1_cell_vs_node_residuals.png` | ref |
| `run_dplurchk_*` | block-DPLUR restart / explicit 回帰スモーク | — | ref |
| `run_regr_single` / `run_regr_multi` | 回帰確認用スモーク | — | ref |
| `run_wmles_smoke_cell` | **WMLES 壁モデルスモーク (cell)**: 層流チャネルで `wallModelLES:1` (ILES, `LESorRANS:0`) + RK3 cfl1、2000 step。壁モデル経路が cell で回るかの機能確認 ([turbulence-wmles-wall-stress plan](../../plans/active/turbulence-wmles-wall-stress.md)) | **完走・NaN 無し・WMLES 経路稼働** (ログ `[WMLES]` カウンタ)。y+≪1 の層流域のため u_τ Newton は設計どおり laminar fallback (τ_w=μu/d) が主。P/T/Ux 物理的。VERDICT=NOT CONVERGED (rms_ro/roe plateau — スモーク範囲では未収束のまま、機能確認のみ) | active (WMLES スモーク) |
| `run_wmles_smoke_node` | 同 **node 版** (`discretization: node`, implicit)。node 取得層 (Normal_Neighbor) + `Tau_Wall` 再スケール経路の機能確認 | **完走・NaN 無し・WMLES 経路稼働**。VERDICT=NOT CONVERGED (全列 falling = 収束途上で打ち切り)。挙動は cell と同傾向 | active (WMLES スモーク) |

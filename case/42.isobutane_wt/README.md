# case/42.isobutane_wt — イソブタン燃焼ガス風洞 (semi-perfect, axis-Mach チェーン)

案件 (2026-08-17): $P_t$=5 MPa, $T_t$=1000 K, C₄H₁₀/空気 φ=0.9 生成物
(CO₂ 16.8 / H₂O 8.6 / O₂ 2.2 / N₂ 72.4 wt%), 出口径 1 m, 入口配管径 1 m, **$M_d$=4.2** (当初 4.0、同日訂正)。
$r_t$=0.12575 m (semi-perfect の $A/A^*(M4.2)$=15.809 から)。
計画: [plans/active/tooling-nozzle-semiperfect-gas.md](../../plans/active/tooling-nozzle-semiperfect-gas.md)。

**熱力学**: 設計 (MOC) は NASA-9/CEA の semi-perfect (frozen 組成、γ 1.31→1.38)。
CFD は forge の TP (単一擬似種 `species_db.yaml`) を意図したが **TP × node 軸対称が
forge 側で発散** (下記 run_0001) → 当面 **CFD は CPG(γ\*=1.309)** (`evaluate.cfd_gas: cpg`)。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_ib_base` | **TP (semi-perfect) CFD の初回投入** (`problem_ib_m4.yaml`, 5 回試行: interp_field CPG 再構成バグ → 背圧 1000→28000 Pa → cfl4→mid 段+cfl2 → axisymMethod 1 → method 0 + axisRFloor 2.8 mm) | 全て発散。最終形は soft 段 step 271 で膨張部中域 (x≈5–6) の軸近傍から。**同一条件の CPG (run_0002) は完走** → **TP × node 軸対称は forge 側の穴** (nodeAxisDirichlet が TP 不可、代替なし)。別セッション (軸対称高度化) へ申し送り | 破棄予定 (**forge 側問題の記録**) |
| `run_0020_ib_R3_tp_newbin` | **TP (semi-perfect) CFD の再挑戦 — 新 binary** (`ea08bcbe`: nodeAxisDirichlet 撤去、軸ノードが真の DOF)。`problem_ib_R3_LU6_Lcmax_tp.yaml`, cfd_gas 無し・axisRFloor 無し | **再び発散** (soft 段 step 190、x≈7.0–7.5 の軸上で P が床 1.46 Pa まで落下)。IC は健全 (T 247–1000 K)。旧 binary の method 0 + Dirichlet 0 (step 271, x≈5–6 軸) と同じ性質 → **TP × node 軸対称の発散は nodeAxisDirichlet と無関係、forge の TP 経路 (EOS/軸ソース項) 固有**。同一メッシュ・IC の CPG は新旧 binary とも完走。申し送り継続 | 破棄予定 (**forge 側問題の記録・第 2 弾**) |
| `run_0002_cpgcheck` | 同一壁 (semi-perfect 設計) を **CPG(γ\*) で CFD** (`problem_ib_cpgcheck.yaml`) — TP×軸対称の切り分け | **完走**: ‖ΔM‖∞ 0.284% $M_d$・出口 ε_M 0.013%/ε_θ 0.0096°。品質 PASS、軸 M は 4000 step 間 8.8e-6 で凍結 (残差は warm 床プラトー)。CFD 経路の基準 | active (**CPG 経路の基準**) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (初回, 破棄) | 設計 semi-perfect × CFD CPG(γ\*) の**熱力学不整合** run | 全 9 点で ΔM≈4.5%・ε_M 3.5%・出口コア M 3.86 — 壁 (A/A*=13.1) と CFD (CPG 1.309 なら 15.3) の食い違いで M4 に届かない。**設計と CFD の熱力学は必ず一致させる**教訓。同名で再投入 (下) | 破棄済 (記録のみ) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (2 回目) | 整合 CPG だが **M4 と M4.2 が混在** (途中でユーザ訂正 4.0→4.2 を YAML に反映したため 0003–0006 が M4、0007–0011 が M4.2) | 個々の値は健全 (0.30–0.68%) だが系列として使わない | 破棄予定 (混在の記録) |
| **`run_0012`〜`run_0019_ib_R*_LU*_Lc*`** | **R/L スタディ CFD (M4.2 統一, 設計・CFD とも CPG γ\*=1.309, `cfd_gas: cpg`)**: R∈{1.5,2,3}×L_c∈{max,8} @L_U=6、L_U∈{4,9} @R2 Lc=max。CPG では出口径 1.09 m (semi-perfect 正本 1.0 m) — 相対比較用 | **R3/L_U6/L_c=max が最良 (‖ΔM‖∞ 0.256%, ε_θ 0.003°)**。L_c=8 は R≥2 でゲート外 (0.52/0.68%)。L_U は軸 M に効かず ε_M に効く (L_U4 で 0.061% vs L_U9 で 0.010%)。全 run 品質 PASS・NaN 0・軸 M 凍結 ~1e-5。`study_cfd_m42.json`、[plan §3.5](../../plans/accepted/tooling-nozzle-semiperfect-gas.md)、[結果ページ](https://claude.ai/code/artifact/07b52b7f-54b1-4cc8-97a5-4ee9af380b2d) | active (**スタディの根拠**) |
| — | 設計のみスイープ 36 点 → `study_design.json` (R=5 全滅、L_c=6 は R=1.5 のみ、L_U=9 は R=1.5 で μ>20) | | ref |

関連 case/41: `run_0076_cpg_axisymmethod1` = CPG でも `axisymMethod:1` は出口軸コーナーで発散
(TP 発散の切り分け副産物、別セッションへの申し送り)。

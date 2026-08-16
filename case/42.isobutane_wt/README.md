# case/42.isobutane_wt — イソブタン燃焼ガス風洞 (semi-perfect, axis-Mach チェーン)

案件 (2026-08-17): $P_t$=5 MPa, $T_t$=1000 K, C₄H₁₀/空気 φ=0.9 生成物
(CO₂ 16.8 / H₂O 8.6 / O₂ 2.2 / N₂ 72.4 wt%), 出口径 1 m, 入口配管径 1 m, $M_d$=4。
$r_t$=0.13812 m (semi-perfect の $A/A^*(M4)$=13.105 から)。
計画: [plans/active/tooling-nozzle-semiperfect-gas.md](../../plans/active/tooling-nozzle-semiperfect-gas.md)。

**熱力学**: 設計 (MOC) は NASA-9/CEA の semi-perfect (frozen 組成、γ 1.31→1.38)。
CFD は forge の TP (単一擬似種 `species_db.yaml`) を意図したが **TP × node 軸対称が
forge 側で発散** (下記 run_0001) → 当面 **CFD は CPG(γ\*=1.309)** (`evaluate.cfd_gas: cpg`)。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_ib_base` | **TP (semi-perfect) CFD の初回投入** (`problem_ib_m4.yaml`, 5 回試行: interp_field CPG 再構成バグ → 背圧 1000→28000 Pa → cfl4→mid 段+cfl2 → axisymMethod 1 → method 0 + axisRFloor 2.8 mm) | 全て発散。最終形は soft 段 step 271 で膨張部中域 (x≈5–6) の軸近傍から。**同一条件の CPG (run_0002) は完走** → **TP × node 軸対称は forge 側の穴** (nodeAxisDirichlet が TP 不可、代替なし)。別セッション (軸対称高度化) へ申し送り | 破棄予定 (**forge 側問題の記録**) |
| `run_0002_cpgcheck` | 同一壁 (semi-perfect 設計) を **CPG(γ\*) で CFD** (`problem_ib_cpgcheck.yaml`) — TP×軸対称の切り分け | **完走**: ‖ΔM‖∞ 0.284% $M_d$・出口 ε_M 0.013%/ε_θ 0.0096°。品質 PASS、軸 M は 4000 step 間 8.8e-6 で凍結 (残差は warm 床プラトー)。CFD 経路の基準 | active (**CPG 経路の基準**) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` | **R/L スタディの CFD** (設計 semi-perfect, CFD CPG。R∈{1.5,2,3}×L_c∈{max,8,(6)} @L_U=6、L_U∈{4,9} @R2 Lc=max) | (投入中) `study_cfd.json` | active |
| — | 設計のみスイープ 36 点 → `study_design.json` (R=5 全滅、L_c=6 は R=1.5 のみ、L_U=9 は R=1.5 で μ>20) | | ref |

関連 case/41: `run_0076_cpg_axisymmethod1` = CPG でも `axisymMethod:1` は出口軸コーナーで発散
(TP 発散の切り分け副産物、別セッションへの申し送り)。

# case/42.isobutane_wt — イソブタン燃焼ガス風洞 (semi-perfect, axis-Mach チェーン)

案件 (2026-08-17): $P_t$=5 MPa, $T_t$=1000 K, C₄H₁₀/空気 φ=0.9 生成物
(CO₂ 16.8 / H₂O 8.6 / O₂ 2.2 / N₂ 72.4 wt%), 出口径 1 m, 入口配管径 1 m, **$M_d$=4.2** (当初 4.0、同日訂正)。
$r_t$=0.12575 m (semi-perfect の $A/A^*(M4.2)$=15.809 から)。
計画: [plans/active/tooling-nozzle-semiperfect-gas.md](../../plans/active/tooling-nozzle-semiperfect-gas.md)。

**熱力学**: 設計 (MOC) は NASA-9/CEA の semi-perfect (frozen 組成、γ 1.31→1.38)。
CFD は forge の TP (単一擬似種 `species_db.yaml`, **`thermoHrefTemp: 298.15` 必須** — 絶対基準では
陰解法 Jacobian の χ_eos が桁違いになり軸で発散 [run_0021–0025 で特定]) で回る (run_0026)。
R/L スタディ (run_0012–0019) は特定前に CPG(γ\*) で実施した相対比較。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_ib_base` | **TP (semi-perfect) CFD の初回投入** (`problem_ib_m4.yaml`, 5 回試行: interp_field CPG 再構成バグ → 背圧 1000→28000 Pa → cfl4→mid 段+cfl2 → axisymMethod 1 → method 0 + axisRFloor 2.8 mm) | 全て発散。最終形は soft 段 step 271 で膨張部中域 (x≈5–6) の軸近傍から。**同一条件の CPG (run_0002) は完走** → **TP × node 軸対称は forge 側の穴** (nodeAxisDirichlet が TP 不可、代替なし)。別セッション (軸対称高度化) へ申し送り | 破棄予定 (**forge 側問題の記録**) |
| `run_0021_ib_tp_constcp` / `run_0022_ib_tp_lowTfreeze` / `run_0023_ib_tp_explicit` / `run_0024_ib_tp_hoopjacfix` / `run_0025_ib_tp_href` | **TP 発散の切り分け系列** (run_0020 の soft 段設定 [1次 cfl0.5 3000 step] を 1 要素ずつ変える): 0021 = 一定 cp 擬似種 (γ 1.309)、0022 = 実 NASA-9 で 200 K 未満 cp 凍結、0023 = 実 NASA-9 陽解法 (timeIntegration 3)、0024 = 軸ホップ Jacobian を一般 EOS の ∂P/∂Q (χ_eos 項) に修正した binary、0025 = **`thermoHrefTemp: 298.15` (sensible-enthalpy datum) + IC を同 datum** | 0021 完走 / 0022 step 191 発散 / 0023 完走 / 0024 step 190 発散 / **0025 完走**。→ **原因 = 陰解法 Jacobian の χ_eos = c² − κh に絶対基準 h (生成エンタルピー込み, −1.8e6 J/kg) が入り桁違いになる**。低温外挿・TP 経路自体・ホップ Jacobian は無罪 (ホップ修正は正しいので残置)。打ち手 = runner の TP 経路に thermoHrefTemp 298.15 を既定化 + IC を同 datum で構成 | 破棄予定 (**原因特定の記録**) |
| `run_0028_ib_R3_tp_inletdir` / `run_0029_ib_R3_tp_pipe2` | **TP 残差床 (rms_ro 3e-2 プラトー) の切り分け**: 場は凍っているが入口面で P/Pt=1.03 (全圧超過) が振動 = TP `inlet_Pressure` (静圧/速度参照 rf=0.5 緩和) の病理。0028 = `inlet_Pressure_dir` に変更、0029 = 配管長 0.5→2 r_t | **0028: step 9059 で発散** (更に悪い)。**0029: 配管 2 r_t でも P/Pt 1.03 振動・残差床 3.6e-2 不変** (内部 0.114%) → **TP の入口 BC 2 種とも大径低速入口 (M≈0.03) で問題**。CPG 入口 BC (同 case 完走・rms_ro 1e-9) にはない → **forge 側 TP 入口 BC の申し送り**。ノズル内部 (軸 M・出口) は凍っており設計評価には使えるが「収束」ではない | 破棄予定 (**申し送りの根拠**) |
| **`run_0026_ib_R3_tp_full`** / `run_0027_ib_R3_tp_cont` | **TP (semi-perfect 設計 + semi-perfect CFD) フル run** (`problem_ib_R3_LU6_Lcmax_tp.yaml`, thermoHrefTemp 298.15, soft→mid→本段 cfl2 12000; 0027 = 継続 +24000) | **完走。‖ΔM‖∞ 0.116% Md・出口 ε_M 0.014%・ε_θ 0.0036°・overshoot +0.02%** — CPG(γ*) 検証 (0.256%) の半分以下 (熱力学が設計と完全整合)。出口軸 M 4.2010、コア M 4.2005、**出口径 0.999 m** (目標 1.0)。品質 PASS・NaN 0。残差は rms_ro 3e-2 プラトー・軸 M 8k→12k で 2.9e-3 とまだ動く (CPG 系列の 8e-6 より深さ不足) → 0027 で継続確認 | active (**TP CFD の基準**) |
| `run_0020_ib_R3_tp_newbin` | **TP (semi-perfect) CFD の再挑戦 — 新 binary** (`ea08bcbe`: nodeAxisDirichlet 撤去、軸ノードが真の DOF)。`problem_ib_R3_LU6_Lcmax_tp.yaml`, cfd_gas 無し・axisRFloor 無し | **再び発散** (soft 段 step 190、x≈7.0–7.5 の軸上で P が床 1.46 Pa まで落下)。IC は健全 (T 247–1000 K)。旧 binary の method 0 + Dirichlet 0 (step 271, x≈5–6 軸) と同じ性質 → **TP × node 軸対称の発散は nodeAxisDirichlet と無関係、forge の TP 経路 (EOS/軸ソース項) 固有**。同一メッシュ・IC の CPG は新旧 binary とも完走。申し送り継続 | 破棄予定 (**forge 側問題の記録・第 2 弾**) |
| `run_0002_cpgcheck` | 同一壁 (semi-perfect 設計) を **CPG(γ\*) で CFD** (`problem_ib_cpgcheck.yaml`) — TP×軸対称の切り分け | **完走**: ‖ΔM‖∞ 0.284% $M_d$・出口 ε_M 0.013%/ε_θ 0.0096°。品質 PASS、軸 M は 4000 step 間 8.8e-6 で凍結 (残差は warm 床プラトー)。CFD 経路の基準 | active (**CPG 経路の基準**) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (初回, 破棄) | 設計 semi-perfect × CFD CPG(γ\*) の**熱力学不整合** run | 全 9 点で ΔM≈4.5%・ε_M 3.5%・出口コア M 3.86 — 壁 (A/A*=13.1) と CFD (CPG 1.309 なら 15.3) の食い違いで M4 に届かない。**設計と CFD の熱力学は必ず一致させる**教訓。同名で再投入 (下) | 破棄済 (記録のみ) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (2 回目) | 整合 CPG だが **M4 と M4.2 が混在** (途中でユーザ訂正 4.0→4.2 を YAML に反映したため 0003–0006 が M4、0007–0011 が M4.2) | 個々の値は健全 (0.30–0.68%) だが系列として使わない | 破棄予定 (混在の記録) |
| **`run_0012`〜`run_0019_ib_R*_LU*_Lc*`** | **R/L スタディ CFD (M4.2 統一, 設計・CFD とも CPG γ\*=1.309, `cfd_gas: cpg`)**: R∈{1.5,2,3}×L_c∈{max,8} @L_U=6、L_U∈{4,9} @R2 Lc=max。CPG では出口径 1.09 m (semi-perfect 正本 1.0 m) — 相対比較用 | **R3/L_U6/L_c=max が最良 (‖ΔM‖∞ 0.256%, ε_θ 0.003°)**。L_c=8 は R≥2 でゲート外 (0.52/0.68%)。L_U は軸 M に効かず ε_M に効く (L_U4 で 0.061% vs L_U9 で 0.010%)。全 run 品質 PASS・NaN 0・軸 M 凍結 ~1e-5。`study_cfd_m42.json`、[plan §3.5](../../plans/accepted/tooling-nozzle-semiperfect-gas.md)、[結果ページ](https://claude.ai/code/artifact/07b52b7f-54b1-4cc8-97a5-4ee9af380b2d) | active (**スタディの根拠**) |
| — | 設計のみスイープ 36 点 → `study_design.json` (R=5 全滅、L_c=6 は R=1.5 のみ、L_U=9 は R=1.5 で μ>20) | | ref |

関連 case/41: `run_0076_cpg_axisymmethod1` = CPG でも `axisymMethod:1` は出口軸コーナーで発散
(TP 発散の切り分け副産物、別セッションへの申し送り)。

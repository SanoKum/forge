# case/42.isobutane_wt — イソブタン燃焼ガス風洞 (semi-perfect, axis-Mach チェーン)

案件 (2026-08-17): $P_t$=5 MPa, $T_t$=1000 K, C₄H₁₀/空気 φ=0.9 生成物
(CO₂ 16.8 / H₂O 8.6 / O₂ 2.2 / N₂ 72.4 wt%), 出口径 1 m, 入口配管径 1 m, **$M_d$=4.2** (当初 4.0、同日訂正)。
$r_t$=0.12575 m (semi-perfect の $A/A^*(M4.2)$=15.809 から)。
計画: [plans/active/tooling-nozzle-semiperfect-gas.md](../../plans/active/tooling-nozzle-semiperfect-gas.md)。

**熱力学**: 設計 (MOC) は NASA-9/CEA の semi-perfect (frozen 組成、γ 1.31→1.38)。
CFD は forge の TP (単一擬似種 `species_db.yaml`, **`thermoHrefTemp: 298.15` 必須** — 絶対基準では
陰解法 Jacobian の χ_eos が桁違いになり軸で発散 [run_0021–0025 で特定]) で回る (run_0026)。
R/L スタディ (run_0012–0019) は特定前に CPG(γ\*) で実施した相対比較。

**M6 拡張 (2026-08-16)**: 同ガス組成・$P_t$=5 MPa で **$T_t$=1550 K・$M_d$=6・出口/入口 1 m** ($r_t$=53.75 mm、$A/A^*$=86.55)。
単一 quintic 軸 M 則は M6 で全点 REJECT (θ_w>μ_w の逆 MOC fold) → **knot 則 (A6)** で成立。全長 ≥ $L_c$+55 $r_t$ (≈5–7 m) は
終端特性線 ($\mu_d$=9.6°) の物理的必然。計画: [tooling-nozzle-axismach-knot-law.md](../../plans/accepted/tooling-nozzle-axismach-knot-law.md)。

## M4.2 ベスト結果の参照・引き継ぎ (他マシン向け, 2026-08-17)

**ベスト = `run_0035_ib_R2_LU6_Lc8_tp` (R2 / L_U6 / L_c8, ‖ΔM‖∞ 0.084 % $M_d$)**、次点
`run_0036_ib_R3_LU6_Lc8_tp` (R3, 0.091 %, ε_M 0.009 % で出口一様性は最良)。差は入口 BC 床と
同程度で **R2/R3 は同等**と読む。CFD 成果物 (`res_*.h5`, メッシュ h5/msh, log) は git 対象外だが、
この 2 run だけ**入力・設計成果物をリファレンスとして git に残してある**:

| ファイル | 内容 |
| --- | --- |
| `problem_ib_R2_LU6_Lc8_tp.yaml` / `problem_ib_R3_LU6_Lc8_tp.yaml` | 問題定義 (設計変数・ガス・評価設定)。**これだけで全て再生成できる** |
| `run_0035_…/wall_design.csv` | 逆 MOC 壁テーブル [x_m, r_m, θ_rad, M_wall] (m 単位)。CAD/別ソルバへ渡す形状はこれ |
| `run_0035_…/target_axis_M.csv` | 軸 Mach 目標 M(x) |
| `run_0035_…/solverConfig.yaml` / `bcondConfig.yaml` / `probe.yaml` / `species_db.yaml` | forge 入力 (単一擬似種 TP, `thermoHrefTemp: 298.15`, 本段 cfl2) |
| `run_0035_…/prepare_info.json` / `metrics.json` / `achieved_vs_target.csv` | 設計メタ (x_A, x_E, アンカー, メッシュ ni/nj) と CFD 結果指標 |
| `run_0035_…/MESH_QUALITY.txt` / `CONVERGENCE_VERDICT.txt` | 品質・収束の記録 |
| `study_cfd_m42_tp.json`, `run_study_m42_tp.py` | 8 点スタディの表と投入ドライバ |

**再生成 (別 PC)**: forge をビルド (`solver_density_cuda/`, [procedures/development-environment.md](../../procedures/development-environment.md)) し、
`design/` に venv を作ってから
```
cd design && ./.venv-opt/bin/python -m forge_design.evaluate.runner_axismach \
    ../case/42.isobutane_wt/problem_ib_R2_LU6_Lc8_tp.yaml ../case/42.isobutane_wt/run_00NN_ib_R2_LU6_Lc8_tp [--prepare-only]
```
`--prepare-only` で設計 + メッシュ (`nozzle.msh`/`nozzle.h5`) + config だけ作れる (メッシュは TFI 構造化を
`design` 側で直接書くので Gmsh 不要、変換に forge の `convertGmshToForge` を使う)。CFD は 12000 step ≈ 20 s
(RTX 3060)。**現行 binary は TP 入口 BC 修正済み** (`run_0050` で残差床 3e-2→9e-7)。run_0031–0038 は修正前の
binary の結果なので、再実行すると残差は深く落ちるが軸 M/出口指標は ~0.01 % 以内で同等 (run_0033 vs run_0050 実測)。

**凝縮計算へ引き継ぐときの注意**:
- **凝縮は `evaluate.tp_species: split_h2o` で回せる (2026-08-17 実装・検証済)**: H₂O 以外を擬似種 `MIXDRY`、H₂O を
  独立種にした 2 種 TP + `evaluate.condensation: {condensation: 1, nCondSpecies: 1, condModel: 1, condKantrowitz: 1, condGrowthModel: 0}`
  (`condGasSpecies` は自動)。実例 `problem_ib_R2_LU6_Lc8_h2o5_split_cond.yaml` → `run_0065`。dry では pseudo と同一
  (run_0063 vs 0064)。Wyslouzil (case/16) で既存 CPG 凝縮結果を再現済 ([計画](../../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md))。
  多成分 TP × 陰解法は cfl ≤2 (段階起動、`thermoHrefTemp` 298.15) で安定。
- 形状・メッシュ・IC は `problem_ib_R2_LU6_Lc8_tp.yaml` から `--prepare-only` で再生成し、`solverConfig.yaml` の
  `physProp`/`condensation` ブロックだけ差し替える。IC は `paste_isentropic_ic` (設計 1D 等エントロピー) が既定。
- 出口静温 ~247 K・静圧 ~20–28 kPa (M4.2, Tt 1000 K) — H₂O 分圧はここから。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_ib_base` | **TP (semi-perfect) CFD の初回投入** (`problem_ib_m4.yaml`, 5 回試行: interp_field CPG 再構成バグ → 背圧 1000→28000 Pa → cfl4→mid 段+cfl2 → axisymMethod 1 → method 0 + axisRFloor 2.8 mm) | 全て発散。最終形は soft 段 step 271 で膨張部中域 (x≈5–6) の軸近傍から。**同一条件の CPG (run_0002) は完走** → **TP × node 軸対称は forge 側の穴** (nodeAxisDirichlet が TP 不可、代替なし)。別セッション (軸対称高度化) へ申し送り | 破棄予定 (**forge 側問題の記録**) |
| `run_0021_ib_tp_constcp` / `run_0022_ib_tp_lowTfreeze` / `run_0023_ib_tp_explicit` / `run_0024_ib_tp_hoopjacfix` / `run_0025_ib_tp_href` | **TP 発散の切り分け系列** (run_0020 の soft 段設定 [1次 cfl0.5 3000 step] を 1 要素ずつ変える): 0021 = 一定 cp 擬似種 (γ 1.309)、0022 = 実 NASA-9 で 200 K 未満 cp 凍結、0023 = 実 NASA-9 陽解法 (timeIntegration 3)、0024 = 軸ホップ Jacobian を一般 EOS の ∂P/∂Q (χ_eos 項) に修正した binary、0025 = **`thermoHrefTemp: 298.15` (sensible-enthalpy datum) + IC を同 datum** | 0021 完走 / 0022 step 191 発散 / 0023 完走 / 0024 step 190 発散 / **0025 完走**。→ **原因 = 陰解法 Jacobian の χ_eos = c² − κh に絶対基準 h (生成エンタルピー込み, −1.8e6 J/kg) が入り桁違いになる**。低温外挿・TP 経路自体・ホップ Jacobian は無罪 (ホップ修正は正しいので残置)。打ち手 = runner の TP 経路に thermoHrefTemp 298.15 を既定化 + IC を同 datum で構成 | 破棄予定 (**原因特定の記録**) |
| `run_0028_ib_R3_tp_inletdir` / `run_0029_ib_R3_tp_pipe2` | **TP 残差床 (rms_ro 3e-2 プラトー) の切り分け**: 場は凍っているが入口面で P/Pt=1.03 (全圧超過) が振動 = TP `inlet_Pressure` (静圧/速度参照 rf=0.5 緩和) の病理。0028 = `inlet_Pressure_dir` に変更、0029 = 配管長 0.5→2 r_t | **0028: step 9059 で発散** (更に悪い)。**0029: 配管 2 r_t でも P/Pt 1.03 振動・残差床 3.6e-2 不変** (内部 0.114%) → **TP の入口 BC 2 種とも大径低速入口 (M≈0.03) で問題**。CPG 入口 BC (同 case 完走・rms_ro 1e-9) にはない → **forge 側 TP 入口 BC の申し送り**。ノズル内部 (軸 M・出口) は凍っており設計評価には使えるが「収束」ではない | 破棄予定 (**申し送りの根拠**) |
| `run_0030_ib_tp_lowTfreeze_href` | datum 修正済みで **200 K 未満 cp 凍結**の効果を再検証 (run_0027 と同 IC・6000 step) | 発散なし、軸 M 差 3e-3、出口 M 同一 (4.2011/4.2010)。収束場の T min 247 K で外挿域に触れない。run_0022 の「凍結でも発散」は datum 未修正下の結果で凍結の判定ではなかった (訂正) | 破棄予定 (記録) |
| **`run_0031`〜`run_0038_ib_*_tp`** | **R/L スタディ TP 正本** (semi-perfect 設計 + TP CFD, thermoHrefTemp 298.15, `problem_ib_*_tp.yaml`; CPG 代替 run_0012–0019 と同 8 点) | **全点 0.5% ゲート内、L_U≥6 は 0.08–0.15%**: R1.5/2/3 Lc=max 0.108/0.117/0.115%、**Lc=8 が最良級 R2 0.084% / R3 0.091%** (CPG の「Lc=8 は R≥2 でゲート外」は熱力学不整合の産物 → 撤回)、L_U4 は 0.386% と TP でも明確に悪い (出口 M 4.216 過大)、L_U9 0.102%。ε_M 0.009–0.028%、ε_θ 0.004–0.022°。全 run 品質 PASS・NaN 0。残差は TP 入口 BC の床 (rms_ro 3e-2、軸 M 変動 3–4e-3) で **VERDICT は NOT CONVERGED** — ~0.1% 以下の設計点間差は床と同程度で有意でない。`study_cfd_m42_tp.json`、[結果ページ](https://claude.ai/code/artifact/07b52b7f-54b1-4cc8-97a5-4ee9af380b2d) §3b | active (**TP スタディの根拠・推奨の正本**) |
| **`run_0026_ib_R3_tp_full`** / `run_0027_ib_R3_tp_cont` | **TP (semi-perfect 設計 + semi-perfect CFD) フル run** (`problem_ib_R3_LU6_Lcmax_tp.yaml`, thermoHrefTemp 298.15, soft→mid→本段 cfl2 12000; 0027 = 継続 +24000) | **完走。‖ΔM‖∞ 0.116% Md・出口 ε_M 0.014%・ε_θ 0.0036°・overshoot +0.02%** — CPG(γ*) 検証 (0.256%) の半分以下 (熱力学が設計と完全整合)。出口軸 M 4.2010、コア M 4.2005、**出口径 0.999 m** (目標 1.0)。品質 PASS・NaN 0。残差は rms_ro 3e-2 プラトー・軸 M 8k→12k で 2.9e-3 とまだ動く (CPG 系列の 8e-6 より深さ不足) → 0027 で継続確認 | active (**TP CFD の基準**) |
| `run_0020_ib_R3_tp_newbin` | **TP (semi-perfect) CFD の再挑戦 — 新 binary** (`ea08bcbe`: nodeAxisDirichlet 撤去、軸ノードが真の DOF)。`problem_ib_R3_LU6_Lcmax_tp.yaml`, cfd_gas 無し・axisRFloor 無し | **再び発散** (soft 段 step 190、x≈7.0–7.5 の軸上で P が床 1.46 Pa まで落下)。IC は健全 (T 247–1000 K)。旧 binary の method 0 + Dirichlet 0 (step 271, x≈5–6 軸) と同じ性質 → **TP × node 軸対称の発散は nodeAxisDirichlet と無関係、forge の TP 経路 (EOS/軸ソース項) 固有**。同一メッシュ・IC の CPG は新旧 binary とも完走。申し送り継続 | 破棄予定 (**forge 側問題の記録・第 2 弾**) |
| `run_0002_cpgcheck` | 同一壁 (semi-perfect 設計) を **CPG(γ\*) で CFD** (`problem_ib_cpgcheck.yaml`) — TP×軸対称の切り分け | **完走**: ‖ΔM‖∞ 0.284% $M_d$・出口 ε_M 0.013%/ε_θ 0.0096°。品質 PASS、軸 M は 4000 step 間 8.8e-6 で凍結 (残差は warm 床プラトー)。CFD 経路の基準 | active (**CPG 経路の基準**) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (初回, 破棄) | 設計 semi-perfect × CFD CPG(γ\*) の**熱力学不整合** run | 全 9 点で ΔM≈4.5%・ε_M 3.5%・出口コア M 3.86 — 壁 (A/A*=13.1) と CFD (CPG 1.309 なら 15.3) の食い違いで M4 に届かない。**設計と CFD の熱力学は必ず一致させる**教訓。同名で再投入 (下) | 破棄済 (記録のみ) |
| `run_0003`〜`run_0011_ib_R*_LU*_Lc*` (2 回目) | 整合 CPG だが **M4 と M4.2 が混在** (途中でユーザ訂正 4.0→4.2 を YAML に反映したため 0003–0006 が M4、0007–0011 が M4.2) | 個々の値は健全 (0.30–0.68%) だが系列として使わない | 破棄予定 (混在の記録) |
| **`run_0012`〜`run_0019_ib_R*_LU*_Lc*`** | **R/L スタディ CFD (M4.2 統一, 設計・CFD とも CPG γ\*=1.309, `cfd_gas: cpg`)**: R∈{1.5,2,3}×L_c∈{max,8} @L_U=6、L_U∈{4,9} @R2 Lc=max。CPG では出口径 1.09 m (semi-perfect 正本 1.0 m) — 相対比較用 | **R3/L_U6/L_c=max が最良 (‖ΔM‖∞ 0.256%, ε_θ 0.003°)**。L_c=8 は R≥2 でゲート外 (0.52/0.68%)。L_U は軸 M に効かず ε_M に効く (L_U4 で 0.061% vs L_U9 で 0.010%)。全 run 品質 PASS・NaN 0・軸 M 凍結 ~1e-5。`study_cfd_m42.json`、[plan §3.5](../../plans/accepted/tooling-nozzle-semiperfect-gas.md)、[結果ページ](https://claude.ai/code/artifact/07b52b7f-54b1-4cc8-97a5-4ee9af380b2d) | active (**スタディの根拠**) |
| **`run_0051_ib_m6_R3_Lc45_MK2.5_tp_inletfix`** / **`run_0061_ib_m6_R3_Lc50_MK2.5_tp`** / `run_0052`〜`run_0060_ib_m6_*_tp` | **M6 スタディ (Tt 1550 K, Pt 5 MPa, 出口/入口 1 m, r_t 53.75 mm; `problem_ib_m6*.yaml`; knot 軸 M 則 A6 + `axis_dx0` 0.03 + `inlet_Pressure` TP 修正 binary)**: R3 × L_c∈{30,40,45,50,60} @M_K 2.5; R∈{2,5} @L_c50; M_K∈{2,3} @R3/L_c50; L_U∈{8,16} @R3/L_c50 (基底 12)。0051 のみ 36000 step (定常性確認)、他 12000 | **全 11 点 0.03–0.08% $M_d$ (ゲート 0.5% の 1/6 以下)**、ε_M 0.002–0.017%、ε_θ 0.002–0.007°、M_exit 6.001–6.004、overshoot ≤0.07%。傾向: L_c 30→60 で θ_max 19.7→11.2°・dM 0.061→0.035%・全長 5.3→6.9 m; R2 が最良 (0.029%)、R5 も可 (0.040%); M_K 2/2.5/3 は同等 (0.040%、M_K2 は ε_M 最小 0.002%); **L_U8 は悪化 (0.078%, ε_M 0.017%)**、L_U16 は L_U12 と同等 (0.049% だが軸 M がまだ 1.8e-3 動く)。全 run 品質 PASS・NaN 0。残差 rms_ro 1e-6〜3e-5 (soft/mid 段で落ち切った warm 床、`check_convergence` は NOT CONVERGED 判定)、設計区間の軸 M は 8k→12k で 2–5e-5 に凍結 (0051 は 36k まで 3e-5)。`study_cfd_m6.json` / `study_design_m6.json` (設計 36 点全合格)、[knot plan](../../plans/accepted/tooling-nozzle-axismach-knot-law.md) | active (**M6 スタディの正本**) |
| `run_0039_ib_m6_R3_Lc45_MK2.5_tp` / `run_0040_..._cont` / **`run_0050_ib_R3_LU6_Lcmax_tp_inletfix`** | **TP `inlet_Pressure` BC の発散切り分けと修正**: 0039 = M6 初回 (旧 binary, 12000 step で 0.041% だが設計区間の軸 M が 2e-3 動く)、0040 = 継続 +24000 (cfl2) で**入口 P/Pt が 1.04→1.14→1.87→5.1→6.2 と指数成長し M_exit 6.21 まで場が壊れる** = TP 入口 BC (静圧参照ブレンド rf 0.5) の不安定が M_in≈0.011 で顕在化。**修正 (rf→0, CPG 分岐と同じ速度参照)** 後の M4.2 再実行 0050: 残差床 3e-2 → **9e-7**、入口 P/Pt 0.9986–0.9993、‖ΔM‖∞ 0.138%・ε_M 0.025% (0033 と同等) | 0039/0040 破棄予定 (**発散の記録**)、0050 active (**修正の根拠**) |
| **`run_0062_ib_m6_R3_Lc45_bsplineM_tp`** | **軸則比較 B (monotone B-spline) の CFD** (`problem_ib_m6_R3_Lc45_bsplineM.yaml`, R3/L_c45/M_d6, A [`run_0051`] と同一条件で `axis_law: bspline_M` のみ差し替え; [軸則比較 plan](../../plans/accepted/tooling-nozzle-axislaw-smoothness.md)) | **完走 (12000 step, PASS/NaN 0)**。‖ΔM‖∞ **0.183%** (A 12k 相当 0.042% の ~4.4 倍)、ε_M 0.023%、ε_θ 0.013°、M_exit 6.0035。B は $M'''$ ジャンプを消す (A の 0.30→B の 0.06) が、その代償として θ_max が 23.5° (A 14.1°) まで増え、`min(μ_w−θ_w)` 余裕が 0.09° (A 1.78°) とほぼ限界まで縮む — **軸則の見た目の滑らかさは壁品質の改善に直結しない**、が結論。C (非負 dν/dx B-spline) は前倒れ対策後も内部特性線網に恒常的な向き反転 (322–496 箇所、解像度で悪化) と壁 spline リンギング (~1°) が残り hard gate 不合格 (詳細は比較 plan)。詳細比較図・表: `compare_axislaw_ABC.json`/`.py`、[結果ページ](https://claude.ai/code/artifact/cf8a4c74-8d1f-44f4-a256-652cc122d00d) | active (**軸則比較 CFD の根拠**) |
| **`run_0063_ib_R2_LU6_Lc8_h2o5_pseudo`** / **`run_0064_…_split`** / **`run_0065_…_split_cond`** | **M4.2 ベスト形状の H₂O 5 % 版 + 凝縮引き継ぎ** (`problem_ib_R2_LU6_Lc8_h2o5_{pseudo,split,split_cond}.yaml`, 組成 N2 0.7528/CO2 0.1743/O2 0.0229/H2O 0.0501; node Euler TP, 段階起動 12000 step; [plan](../../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md)) — 0063 = 単一擬似種 MIX、0064 = `tp_species: split_h2o` (`[MIXDRY,H2O]`, dry)、0065 = split + H₂O 凝縮 Kw+HK (`condGasSpecies 1`) | **pseudo と split dry は同一** (‖ΔM‖∞ 0.062 %/ε_M 0.0117 %/ε_θ 0.0194° が完全一致、軸 M 差 max 1.7e-5)。**凝縮 ON**: onset x≈6.7 r_t (T≈250 K)、軸 g_max **0.0126** (H₂O の 25 %)、潜熱で出口軸 T 245→278 K・P +4.5 %・**出口軸 M 4.20→3.94** (‖ΔM‖∞ 2.2 %、ε_M 7 % — 凝縮域が半径方向に不均一)。全 run 品質 PASS・NaN 0。残差は warm 床 (rms_ro ~1e-6, NOT CONVERGED 判定; 0065 は 2.7 桁降下)。比較図 `compare_h2o5_split.py` → `compare_h2o5_split.png`/`.json` | active (**凝縮引き継ぎの正本**) |
| — | **軸則 D 案 (内部補間点 1 点 C⁴) の設計比較** (CFD なし): `compare_axislaw_D.py` → `compare_axislaw_D.json` / `compare_axislaw_D_summary.csv`。L_c 60→30 continuation、n_axis 600/1200/2400。D 不採用 ([plan](../../plans/accepted/tooling-nozzle-axislaw-onepoint.md), [結果ページ](https://claude.ai/code/artifact/b3142a8a-0c48-47d5-b7bc-f4653cc939c0)) | | ref |
| — | **A knot の $M_K$ 詳細感度** (CFDなし): `sweep_axislaw_A_MK.py` → `.json` / `_summary.csv` / `.png`。R3、L_c 30–60、M_K 1.5–4.0/0.1、n_axis 600 + 代表1200/2400。45は2.6、50は2.7で壁角・marginが小改善するが topology margin は6/13%低下。既定2.5を維持 ([knot plan](../../plans/accepted/tooling-nozzle-axismach-knot-law.md)) | | ref |
| — | **A knot の最短ロバスト軸長探索** (CFDなし): `optimize_axislaw_A_shortest.py` → `.json` / `_summary.csv` / [`report_axislaw_A_shortest.html`](report_axislaw_A_shortest.html)。hard gate・単峰・μ−θ≥1°・最小sin角≥0.02を制約に$x_F$最小化。600点で35–40を粗→0.01刻み、境界のみ1200点。**最短合格 $L_c=39.05,M_K=2.76$** ($x_F=94.34865r_t$, θ_max 16.1281°、margin 1.00095°、flip/交差0)。図: [`axis Mach`](best_axislaw_A_shortest_axis_mach.png) / [`Mach contour`](best_axislaw_A_shortest_mach_contour.png)、図示場 `best_axislaw_A_shortest_moc_field.npz`。2400点・CFDなし ([plan](../../plans/accepted/tooling-nozzle-shortest-robust-axis.md)) | | ref |
| — | 設計のみスイープ 36 点 → `study_design.json` (R=5 全滅、L_c=6 は R=1.5 のみ、L_U=9 は R=1.5 で μ>20) | | ref |

関連 case/41: `run_0076_cpg_axisymmethod1` = CPG でも `axisymMethod:1` は出口軸コーナーで発散
(TP 発散の切り分け副産物、別セッションへの申し送り)。

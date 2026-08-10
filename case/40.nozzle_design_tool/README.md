# case/40.nozzle_design_tool — ノズル設計ツール (forge_design) の検証ケース

`design/forge_design` (ノズル設計ツール — [親計画](../../plans/active/tooling-nozzle-design-tool.md) /
[Phase 0 子 plan](../../plans/active/tooling-nozzle-phase0-foundation.md)) の評価パイプライン検証用ケース。
問題定義 YAML (`problem_bell_smoke.yaml`) から「区分構成壁 → 構造化 quad msh4.1 →
`convertGmshToForge` → 品質ゲート → 準 1D 等エントロピー IC → forge (SLAU+SST+陰解法) →
収束判定 → 推力メトリクス」を 1 コマンドで実行する:

```bash
PYTHONPATH=design python3 -m forge_design.evaluate.runner \
    case/40.nozzle_design_tool/problem_bell_smoke.yaml \
    case/40.nozzle_design_tool/run_NNNN_<slug> [--steps N] [--prepare-only]
```

条件はいずれも ③ベル (Pt 4 MPa / Tt 1500 K / 背圧 20 kPa / rt 10 mm / ε 9 /
θa 20° / L/rt 7)、メッシュ 221×65 (14,080 cells)、品質 VERDICT PASS (AR max 19 / skew max 0.42)。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_bell_smoke` | E2E 初回スモーク。SST の roK/roOmega を IC 未設定 (既定 0) | step 2–6 で rms_roK/roOmega に NaN → EOS 床洗浄後プラトー (既知指紋の再現)。パイプライン自体は完走し η=0.981 | 破棄予定 (NaN 指紋の記録) |
| `run_0002_bell_smoke_sstic` | IC に roK/roOmega (k=1, ω=18000) を追加、壁第一セル y+~70 (壁関数) | **NaN なし**。残差 1.1–1.7 桁低下後プラトー (`CONVERGENCE_VERDICT.txt`: NOT CONVERGED/stalled — cell モード atomicAdd ノイズ床水準)。**推力は 4000→12000 step で 0.002% ドリフト** (F=1961 N, η=0.9790, ṁ=1.299 kg/s, `metrics.json`)。`residual_history.png` あり | active (Phase 0 E2E 基準) |
| `run_0003_bell_L6` | dv 応答確認: L/rt=7→6 (`problem_bell_L6.yaml`) | η=0.9708, F=1944.7 N (短縮 → 発散損失増で η 低下 ✓)。品質 PASS・NaN なし・プラトー同傾向 | active (dv 感度の記録) |
| `run_0004_bell_L9` | dv 応答確認: L/rt=7→9 (`problem_bell_L9.yaml`) | η=0.9837, F=1970.5 N (延長 → η 単調増 ✓)。ṁ は 3 run で 0.03% 一定 (スロート同一の整合) | active (dv 感度の記録) |
| `run_0005_bell_smoke_node` | node 再計算の初回 (冷間 IC + implicit + 細壁 1e-3) | step 2 から全列 NaN。後の切り分けで真因 = 壁第一セル過細 + レシピ不一致と確定 → [調査 plan](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) §2.5 | 破棄予定 (発散の記録) |
| `run_0006_bell_node_warm` | **node 初収束**: 実証レシピ (壁 5e-3・warm start from run_0002・explicit RK3 cfl0.1・1次・katoLaunder) | **VERDICT: PASS / ALL PASS** (全列 3.0–4.6 桁低下)。η=0.9809, F=1964.8 N (cell run_0002: η=0.9790 と 0.2% 差) | active (node 基準) |
| `run_0007_bell_node_2nd_imp` → `_2nd_expl` | 2 次化・陰解法化の試行 (run_0006 収束場から) | **いずれも step 4–10 で発散** — node×SST の陰解法と 2 次精度は未解決のソルバ課題 (case/29 の実績も 1 次 explicit のみ)。plan §2.5 参照 | 破棄予定 (発散の記録) |
| `run_0008_bell_node_runner` | runner 実装後の E2E 一気通貫確認 (`--warm-from`, 24000 step) | **VERDICT: PASS / ALL PASS** (3.4–5.1 桁・全列 falling)。η=0.9817, F=1966.4 N — node レシピが runner 既定で完結 | active (旧 node 1次 explicit 基準) |
| `run_0009_node_sst_imp_repro` | 課題1再現試行: 1次+陰解法 cfl_pseudo 2 (bndFirstOrder あり), 300 step | 「step 4 発散」は再現せず (300 step 健全)。旧発散記録は 2 次との交絡と判明 | 破棄予定 (切り分け記録) |
| `run_0010_node_2nd_repro` | 課題2再現試行: 2次+explicit cfl0.15 (bndFirstOrder あり), 300 step | 300 step 健全。反実仮想 (bndFirstOrder 除去, scratchpad) は step 8–10 で発散 = **課題2の因子は bndFirstOrder 欠落と確定** | 破棄予定 (切り分け記録) |
| `run_0011_node_2nd_imp` | 2次+陰解法 cfl2 の短尺確認 (300 step) | 300 step 健全 (全列 1.4–2.0 桁低下) — ただし長尺は下記 0012 で発散 | 破棄予定 (切り分け記録) |
| `run_0012_node_2nd_imp_full` | 2次+陰解法 cfl2 の 12000 step (軸修正**前**) | **step 10612 で発散** (roK がベル部近軸で e-fold ~2000 step の指数成長 → roOmega inf)。課題1の実態 = 遅発性 | 破棄予定 (発散の記録) |
| `run_0013_node_imp_long` | 1次+陰解法 cfl2 の 12000 step (軸修正**前**) | **step 7875 で発散** (同モード)。= 課題1は convMethod 非依存、軸行真空化×SST k シートが真因 ([plan §2.6](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | 破棄予定 (発散の記録) |
| `run_0014_node_2nd_expl_long` | 2次+explicit cfl0.15 の 12000 step (bndFirstOrder あり) | **完走・全列 2.5–2.8 桁 falling** = 課題2は bndFirstOrder で解決済みの確証 | active (課題2解消の証拠) |
| `run_0015_node_imp_axisdir` | **nodeAxisDirichlet=1** + 1次+陰解法 cfl2, 12000 step | 完走。rms_ro 1.4e-9 / roK 2.5e-7 (floor)、roK 成長モード消滅。軸床ピン 0 ノード。η=0.9816 (1次 explicit 0.9817 と一致) | active (軸修正の効き) |
| `run_0016_node_2nd_imp_axisdir` | **最終目標: 2次+陰解法 cfl2 + 軸 Dirichlet** (runner E2E, 12000 step) | **VERDICT: ALL PASS**・quasisteady ALL STEADY。η=0.9907, F=1984.4 N, ṁ=1.2923。**12000 step ≈ 20 秒** (explicit レシピ ~7 分の ~20 倍) | active (**node 生産基準**) |
| `run_0017_cell_kato` | cell 対照: run_0002 + katoLaunder (warm from run_0002) | η=0.9790 = run_0002 と同一 → **node−cell の η 差 +1.2% は katoLaunder では説明されない** (2次化の node/cell 離散差) | active (η 帰属の対照) |
| `run_0018_bell_L6_node` | L スイープ node 取直し: L/rt=6 (壁 5e-3 修正版 yaml, warm from run_0016) | 完走 (4.0–6.1 桁低下, roOmega は床で振動 2.9e-3↔5.6e-3)。**η=0.9835** | active (dv 感度 node) |
| `run_0019_bell_L9_node` | 同 L/rt=9。cross-geometry warm は cfl4 直投入で発散 → **段階起動** (soft explicit 3000 step → 陰解法 12000 step) | **VERDICT: PASS**。**η=0.9957** | active (dv 感度 node) |
| `run_0020_node_runner_cfl4` | 2次+陰解法 **cfl_pseudo 4** + bndFirstOrder の E2E (旧既定) | **VERDICT: PASS / ALL PASS** (3.2–4.6 桁)。η=0.9907 (cfl2 の run_0016 と一致) | active (bndFirstOrder あり世代の記録) |
| `run_0021_node_2nd_imp_nobfo` | **bndFirstOrder 撤去の主検証**: 全域 2次+陰解法 cfl4+軸 Dirichlet | 12000 step NaN なし。rms_ro ~1e-7 プラトー (入口/収縮部の微小リミットサイクル, ro 相対 ~1e-5)・SST 4.1–4.2 桁。quasisteady **ALL STEADY**。η=0.9896 | active (全域 2 次の主検証) |
| `run_0022_node_runner_nobfo` | **現 runner 既定** (全域 2次+陰解法 cfl4, bndFirstOrder なし) の E2E 基準 | run_0021 と同等 (η=0.9896, ALL STEADY, プラトー同性状) | active (**runner 既定基準**) |
| `run_0023_bell_L6_nobfo` | L/rt=6 を現既定で取り直し (warm from run_0021) | **VERDICT: PASS (converged)**。η=0.9822 | active (dv 感度 node 現行) |
| `run_0024_bell_L9_nobfo` | L/rt=9 を現既定で取り直し (段階起動: soft explicit 3000→陰解法 12000) | rms_ro 9.2e-9 (実質床, flat 判定)・ALL STEADY。η=0.9942 | active (dv 感度 node 現行) |
| `run_0025_node_wallsmooth` | 壁列エントロピー市松の種切り分け (warm start 壁 T 平滑化) | res_0 は市松 5K に清浄化されるが **12000 step 後に run_0022 とビット同一の 206K 市松が再生成** = 市松は定常解固有 ([plan §2.8](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | active (壁市松の証拠) |
| `run_0026_node_su2axis` | **SU2 流軸対称 (`axisymMethod: 1`)** の一次検証 (全域2次+陰解法 cfl4, nodeAxisDirichlet **なし**, warm from run_0021 収束場) | **軸行が真に健全** (床0・kシート消滅・軸線滑らか)。12000 step 完走・η=0.9904・ALL STEADY。ただし implicit は喉部近軸 limit cycle で rms_ro ~9e-5 プラトー (explicit は 3e-7 到達 = 空間健全) → [SU2化 plan](../../plans/active/axisymmetric-su2-source-formulation.md) | active (axisymMethod 1 基準) |
| `run_0029_node_adiabfix` | **断熱壁熱流束リーク修正の反実仮想** (現既定構成, 修正のみ, warm from run_0022) | **壁 T 市松 mean 206→43.7K (−79%)・壁温 1631→1195K**。η=0.9884 (−0.12%)。残差は `NOT CONVERGED` プラトーだが、4000/8000/12000 の η・ṁ・壁温は手計算で不変。以後の node 断熱壁 run はこの物理が既定 | active (**断熱修正基準**) |
| `run_0030_cell_5e3` | **cell を現行 5e-3 メッシュで取得** (壁温三者比較用, cold start, `problem_bell_smoke_cell.yaml`) | η=0.9813, ṁ=1.3035。壁温 mean 1185.8K = node (1195.9K) と同水準 → **壁温 −240K は node 固有でなく forge 共通** | active (cell 同一メッシュ基準) |
| `run_0028_su2_sst` | **SU2 v8.5 クロスチェック** (同一メッシュを gmsh で .su2 化, 同一 BC, RANS-SST/ROE/MUSCL/Euler implicit+FGMRES, 手順書テンプレ) | 15000 iter + 固定 CFL10 継続 5000 iter。全残差は未収束プラトー (`rms[RhoE]≈10^1.06`, `rms[ω]≈10^0.30`)。軸健全・壁 T 市松なし、η_CF=0.9573。入口コア乱流量は SU2 `k=0.975/ω=17160` vs forge `1/18000` でほぼ一致する一方、SU2 low-Re と forge automatic 壁処理、Roe/SLAU、SST 生産補正、軸離散が未整合の粗比較 | active (SU2 対照) |
| `run_0030_su2_sst_wallfunc` | run_0028 固定 CFL 最終場から継続し、SU2 `STANDARD_WALL_FUNCTION` だけを追加 | 第一内点が forge に接近 (`U=925 vs 843m/s`, `k=9.93e3 vs 8.54e3`, `μt/μ=8.3 vs 5.6`)。ṁ=1.2934もforge 1.2928に一致、診断 η=0.9725 (forge 0.9884)。それでも壁 T は面積重み **1422K**、forge は1193K。SU2 は Crocco–Busemann 断熱回復温度を明示設定するが forge SST 壁関数には熱閉包がないことを確証。`wall_temperature_compare.png`; 残差は `NOT CONVERGED` の診断run | active (熱壁関数の対照) |
| `run_0031_su2_sst_constk` | run_0028 固定 CFL 最終場から継続し、分子熱伝導率を forge と同じ 0.0257 固定へ変更 | SU2 壁 T は1426→**1441K** (+15K) で、forge の低温差を説明しない。forge の Sutherland μ + 固定 λ はベル第一内点で分子 Pr≈1.50–1.87となり air (`Pr≈0.72`) と非整合。残差は `NOT CONVERGED`; `residual_history.png` | active (熱物性の対照) |
| `run_0027_node_rfloor` | **r 床 (`axisRFloor: 3.0e-4`, ユーザ提案)** の検証 (method 0, nodeAxisDirichlet **なし**, 全域2次+陰解法 cfl4) | 軸健全 (床0)・η=0.9906・rms_ro ~1e-7 プラトー。床帯縁に有界の k 帯 (~230; scratchpad soak +12000 step でビット不変の平衡)。素朴な面床は閉性破れで発散 → 離散閉性面積 (A_closure) で解決 | active (axisRFloor 基準) |
| `run_0032_cell_yp1_lowre` | y+≈1 メッシュ調整 1 発目 (`wall_first_frac 2e-4`, 221×81, AR 191 PASS)。cell wf=0+constPr, warm from run_0030 | bell y+ mean 1.84 (狙い未達)・T_w 1378K・η 0.9780・ALL STEADY | ref (frac 感度点) |
| `run_0033_cell_yp1_lowre` | **y+≈1 確定メッシュ** (`frac 1e-4`, AR 381 PASS)。cell `wallTreatmentSST:0`+`thermCondMethod:1, prandtlLam:0.72` (Prt は既定 0.85), warm from run_0032 | bell y+ mean 0.94・T_w 1368.2K・η 0.9779・ṁ 1.2993・ALL STEADY (残差はプラトー) | active (Prt 感度対照) |
| `run_0034_node_yp1_lowre` | node 同条件 (wf=0, Prt 0.85)。**発散切り分けで 2 バグ修正** (双対 CV 重心/体積の float32 桁落ち→wall_dist ガベージ / interp_field の node 照会点) — [plan §2.9](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)。段階起動 (陰 cfl0.5 1次 3000→陰 cfl4 2次, stageA/ 同梱)。warm は node 場 (run_0029) 必須 | 全列 2.8–3.5 桁低下 (rms_ro 3.4e-9)・ALL STEADY。T_w 1365.9K = cell と 2.3K 差・τ_w ~1% 差・ṁ 1.2959 (5e-3 世代の wf=0 ṁ −4.5% 異常は消滅) | active (Prt 感度対照) |
| `run_0035_su2_yp1_lowre` | **SU2 v8.5 low-Re を同一 y+1 メッシュで** (run_0028 cfg 流用: adaptive 15000 + 固定 CFL10 5000 iter, PRANDTL_LAM 0.72 / PRANDTL_TURB 0.9 既定) | rms[Rho] 10^-5.95 (固定 CFL 相でも全列低下継続)。bell y+ mean 0.93・**T_w 1414.2K**・τ_w forge±1–3%・ṁ 1.3017・η 0.9796。**チェンバで T_w 1576K = Tt+76K の非物理超過あり** (断熱壁は Tt 超え不可; `wall_temperature_compare_yp1.png`) | active (**SU2 y+1 対照**) |
| `run_0036_cell_yp1_prt09` | **cell y+1 正式基準**: run_0033 + `turbulentPrandtl: 0.9` (SU2 既定に整合; Prt 0.85→0.9 で T_w +20.5K, scratchpad D5) | **T_w 1388.9K・η 0.9779・ṁ 1.2993**・ALL STEADY | active (**Step1 cell 基準**) |
| `run_0037_node_yp1_prt09` | **node y+1 正式基準**: run_0034 + `turbulentPrandtl: 0.9` | **T_w 1387.2K (cell と 1.7K 差)**・ṁ 1.2959・ALL STEADY。η は出口積分 0.9905 だが**出口列アーティファクトで過大** — 内部列積分 0.970–0.977 ([plan §2.10](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | active (**Step1 node 基準**) |
| `run_0038_node_5e3_tawwf` | **Step2 熱的閉包の主検証**: node 5e-3 wf=1 + `sstThermalWallFunction: 1` (+constPr/Prt0.9), warm from run_0029 | **bell 壁温 1417.9K = SU2 壁関数 1422K と 4.1K 一致** (OFF 1196K→+222K)。η 0.9884・ṁ 1.2928 = OFF 基準と一致・ALL STEADY。※初版の「状態ピン」設計は暴走 (壁温 1832K) し弱閉包に改訂 — [plan](../../plans/accepted/turbulence-sst-thermal-wall-function.md) §4 | active (**熱的閉包 node 検証**) |
| `run_0039_cell_5e3_tawwf` | 同 cell 5e-3 (ghost 閉包→弱閉包に改訂後) | 壁温 1489.5K = **+70–90K 過大** (cell 代表点=第一セルの T が node/SU2 同高さ比 ~100–160K 高い — cell wf=1 BL 熱監査 follow-up)。初版 ghost Dirichlet は Tt 飽和で却下 | active (熱的閉包 cell の限界記録) |
| `run_0040_node_yp30_tawwf` | **Step3 壁関数系列 node** (`problem_bell_yp30.yaml` frac 1e-2, bell y+ mean 98, AR 3.8 PASS) wf=1+熱的閉包+constPr/Prt0.9, warm from run_0038 | 壁温 1405.2K (SU2 1418.9K と 13.7K 差)・ALL STEADY・全列 2.8–5.5 桁低下。η 出口積分 0.9835 (内部列 0.971)・ṁ 1.2864 | active (**Step3 node**) |
| `run_0041_cell_yp30_tawwf` | 同 cell (`problem_bell_yp30_cell.yaml`, bell y+ mean 71), warm from run_0039 | 壁温 1387.2K (±20K 市松)・η 0.9802・ṁ 1.3029・ALL STEADY。**cell 代表点バイアスは y+~70 では消滅** (5e-3 のバッファ層代表点固有) | active (Step3 cell) |
| `run_0042_su2_yp30_wallfunc` | **SU2 STANDARD_WALL_FUNCTION を同一 y+30 メッシュで** (low-Re 15000 → wf 継続 5000+10000 iter) | 壁温 **1418.9K**・η 0.9673・ṁ 1.2868。+10000 iter で全量不変 (準定常; rms[Rho] −4.5 プラトー, `residual_history.png`)。`wall_temperature_compare_yp30.png` | active (**Step3 SU2 対照**) |
| `run_0043_node_yp1_outletfix` | **node 出口列欠陥の根治検証** (run_0037 と同一入力の A/B): `outlet_statPress` の bvar `Psb` 動的化修正 ([plan §2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | 出口列 sag **19.4%→0.39%** (超音速行 0.40%)。**η 出口積分 0.9754 = 内部列 0.9747–0.9755** (アーティファクト解消)。ṁ 1.2959 列間ばらつき消滅・ALL STEADY。壁温は run_0037 と ≤0.85K 差 (出口隣接除く) | active (**node y+1 正 (η 込み)**) |
| `run_0044_node_yp30_outletfix` | 同修正の y+30 wf 系 A/B (run_0040 と同一入力) | sag 13.6%→0.17%。**η 出口積分 0.9835 (偽) → 0.9687 = 内部列 0.9679–0.9687**、SU2 wf (run_0042: 0.9673) と +0.14%。全列 2.8–5.1 桁低下・ALL STEADY。ṁ 1.2864 不変 | active (**Step3 node 正**) |
| `run_0045_node_yp1_outletfix_cont` | run_0043 の res_12000 から +12000 step 継続 (準定常確認) | 場は run_0043 と不変 (出口 P 最大 6e-4)。残差は低レベルリミットサイクル (rms_ro 1.1e-8, rms_roUx 3.5e-5) で ALL STEADY | active (継続確認) |

**軸処理 3 方式の比較 (全域 2 次+陰解法 cfl4, bndFirstOrder なし)**:

| 方式 | 軸ノード | row1 k 帯 | 収束 (rms_ro) | η_CF |
| --- | --- | --- | --- | --- |
| `nodeAxisDirichlet: 1` (現 runner 既定) | 解かずミラー | 3.4 | ~1e-7 プラトー (run_0022; 1 次 or bndFirstOrder 世代なら 1e-9 PASS = run_0015/0016) | 0.9896 |
| `axisymMethod: 1` (SU2 流) | 真に解く | **1.35** | implicit ~1e-4 プラトー (explicit 3e-7) | 0.9904 |
| `axisRFloor: 3e-4` (r 床+閉性補正) | 解く (環状近似) | 230 (有界, 帯定義: y∈[0.5,1.5]mm × x∈[30,60]mm の max k) | ~1e-7 プラトー | 0.9906 |

※ 収束値は check_convergence の実測 (全域 2 次では 3 方式とも NOT CONVERGED/プラトー判定であり、
計量は quasisteady STEADY で担保)。η は同一 warm 系統・同一メッシュだが cell (run_0002) は壁 1e-3
世代メッシュのため単純比較不可。**node−cell 差 +1.1% は出口列アーティファクトと判明 (2026-08-11,
[plan §2.10](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md))**。
2026-08-05 の断熱壁修正 (run_0029) 以後は node の η が ~0.12% 下がる (漏れ加熱の除去)。


**dv 感度まとめ (η–L トレードオフの応答確認)**:
L/rt = 6 / 7 / 9 → η_CF は cell (2026-08-03 前半, 壁 1e-3 メッシュ) 0.9708 / 0.9790 / 0.9837、
**node 全域 2次+陰解法 (現既定, run_0023 / 0021 / 0024, 壁 5e-3)** 0.9822 / 0.9896 / 0.9942
(bndFirstOrder あり世代 run_0018/0016/0019 は 0.9835 / 0.9907 / 0.9957)。いずれも物理的に
正しい単調応答で、node−cell はほぼ一定の +1.1〜1.2% オフセット (トレンド保存)。オフセットの
帰属: katoLaunder ではない (run_0017)・軸修正でもない (修正前 run_0012 も 0.991)。
**2026-08-11 訂正: この +1.1〜1.3% は node 出口列の系統的不整合による出口積分アーティファクト**
([plan §2.10](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) — 内部列で再積分すると
node η は cell と整合する。node の η 出口積分値は修正まで生産値に使わない。

**node モードの状態 (2026-08-04 更新)**: 残っていた 2 課題は解決した
([plan §2.6–2.7](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md))。
2 次発散・陰解法の遅発性発散は**いずれも node 軸対称の軸行真空化が種** (ベル部で軸 CV が
radial 圧力平衡からデカップルし EOS 床まで過膨張 → 偽せん断 → SST k シート) — ソルバ新機能
**`nodeAxisDirichlet: 1`** (軸ノードを radial 隣接からの対称 Dirichlet に置換) で根治。
当初レシピに入れていた `bndFirstOrder: 1` (境界隣接 1 次化) は軸病変のマスキングであり、
壁境界層の精度を損なうため**撤去** (ユーザ指示 2026-08-04)。**現 runner 既定 =
全域 2 次 (convMethod 1) + 陰解法 (blockDPLUR, cfl_pseudo 4) + nodeAxisDirichlet**。
1 評価 12000 step ≈ 20 秒 (旧 explicit レシピ比 ~20 倍)。残差は形状により入口/収縮部の
微小リミットサイクルで ~1e-7 プラトーになるため (cell 全域 2 次と同格)、**計量は
`check_quasisteady.py` STEADY を確認して使う**。warm start は引き続き必須で、
**同メッシュ node 場からが最良**。異形状/cell 場からの warm は cfl4 直投入で初手発散する
ことがあり、その場合は段階起動 (soft explicit 3000 step → 陰解法、run_0019/0024 の手順)。
run_0001–0004 は壁 1e-3/2e-3 世代のメッシュ (現 problem yaml は全て 5e-3 に統一済み)。

注記: 残差プラトーの深掘り (真の定常収束品質) は Phase 2 (Rao 照合) で扱う。η の妥当性
(0.98 前後はベルノズルとして自然なオーダー) も Phase 2 の照合が正式判定。



**壁温・η_CF の生産値確定 (2026-08-11, 壁温真値 3 段階キャンペーン完了)**:
y+≈1 low-Re 三者基準 (run_0033–0037)・熱的壁関数 `sstThermalWallFunction` (run_0038–0039,
[plan](../../plans/accepted/turbulence-sst-thermal-wall-function.md))・y+≈65–300 壁関数系列
(run_0040–0042) による挟み込みの結論:

| 量 (③ベル Pt4MPa/Tt1500K/ε9) | 生産値 | 根拠 |
| --- | --- | --- |
| **ベル部壁温 T_w (断熱)** | **1400 ± 15 K** (理論 T_aw ≈1400) | y+1 low-Re: forge 1387–1389 / SU2 1414。wf+閉包: node 1405–1418 / SU2 1419–1422。旧 wf 値 1193K は熱閉包欠落、旧 cell 1767K は熱物性誤り |
| **η_CF** | **0.978 ± 0.003** | 最良解像 (y+1) で cell 0.9779 / SU2 0.9796 / node 0.9754 (出口欠陥修正後 run_0043, 出口積分=内部列)。壁関数系列は 0.967 (SU2)〜0.981 (cell) に散り、メッシュ/壁処理依存が残る |
| ṁ | 1.29–1.30 kg/s (±0.5%) | 全系列。y+30 wf 系 (SU2/node) は −1% 側 |

**注意**: ① node 出口列欠陥は **2026-08-11 に根治済み**
([plan §2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md), run_0043–0045) —
node の η_CF は出口積分値をそのまま使ってよい (修正後の正値: y+1 low-Re 0.975 / y+30 wf 0.969)。
修正前 run (run_0042 以前) の node η 出口積分値 (0.9835–0.9905 系) は +1.3% 過大のため使わない。
② 壁温は `wallTreatmentSST: 1` +
`sstThermalWallFunction: 1` + `thermCondMethod: 1, prandtlLam: 0.72` +
`turbulentPrandtl: 0.9` を生産構成とする (node 一次)。③ 格子系列は壁処理レジームが
異なるため Richardson/GCI は非適用 (単調収束系列でない)。

**y+≈1 low-Re 三者基準解 (2026-08-11, Step 1 完了)**: 専用メッシュ 221×81 `wall_first_frac 1e-4`
(AR 381 / skew 0.42 PASS, bell y+ mean ≈0.93) で forge cell / forge node / SU2 の低 Re 三者比較を実施
(`problem_bell_yp1*.yaml`, run_0033–0037)。熱物性は SU2 と整合 (`thermCondMethod: 1, prandtlLam: 0.72`,
`turbulentPrandtl: 0.9` — **Prt 0.85→0.9 で壁温 +20.5K の強感度**を発見し 0.9 を正式条件とする)。

| 量 (bell 部) | forge cell (run_0036) | forge node (run_0037) | SU2 (run_0035) | 理論 T_aw |
| --- | --- | --- | --- | --- |
| T_w 面積平均 [K] | 1388.9 | 1387.2 | 1414.2 | ~1400 |
| T_w @x=30mm [K] | 1388.5 | 1391.1 | 1416.6 | 1406 |
| τ_w @x=30mm [Pa] | 1707 | 1754 | 1756 | — |
| ṁ [kg/s] | 1.2993 | 1.2959 | 1.3017 | — |
| η_CF | 0.9779 | 0.9905→**内部列 0.970–0.977** | 0.9796 | — |

結論: ① cell/node は壁温 1.7K 差・τ_w ~1–3% 差で相互整合し、理論 T_aw 1400K の −1% 以内。
② SU2 は bell で +1% 高めだが**チェンバで Tt+76K の非物理超過**があり厳密な正解ではない —
壁温の真値帯は **~1390–1415K (理論 1400K ± 1%)** とみなす (wf=1 の 1193K に対し熱的閉包欠落
−210K が本質、Step 2 で対処)。③ η_CF の低 Re y+1 値は **cell 0.978 / SU2 0.980 / node 内部列
0.975** に収束 — 従来の壁処理依存幅 0.955〜0.988 の下限 0.955 は「未解像 low-Re」、上限 0.988 は
「node 出口アーティファクト」でありいずれも真値でない。④ 発見・修正したソルバ/ツールバグ 3 件
(双対 CV 幾何 float32 桁落ち・interp_field node 照会点・node 出口列不整合 [根治済み 2026-08-11,
run_0043–0045]) は
[plan §2.9–2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) 参照。
比較図: `run_0035_su2_yp1_lowre/wall_temperature_compare_yp1.png`。

**壁温の三者比較と原因 (2026-08-11 更新)**: SU2 low-Re は解析回復温度に近く (1431K)、forge は
cell/node とも約−240K (1186/1196K)。判別実験 D2 で forge `wallTreatmentSST=0` は同じ SU2 low-Re と
壁温1K差・η 0.25%差へ整列する。ただし両者とも `y+≈19–40` でlow-Reとして未解像なので、これだけで
`wf_pk` を欠陥と断定してはならない。追加対照 `run_0030_su2_sst_wallfunc` では、x=30mm第一内点の
`μt/μ=6.88` (forge wf=1 cellは3.80)・T=1015K (forge 997K) と近いまま、SU2はCrocco–Busemann
熱閉包により壁温1427Kを得た。したがって**壁面出力の約240K差の直接原因はforge SST壁関数に
断熱回復温度の熱閉包がないこと**であり、「BL全体のμt不足だけが原因」という帰属は棄却する。
一方 η 0.955〜0.988 の真値と壁応力側の妥当性は、`y+≈1` low-Re と `y+>30` wall-function の
格子系列で未決。元のD1–D3記録と限界は
[調査メモ](../../notes/sessions/wall-temperature-three-way-analysis.md)を参照。

**レビュー対応の実装 (2026-08-11)**: ① `check_quasisteady.py` の step パースバグ修正
(res_outlet_*/拡張子数字の混入 — 既報 ALL STEADY は修正後ツールで再確認済み・結論不変)。
② **`Taw_diag` 診断出力**を SST 壁関数に追加 (Crocco 型 T_aw = T_rep + Pr^{1/3}·U_t²/(2cp)):
現既定 node 構成で **ベル部 mean 1420.2K = SU2 壁関数の壁温 1422K と 2K 差** — 熱的閉包を
適用すれば SU2 と一致することをフォワード確認。壁状態への適用は保存整合の plan 化後。
③ **`thermCondMethod: 1`** (constant-Pr: k=μ(T)cp/prandtlLam, 既定 0=従来) を追加 —
スモーク (scratchpad D4) で壁温 −25K シフト (SU2 の逆向き +14.7K と符号整合、主因でないことも一致)。


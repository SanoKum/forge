# case/46.sern_design — ⑤ SERN 設計チェーンの検証ケース

計画: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)。
問題定義 YAML (`problem_*.yaml`) → `design/forge_design/evaluate/runner_sern.py` で
逆設計 → 2 バンド構造メッシュ (スリットカウル) → forge 平面 2D → 力係数 (`metrics.json`)。

実行例 (design/ で):

```bash
PYTHONPATH=. .venv-opt/bin/python -m forge_design.evaluate.runner_sern \
    ../case/46.sern_design/problem_smoke_euler.yaml ../case/46.sern_design/run_NNNN_<slug>
```

境界タグ: inlet_nozzle=1 / inlet_ext=2 / outlet=3 / ramp=4 / cowl_in=5 / cowl_out=6 / bottom=7 / top_out=8。
壁は `outputHDFflg: 1` で `res_wall_<id>_<step>.h5` を吐き、力係数はそこから積分する
(`metrics/sern_forces.py`)。規約: 推力 = 壁力の −x、揚力 = +y、モーメント = 頭上げ正。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_smoke_euler_cell` | S4(a) スモーク初回。**TE ノードを重複させたメッシュに幅 0 の隙間 (未タグ境界辺 2 本)** があった版。段階起動 soft 3000 + 本段 6000 (cfl 4) | C_T 0.9660 / C_L 0.1426 / C_M −0.938 (run_0002 と同値; 隙間の影響なし)。残差 1.4 桁プラトー | 破棄予定 (run_0002 で置換) |
| `run_0002_smoke_euler_cell` | S4(a) スモーク: S1 テスト 4 の設計 (M_in 2.5, ランプ 15°, カウル 5°/1H, M_c 3.9, f 0.45, p_ext/p_in 0.05, M∞ 6) を cell Euler (slip) で評価し MOC と照合。修正メッシュ (TE 共有ノード) | **C_T 0.9660 (MOC 0.9666, −0.05 %) / C_L 0.1427 (0.1392) / C_M −0.939 (−0.946)**、力係数は STEADY (500 step 毎の出力で 1e-5 以内)。残差は `NOT CONVERGED` (1.4 桁プラトー、せん断層・カウル衝撃の cell 床)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png`, `residual_history.png`, `MESH_QUALITY.txt` PASS (AR 120, skew 0.29) | active (ref) |
| `run_0003_nasa_Lc20_euler` | S4(b) NASA TM X-71972 傾向照合: 平板ランプ 20°/18.54H, カウル 6°, **カウル長 2.0H**, M∞ 10 (q 71850 Pa), γ 1.3, `geometry.mode: straight` | C_T 0.9605 / C_L +0.035 / C_M −0.601 (内面のみ C_T 0.9640 vs MOC 0.9628)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0004_nasa_Lc312_euler` | 同上、カウル長 **3.12H (NASA 基準形)** | C_T 0.9704 / C_L −0.000 / C_M −0.685 (内面のみ 0.9760 vs MOC 0.9754)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0005_nasa_Lc45_euler` | 同上、カウル長 **4.5H** | C_T 0.9720 / C_L −0.014 / C_M −0.641 (内面のみ 0.9800 vs MOC 0.9795)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0006_nasa_cowl3deg_euler` | 同上、カウル長 3.12H、**カウル角 3°** | C_T 0.9786 / C_L −0.031 / C_M −0.904 (内面のみ 0.9797 vs MOC 0.9791)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |
| `run_0007_nasa_cowl12deg_euler` | 同上、カウル長 3.12H、**カウル角 12°** | C_T 0.9314 / C_L +0.120 / C_M −0.287 (内面のみ 0.9662 vs MOC 0.9656; カウル外面の衝撃圧が推力 −0.035)。力係数 STEADY、残差プラトー (`NOT CONVERGED`)。`metrics.json`, `mach_field.png`, `wall_pressure_cfd_vs_moc.png` | active (ref) |

| `run_0008_smoke_sst_node` | RANS 化 (S3): node + SST (壁関数), 本段 cfl 4 | **発散** (本段 step 3, ω が inf) | 破棄予定 |
| `run_0009_smoke_sst_node_cfl1` | 同上、本段 cfl 1 | **発散** (step 7, カウル角部直下流 x≈0.13H の内面側で roOmega inf。node+SST 固有、未解決) | 破棄予定 |
| `run_0011_smoke_sst_cell_cfl2` | RANS 化: **cell + SST (wallTreatmentSST 1, 第一セル 4e-3 H, y+~20–40)**, 本段 cfl 2, 8000 step | **C_T(p) 0.9685 / C_T(p+τ) 0.9602 (摩擦 −0.0085) / C_L 0.154 / C_M −0.980**、力係数 STEADY。残差 rms_roK/roOmega は増加後プラトー (`NOT CONVERGED`)。ランプ壁圧に 3% の偶奇鋸歯。`deltastar_vs_euler.png`, `wall_offset_deltastar_vs_euler.json` | active (ref) |
| `run_0012_smoke_sst_cell_dstar` | S5 δ* 補正 (単独抽出版 δ*: 膨張扇の ρu 変化を欠損と誤認し δ*/H 0.24 の過大値) | C_T 0.958 — δ* が誤りで無効 | 破棄予定 |
| `run_0013_smoke_sst_cell_dstar_euler` | **S5 δ* 一発補正**: run_0011 と run_0002 (Euler) の質量流束欠損差から δ*(x) を抽出 (ランプ 0.009→0.11 H, カウル 0.006→0.013 H) し壁を法線オフセット、cell SST で再評価 | **C_T(p) 0.9685 (不変) / C_T(p+τ) 0.9602 / C_L 0.146 (無補正 0.154 → Euler 0.143) / C_M −0.930 (−0.980 → Euler −0.939)**。壁圧が MOC の破線に全域で乗る (`wall_pressure_cfd_vs_moc.png`) | active (ref) |
| `run_0010_moo_euler_2op/` | **S6 多作動点 MOO** (cell Euler, 作動点 cruise M∞6/p 0.05 p_in [w 0.6] + accel M∞4/0.15 [w 0.4], dv 5 個, LHS 12 + EHVI infill 3×2) | 18 評価 / 17 PASS / 1 INFEASIBLE、HV 1.179。パレート 8 点 (L 3.9–12.5 H, C_T,w 0.960–0.977)。`pareto.json`, `pareto.png`, `ledger.jsonl`, 各点 `doe_NNN_<op>/`, `inf_II_J_<op>/` | active (ref) |

| `run_0014_smoke_sst_node_thick` | node+SST 再挑戦 1: カウル板厚 0.2 % H (入口側で 0 に絞る台形) + interp 移植 | **soft 段 step 4 で発散** — 入口角 (x=−L_up, y=0) の 1 ノードが inlet×2 + wall×2 の 4 境界を持つ形になったため | 破棄予定 |
| `run_0015_smoke_sst_node_thick_idx` | node+SST 再挑戦 2: 板厚 0.2 % H (入口から一定、TE 手前で 0) + **stage 間の場移植を index コピーに変更** | **完走**: C_T(p) 0.9691 / C_T(p+τ) 0.9626 (摩擦 −0.0065) / C_L 0.155 / C_M −0.978、力 STEADY。rms_roOmega は本段開始直後から 3e18 一定 (壁ノード ω ピン留めの診断値、場は健全) | active (ref) |
| `run_0016_smoke_sst_node_t0_idx` | node+SST 再挑戦 3: **板厚 0 (元のスリット) + index コピー移植** — 真因の切り分け | **完走**、run_0015 と同値 (C_T(p) 0.9691 / C_L 0.155 / C_M −0.981) → **真因は interp_field の最近傍移植が座標一致の双子壁ノードを同じ元ノードに写していたこと** (板厚は不要) | active (ref) |
### S4(b) NASA TM X-71972 傾向照合のまとめ (2026-09-04, run_0003–0007, 図 `nasa_trends.png`, 表 `nasa_trend_table.py`)

- **内面 (ランプ + カウル内面) の力は forge Euler と MOC が全 5 形状で C_T +0.0006〜+0.0012、C_M 0.01〜0.05 以内で一致**。差の残りは
  カウル外面 (M∞10 の外部流がカウル角部で圧縮される衝撃圧) で、MOC は扱わない。カウル角 12° では外面だけで C_T −0.035。
- **カウル長** (θ_c 6°): 2.0 → 3.12 → 4.5H で C_T 0.9605 → 0.9704 → 0.9720。短縮で推力が大きく落ち、延長の追加利得は小
  — NASA の結論と同じ。ピッチモーメントは基準点依存: 入口下端基準では −0.60 / −0.69 / −0.64 だが、機体 CG 相当の前方基準
  (x_ref = −20H) では **−1.31 / −0.68 / −0.36** となり、NASA の「短いカウルで大きな頭下げ、延長で頭上げ増分」を再現する。
- **カウル角** (L_cowl 3.12H): 3° → 6° → 12° で C_T 0.9786 → 0.9704 → 0.9314 (単調減少、NASA の「角を開くと推力劣化」と同じ)。
  C_M は −0.90 → −0.69 → −0.29 (入口基準) で角を開くほど頭下げが緩む。NASA が 6° を最良としたのは飛行域全体のトリム舵角で、
  1 作動点の本比較では順位付けできない (S6 の多作動点束ねで評価する量)。
- 残差は全 run で 1.3〜1.5 桁プラトー (せん断層・衝撃の cell 床)。力係数は 500 step 出力で 1e-5 以内 (STEADY)。

### S5 / S6 のまとめ (2026-09-04)

- **RANS 化**: node + SST は本段でカウル角部直下流 (内面側 x≈0.13H) から ω が発散 (cfl 4/1 とも)。cell + SST (壁関数) は
  cfl 2 で安定し、圧力推力は Euler 比 +0.25 %、摩擦で −0.85 % (twall は forge の「流体に働く traction」規約、
  `viscousFlux_d.cu` L527 で符号確認)。node 側は未解決の課題として plan に残す。
- **δ* 一発補正 (S5)**: 単独の 0.99 max 縁判定は SERN の非一様コア (膨張扇) を欠損と誤認して 20 倍過大になる。
  **Euler 場を基準にした質量流束欠損** (`deltastar_sern_vs_euler`) で妥当な δ* (ランプ末端 0.11 H) が取れ、
  オフセット壁の RANS は壁圧・揚力・モーメントを非粘性設計値に戻した。推力は δ* に鈍感 (0.9685 不変)。
- **MOO (S6)**: 1 点 ≈ 135 s (2 作動点、GPU 共有時)。パレートは「長いランプほど C_T,w 高い」の単調前線
  (3.9H で 0.960 → 12.5H で 0.977) で、短い側は θ_r0 14–16°・M_c 3.2 付近、長い側は M_c 3.9・f 0.47。
  C_M (−20H 基準) は前線上で −0.19〜−2.6 と大きく変わり、カウル角 2° 台が頭下げ最小 — C_M を制約に入れると前線の選択が変わる。

### node + SST 発散の真因と解決 (2026-09-04, run_0014–0016)

- run_0009 の本段初期場 (soft 段末尾を `interp_field.py` で移植したもの) を調べると、カウル内面と外面の壁ノード (座標が一致する
  スリットの双子) が全 station で**同一の圧力**を持っていた (排気側 15 kPa のはずが外部流側の 2 kPa など)。合成場で
  `interp_field.py` を試すと双子 134 station 全部が同じ元ノードに写った = **最近傍補間が座標一致ノードを区別できない**。
  壁ノードと隣接内部ノードの 6 倍の圧力段差が 2 次で step 7 に爆発した。solver 側の欠陥ではない。
- 対策: stage 間の同一メッシュ移植を `runner_sern.restart_by_index` (VALUE の index コピー) に変更。板厚 0 のスリットのまま
  node+SST が完走 (run_0016)。板厚オプション (`mesh.cowl_thickness`) は残すが必須ではない (run_0015 と同値)。
- node の twall は「壁ノードに働く力」(cell は「流体に働く力」) で符号が逆。`sern_forces` は離散化で規約を切替 (摩擦 node −0.0065、
  cell −0.0085)。
- 残る node 固有の観察: (i) rms_roOmega が本段で 3.3e18 一定 (壁ノードの ω ピン留め残差が混入する診断値。場は STEADY で ω max 1.5e7 は
  外部流側の壁ノード)、(ii) 入口角 (inlet+wall) と後縁・ランプ後縁の**単ノードの圧力が外れる** (cowl_in 入口角 2.2 p_in、cowl_out 入口角
  0.36 p_in)。力積分への影響は小 (node/cell の C_L 差 0.001) だが、壁圧分布を読むときは端点を除く。

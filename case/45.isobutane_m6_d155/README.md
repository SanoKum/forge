# case/45 — イソブタン燃焼ガス M6 風洞・出口径 1.55 m の最短全長スタディ

Pt 5.5 MPa / Tt 1600 K / φ=0.9 燃焼ガス (semi-perfect NASA-9) / M_design 6 /
**出口径 1.55 m は δ\* 込み物理壁で定義**。全長 (スロート→物理出口) の最短化と
粘性壁 (δ\* 物理壁 + SST NS) までの検証。計画:
[`plans/active/design-isobutane-m6-d155.md`](../../plans/active/design-isobutane-m6-d155.md)。

- 最短探索: `search_shortest.py` → `search_shortest.json`
  (基準 = shortest-robust study: margin≥1°/topo≥0.02/hard gate/単峰、n1200 確認)。
  **勝者 R2 / L_c 39.3 / M_K 2.7 → x_F = 95.104 r_t** (R3: 39.7/2.8 → 95.433)。
- r_t = (0.775 − δ\*_exit)/(r_F/r_t), 初期推定 0.0771 m (δ\*_exit ≈ 0.66 r_t)。
  A/A\*(M6, Tt1600) = 88.21。
- 凝縮事前見積り: 出口 T 233.5 K vs Tsat(H₂O) 263.9 K = **−30 K 過冷却** →
  dry 本体 + 採用点の凝縮 ON 再評価 (方針 a)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_euler_shortest_dry` | 最短点 (R2/L_c39.3/M_K2.7) の node Euler TP dry (`problem_d155_R2_Lc39.3_tp.yaml`, split_h2o, 段階起動 12000 step) | **dM_max 0.036 % M_d・軸出口 M 6.00003・overshoot 0.034 %・ε_M_rms 0.0057 %・ε_θ 0.0054°・コア 64/65 面**。品質 PASS (AR OK/skew 0.44)・NaN 0・quasisteady ALL STEADY・軸 M 凍結 8k→12k 5.6e-5。残差は warm 床 (rms_ro 1.1e-6, `check_convergence` NOT CONVERGED 判定 = M6 系列の既知パターン)。`metrics.json`/`residual_history.png` | active (**最短点 Euler 検証の正本**) |
| `run_0002_ns_coarse` | NS 中継 (y+~50, `problem_d155_ns_coarse.yaml`, 物理壁 δ\* v1 相関, IC=run_0001) | 完走 12000 step・NaN 0。残差 2–3 桁降下 (ro/roe は低レベルでプラトー = 中継品質として想定内)。出口壁半径 0.7771 m (v1 相関 δ\* ≈ 0.71 r_t) | active (中継) |
| `run_0003_ns_v1` | NS 本計算 v1 (y+~1.4, `problem_d155_ns.yaml`, 相関 δ\*, 24000 step, IC=run_0002) | 完走・NaN 0・品質 SOFT-PASS。**軸出口 M 6.00098・dM_max 0.327 % M_d (ゲート内)**・ε_M_rms 0.075 %。残差 2 桁降下 warm 床。δ\* v3 抽出元 (`dstar_v3.csv`: 相関は中流過大 median 比 0.909, x14 で 0.669) | active (v1 記録・v3 抽出元) |
| `run_0004_ns_v3` | NS 本計算 v3 (CFD 抽出 δ\* 全域採用, IC=run_0003) | 完走・NaN 0・品質 SOFT-PASS・16k→24k 軸 M 凍結 1.2e-4。**δ\* 固定点確認 (2巡目比 median 1.001 [0.981,1.005])**。出口面軸 M **5.9856 (−0.24 % M_d)**、x_E ディップ 5.956 (law 側残差)、一様区間 5.975–5.996 のうねり。metrics の `M_axis_exit` は x_E 評価値であることに注意 | active (**dry NS 固定点・トリム根拠**) |
| `run_0005_euler_trim` | **粘性トリム版** Euler (`problem_d155_trim_tp.yaml`: Md 6.0144, L_c 39.7/M_K 2.7 [トリム後最短 x_F 96.00 r_t], r_t 0.076511) | dM_max 0.036 % M_d・x_E で M 6.01445 (=目標)・ε_M_rms 0.0056 %・コア 64/65。NaN 0 | active (トリム Euler) |
| `run_0006_ns_trim_v1` | トリム版 NS v1 (相関 δ\*, `problem_d155_trim_ns.yaml`, IC=run_0004 cross-mesh) | 完走・NaN 0・凍結 2.7e-4。出口面軸 M **5.9810** (= 旧 v1 5.9667 + トリム 0.24 % — トリムは線形に作用)。δ\* v3' 抽出元 | active (トリム v1・抽出元) |
| `run_0007_ns_trim_v3` | トリム版 NS v3 (CFD 抽出 δ\*, IC=run_0006) | 完走・NaN 0・凍結 1.2e-4・**δ\* 固定点 (2巡目比 median 0.999)**。**出口面軸 M = 6.00079 (+0.013 % M_d)** — トリム狙い通り。quasisteady machmax/pmax STEADY (shock 指標は無衝撃のため対象外)。出口壁半径 0.77591 m → 最終 r_t 補正 −0.12 % で 0.775 m 合わせ | active (**dry 最終形の正本**) |
| `run_0008_ns_trim_cond` | 凝縮 ON restart (`problem_d155_trim_ns_cond.yaml`: Kw+HK condModel1+Kantrowitz, 蒸発 ON, IC=run_0007, 12000 step) | 完走・NaN 0・**STEADY** (4k/8k/12k で M_exit 差 5e-4)。軸 onset x≈69 r_t、出口 g 0.20 % (H₂O の 2 %)、**出口軸 M 5.9273 (−1.2 %)**・試験区間に M 低下勾配 (x60→96 で 6.00→5.93)。dry の軸は x≈24 r_t (M5.5) で飽和線越え S≈14 (`axis_values.csv` の Tsat_post) | active (**凝縮評価の正本**) |

## 結論 (2026-08-31)

- **最短全長: スロート→物理出口 7.337 m** (x_F 96.00 r_t, r_t 76.42 mm)、入口配管端→出口 8.292 m。うち E→F 一様化区間 ≈4.36 m は物理の床。
- **出口軸 M (dry): 6.0008 (+0.013 %)** — Md 6.0144 トリム (NS 固定点の −0.24 % 欠損を打ち消し) で達成。
- **凝縮リスク**: Tt 1600 K では dry 不成立 (S≈14)。凝縮 ON で出口 M 5.927・軸勾配残り。凝縮フリーは Tt ≳ 1820 K。
- 結果ページ: https://claude.ai/code/artifact/f1cfbdf8-2415-4bfc-ae2d-156793569bd0

注: `run_0006_ns_trim` (旧設計 run_0004 の δ\* CSV を新設計に流用するショートカット) は物理壁フィルタ不合格 (スロート下流非単調) で prepare 段階に失敗し**削除済み**。トリム版も正規フロー (v1→抽出→v3) で回す。

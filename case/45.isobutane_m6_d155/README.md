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

## SU2 クロスチェック (境界層厚さ, 2026-09-01)

δ\* 抽出の軸対称化 (円環重み) と合わせて、**同一メッシュ・同一 BC・同一ガス (CPG γ\*=1.27354)** で forge node SST と SU2 v8.5 axisym SST の境界層厚さを比較する ([procedures/su2-cross-check.md](../../procedures/su2-cross-check.md) 準拠)。

| run | 内容 | 結果 | 状態 |
| --- | --- | --- | --- |
| `run_0009_cpg_euler` | CPG 版 Euler (`problem_d155_cpg_euler.yaml`) | 完走・NaN 0。CPG は A/A*=159 (semi-perfect 88.2) でノズルが伸びる (出口壁半径 1.03 m) — 比較チェーン専用で製品形状ではない | active (比較チェーン) |
| `run_0010_cpg_ns_coarse` | CPG NS 中継 (y+~50) | 完走・NaN 0 | active (中継) |
| `run_0011_cpg_ns` | **CPG NS 本計算** (wall_first 6.5e-5=y+~2 [4.5e-5 は CPG 長尺で AR1403 FAIL→緩和], 相関 δ\* 物理壁, 24000 step) — 比較の forge 側 | 完走・NaN 0・品質 PASS。残差 2 桁 warm 床 (既知パターン) | active (**比較 forge 側**) |
| `run_0016_chamber_cpg` | **チャンバー A/B** (Codex レビュー準拠): 完全形状 (入口管+収縮+ノズル) + 出口後方 4D_e×外径 3D_e チャンバー+フランジ壁の 4 ブロック transfinite (`gen_chamber_mesh.py` → `mesh_chamber/`)。品質 PASS (AR296)。静止雰囲気 1389 Pa IC + 段階起動 18000 step で plume 確立 (ジェット軸 M5.97, 外周静止)。**probe (cfl4 OK/8@181/16@21) は素メッシュと同一挙動 = 出口境界は無罪**。**line-implicit 併用の 4 象限も全て同一** (coarse×line 8@216/16@15, chamber×line 8@194/16@15; チャンバーのライン被覆 81.9%・噴流域 point fallback) — 壁法線×出口境界の直積消去で streamwise 内部モード律速が確定 | active (**A/B の正本**) |
| `run_0012_su2_sst` | **SU2 v8.5 RANS-SST** (同一メッシュ msh→su2 変換, axisym, ROE+MUSCL, SST V2003m, Sutherland/Pr0.72/Prt0.9, 40k+20k iter) | 収束: rms[Rho] −6.5・rms[RhoE] −2.3・出口 massflux/Mach ドリフト 0.003 %/5k (手順書合格実績以上)。`history.csv` は継続分のみ | active (**比較 SU2 側**) |
| `run_0013_cpg_ns_plainsst` | forge 素 SST 化 A/B (`dilatationCorrection: 0, katoLaunder: 0` — prepare 後に solverConfig 直接編集。YAML への独自キー追記は黙って無視されるため不可), IC=run_0011, 12000 step | 完走・NaN 0。**素 SST にすると SU2 と一致: δ99 ≤3 %・θ ≤0.4 %・δ* ≤1.3 %** | active (**帰属 A/B の正本**) |
| `run_0014_cpg_ns_dilat0_kl1` | 帰属分離 (dilatation OFF / KL ON) | run_0013 と ≤0.7 % 差 → **Kato-Launder は無関係、差は dilatationCorrection 単独** | active (帰属分離) |

### SU2 クロスチェックの結論 (2026-09-01)

- **forge 生産設定 SST は SU2 比で境界層が薄い**: δ99 −17〜21 %・θ −16 %・δ* −5 % (4 ステーション一貫)。
- **帰属は `dilatationCorrection: 2` (圧縮性生産項の正確形: deviatoric trace 除去 + 等方項 −⅔ρk∇·u) 単独**。素 SST 化した forge は SU2 と δ99 ≤3 %/θ ≤0.4 %/δ* ≤1.3 % で一致 = **ソルバ・離散化は無罪、意図した乱流モデル差**。Kato-Launder の寄与 ≤0.7 %。
- ノズルは全域 ∇·u>0 (強膨張) のため等方項が k のシンクとして常時働く。急膨張の乱流抑制は物理的に実在する効果で forge 形は正当だが、**SST のモデル形式差として境界層厚に ±16 % 級の不確かさ**があると認識すること。δ* への影響は ±5 % (出口 M ±0.1〜0.15 %) に留まり、Md トリムは forge 自身の NS で較正しているため製品性能の結論は不変。
- 比較図・数値: `compare_bl_su2.{png,json}` (3 者重ね描き)、40k 時点の記録 `compare_bl_su2_iter40k.json`。

## 陰解法 cfl_pseudo 上限スイープ (2026-09-02)

`run_0015_cflsweep/` (probe 22 本, `sweep_cfl_implicit.py` + `summary.json`。中間 res は削除済み・各 probe の config/residual/NaN ダンプは保持)。warm 場 (run_0007 TP / run_0011 CPG) restart 2000 step + coarse IC 収束レース 4000 step。

- **素の上限 (implicitRelax=1): TP は cfl 2 (3 で発散)、CPG は 3 (4 で発散)**。既定挙動としては単桁前半で頭打ち。
- **原因は低マッハ調整ではない** (仮説棄却): ① 現行 config は `lowMachPrecond: 0` で低マッハ処理は不活性、② 発散箇所は亜音速チェンバでなく**下流端の超音速壁 BL 内** (壁から ~0.003 r_t の対数層)、③ `lowMachPrecond: 2` を入れると逆に悪化 (step 97→13)。
- **終端症状 = 陰的更新の P アンダーシュート → EOS 床洗浄 → NaN** (ω 爆発は増幅表示)。**帰属の最終確定 (2026-09-02, チャンバー A/B 完了)**: 発散種の位置はメッシュ依存 (fine=壁∩出口コーナー / coarse=出口手前コア / 高cfl=スロート軸上) で、①背圧値 100/10 Pa・②逆流 Tt 1500K・③出口近傍局所 dt キャップ・④**チャンバー付加 (出口境界を 8 m 後方へ, run_0016)** の全てに**不感** (発散 cfl・step ほぼ不変: chamber cfl4 OK/8@181/16@21 vs 素 cfl4 OK/8@215/16@17)。⑤ **乱流凍結 (`FORGE_FREEZE_TURB=1`, k-ω 更新停止・μt 固定) でも発散 step が 1 step 単位で一致** (cpg@98/28/18, tp@53/24 = baseline と同一) — SST は完全に無関与で、ω 先行爆発は凍結方程式の残差が壊れた流れ場を映す純粋な症状と**証明**。⇒ **上限の律速 = 平均流のみの内部 streamwise defect-correction 反復不安定** (6 系統の消去で確定)。根治候補は streamwise 第 2 ライン族 (ADI) / LU-SGS 流下順序のみ。SST 連成陰化は不要と事前確定。
- **切り分け**: `nStepInner` 10/20 無効 (DPLUR 内部収束でない)・`implicitRelaxSST` 0.5 無効 (SST 方程式でない) — **`implicitRelax` (流れ方程式の緩和) だけが効く**。
- **緩和込みの上限**: relax0.7 → cfl 8 (12 で発散)、relax0.5 → cfl 16+。積 cfl×relax ≈ 6-8 で飽和。
- **収束レース**: 実効速度も cfl×relax でスケールし **cfl8+relax0.7 が最速** (4000 step で rms_roUx 現行 cfl1 比 1/3、rms_roOmega 1/4)。cfl16+relax0.5 は頭打ち。図 `run_0015_cflsweep/race_residuals.png`。
- **生産推奨**: CPG `cfl 8 + implicitRelax 0.7` / TP `cfl 4-6 + implicitRelax 0.7` (warm 検証済上限 6 に対しマージン)。現行 cfl1 比で 3〜4 倍速の見込み。
- **line-implicit 追試 (2026-09-02, probes7–10, [plan](../../plans/accepted/time_integration-line-implicit.md))**: 壁法線ライン block-Thomas (`lineImplicit: 1`) を実装して検証。ライン構築は完璧 (1250 本・被覆 100 %) だが **cfl 上限は不変 (積 ≈6–8)** — cfl8-line の発散は出口域全断面で、**律速は壁法線でなく streamwise lag** と判明 (壁法線剛性は粘性対角が既に吸収)。効果は同一設定の収束 −14 %/step のみ。実装過程で lu5 ピボット罠 (LASWP 先行必須)・device printf 引数上限罠を踏んで記録。
- **追試 (2026-09-02, probes5/6)**: commit 時の正値性ガード (`updateGuardAlpha`, 新実装 opt-in) は**上限を上げない** (発散遅延のみ。α=0.5 は毎 step 半減を許すため床への歩行が続く)。**lowMachPrecond=2+ガード併用も不成立** — P が床に触れないまま SST チャネル (ω 2.8e22) で爆発。⇒ 真因は「defect-correction 反復 (近似 Jacobian + SST 隔離更新) の不安定」で、床 NaN は終端症状。`implicitRelax` が効くのは反復スペクトル半径を縮めるから。詳細 [plan](../../plans/accepted/time_integration-update-positivity-guard.md)

## 結論 (2026-08-31)

- **最短全長: スロート→物理出口 7.337 m** (x_F 96.00 r_t, r_t 76.42 mm)、入口配管端→出口 8.292 m。うち E→F 一様化区間 ≈4.36 m は物理の床。
- **出口軸 M (dry): 6.0008 (+0.013 %)** — Md 6.0144 トリム (NS 固定点の −0.24 % 欠損を打ち消し) で達成。
- **凝縮リスク**: Tt 1600 K では dry 不成立 (S≈14)。凝縮 ON で出口 M 5.927・軸勾配残り。凝縮フリーは Tt ≳ 1820 K。
- 結果ページ: https://claude.ai/code/artifact/f1cfbdf8-2415-4bfc-ae2d-156793569bd0

注: `run_0006_ns_trim` (旧設計 run_0004 の δ\* CSV を新設計に流用するショートカット) は物理壁フィルタ不合格 (スロート下流非単調) で prepare 段階に失敗し**削除済み**。トリム版も正規フロー (v1→抽出→v3) で回す。

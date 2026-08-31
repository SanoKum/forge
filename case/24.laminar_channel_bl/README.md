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
| `run_isoT_cond_{cell,node}` | **等温壁検証① (純伝導・過渡)**: `mesh/conduction.geo` (上壁 300K/下壁 350K `wall_isothermal`, 側面 slip, 静止空気)。一様 300K IC から implicit 20k step | 音速基準擬似 dt では拡散時定数 (~5s) に届かず**過渡途中** (判定不能)。→ 厳密解シード方式 (isoT_condL) へ切替 | 削除済み(2026-08-31) (方式転換) |
| `run_isoT_condL_{cell,node}` | **等温壁検証② (厳密解シード)**: 線形 T 厳密解を IC に張り「不動点として保持されるか」を判定 | cell: ドリフト max 0.20 K・静穏。**node: slip 側面上に市松状スプリアス接線流 0.47 m/s が定在** (壁ノードは 0)、その移流で T 最大 2.2 K 汚染 → slip 境界が犯人と切り分け (isoT_condP で確定) | ref (node slip バグの証拠) |
| `run_isoT_condP_{cell,node}` | 同③: **側面を periodic に変更** (slip 排除の A/B)。10k+40k step で固定点確認 (10k→40k で 4 桁一致) | node のスプリアス流 0.47→**3.6e-4 m/s**・T ドリフト 0.11 K = **slip が犯人と確定、等温壁自体は node でも機能**。**中央 50% の熱流束は両モードとも厳密 (−0.03/−0.05%)**。壁隣接のみキンク: cell は壁 BC-第1点勾配 +94% ≈ **係数 2 バグの指紋** (→ ghost T 修正)、node は壁ノードが CV 平均値 (350→349.894) に緩み第1スペーシング −24% | ref (等温壁検証の主根拠) |
| `run_isoT_condR_{cell,node}` | 同④: ny=64→128 細分化 (キンクの次数判定) | cell キンク 0.20→0.099 K = **壁距離に比例して半減** (抵抗 2 倍バグと整合)。node は擬似 dt 減で 50k step では固定点未達 (参考値) | ref (次数判定) |
| `run_isoT_condF_{cell,node}` | 同⑤: **ghost T 鏡像外挿修正後** (`boundaryCond_d.cu` `wall_isothermal_d`: T_ghost = 2Tw−Tc) の再確認 | **cell: 線形厳密解が機械精度の不動点に** (max\|ΔT\|=0.0000 K, 壁 q_w +0.02%, max\|U\|=2.7e-7)。node: 修正と独立で不変 (想定どおり)。回帰スイート 4/4 PASS | ref (修正検証) |
| `run_isoT_condL_node_m5` / `run_isoT_condL_node_1st` | **slip 市松の次数切り分け A/B**: 同一バイナリ (T ピン込)・cfl_pseudo 5 で MUSCL vs 1 次風上 (`convMethod: 0`)、slip 側面のまま 10k step | 側面列 max\|U\| **0.476 vs 0.446 m/s・符号反転率 0.48 vs 0.41 でほぼ不変** → 市松は再構成次数でなく **slip 閉包自体** ([調査ノート](../../notes/investigations/node-slip-tangential-density-spurious-flow.md) 追試) | ref (slip 診断 A/B) |
| `run_couette_cell` | **粘性加熱 Couette (cell・厳密解)**: conduction メッシュ periodic 側面、両壁等温 300K・上壁 `Ux: 100` (移動壁 ghost 一般化 2U_w−U_c の検証)。厳密解 (u 線形・T 放物線 max 0.94 K・q_w=μU_w²/2H=9.05 W/m²) シード、implicit cfl50 20k step | **厳密一致**: u ドリフト 0.012%Uw・T ドリフト 0.0004 K・**q_w 下/上 +0.17/−0.02%・τ_w −0.02%**。移動壁 ghost + 粘性仕事項 τ·U_b + 発熱-壁排熱バランスを対流条件で検証。回帰 4/4 PASS (U_w=0 経路ビット同一) | ref (移動壁+粘性加熱検証) |
| `run_poiseuille_bf_{cell,node}` | **体積力 config の厳密解検証**: `bodyForce: [7.24,0,0]` 駆動 Poiseuille (periodic 側面・両壁等温 300K)。厳密解 (u 放物線 u_max=5 m/s・T 四次 ΔT=6.3e-3 K・q_w=0.1207 W/m²・τ_w=fH/2) シード 20k step | **cell: u 0.30%・T 1.5e-4 K・τ_w 一致 (評価 FD バイアス −2% 込みで厳密)**。node: q_w −0.0%・第 1 内点 u −20% = 既知 W-I 閉包 O(Δ) バイアス (体積力は無罪)。回帰 4/4 PASS (off ビット不変) | ref (bodyForce 検証) |
| `run_isoT_condN_node` | 同⑥: **node 壁ノード T ピン実装後** (状態ピン + res_roe ゼロ化 + block-DPLUR エネルギー行 decouple `iso_wall_flag`) の検証。explicit RK3 / implicit cfl_pseudo 5 / 同 20 (プローブ, 破棄済) | **壁ノード T = BC 厳密 (350.0000/300.0000)・q_mid +0.00%**。explicit 3000 step 健全。implicit は **decouple 無しだと数 step で発散** (Jacobian 不整合)・decouple 込みで cfl_pseudo 5 安定 (有界リミットサイクル) / **20 は発散 = 実用上限 ~5**。第 1 スペーシング勾配は −24%→−15% (残りは W-I 弱形式閉包の O(Δ) バイアス・別課題)。**2026-07-20 後段: nodeWallDirichlet_d への移設リファクタ後の再実行 (一時 run, 破棄済) で物理指標一致** (壁 T 厳密・q_mid ±0.000%・ドリフト 0.087 K、差はコンパイル単位変更の ulp 級のみ) | ref (node T ピン検証) |

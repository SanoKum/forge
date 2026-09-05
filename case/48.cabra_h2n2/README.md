# case/48.cabra_h2n2 — Cabra H₂/N₂ 浮き上がり火炎 (vitiated coflow, 燃焼加熱器の検証ケース)

Cabra et al. (UC Berkeley; NASA/CR-2004-212887 Table 6.1, [`papers/combustion/`](../../papers/combustion/)): 中心ジェット d=4.57 mm
(管 OD 6.35 mm) から H₂ 25 %/N₂ 75 % (体積)・305 K・バルク 107 m/s (Re 23,600) を、希薄 H₂/空気予混合火炎の燃焼ガス coflow
(D=210 mm, 1045 K, 3.5 m/s, X O₂ 0.15 / H₂O 0.099 / N₂ 0.751) へ噴射。自着火で浮き上がり高さ H/d≈10、火炎長 30 d。
公開データ (TNF): 中心軸 (z/d=1–34) と半径分布 (z/d=1,8,9,10,11,14,26) の Favre 平均 T・主要種・OH・NO。ここでは NASA 報告書の図を
目視で読み取った `exp_centerline.csv` / `exp_radial_T.csv` を参照にしている (±3 %)。加熱器の「高温 vitiated 流への H₂ 噴射→自着火」と同じ物理で、
**低マッハ前処理 + 反応** の経路と着火位置予測を検証する。計画: [`plans/active/chemistry-finite-rate-h2.md`](../../plans/active/chemistry-finite-rate-h2.md)。

## セットアップ

- `mesh/make_cabra_mesh.py`: 軸対称平面構造 quad (node)。ジェット管を上流 20 mm 延長 (管内 no-slip、入口は 1/7 乗則プロファイル)、
  リップ 0.89 mm、外周 R=100 mm slip、x=0–250 mm。60k ノード、品質 PASS (AR max 70, skew 0)。
- `setup_cabra_case.py run_dir [--chem 0|1] [--ji 5] [--tci 0|1] [--cfl] [--conv] [--relax] [--nstep] [--restart]`: Cantera で入口状態、
  BC (`inlet_uniformVelocity`×2 + profile, `outlet_statPress`, `wall`, `slip`, `axis`)、config (node 軸対称, `nodeInletCornerWall: 1`,
  `lowMachPrecond: 2`, implicit, SST, TP 9 種, chemistry `jacobianInterval 5`)、IC (coflow 一様 + ジェット柱)。
  切り分け用オプション: `--single 1` (N₂ 単一種), `--jetcof 1` (ジェット組成=coflow), `--jetn2 1` (純 N₂ ジェット), `--iccol 0` (ジェット柱 IC なし=ピストン起動),
  `--planar 1` (平面・軸→slip), `--precond 0|1|2|3`, `--coupling 0|1|2`, `--far slip|outlet_statPress`, `--sfr`, `--nojet`。
- `compare_cabra.py run_dir step out.png`: 中心軸 T/Y と半径 T を実験と重ねる。

## 計算 run 一覧

このケースで **forge の低マッハ node 経路のバグを 3 件特定・修正** した (2026-09-05, 詳細は
[plan lowmach 変更ログ](../../plans/accepted/time_integration-lowmach-preconditioning.md) と
[plan outlet §2.12](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)):
(A) node × `lowMachPrecond>=2` の前処理 block カーネルが境界ノード (軸・壁・出口・遠方) を凍結、
(B) TP 亜音速 `outlet_statPress` の γ 混用で低マッハ出口が一様逆流、
(C) 多成分 × 定常擬似時間 × precond で組成前線と密度前線の擬似速度が不整合 → `speciesPrecondDt` (既定 1)。
`run_0001`–`run_0035` は (A)(B)(C) すべてを含む旧バイナリ、`run_0036`–`run_0046` は (A) 修正後、`run_0047`–`run_0061` は
(A)(B) 修正後、`run_0062` 以降が (A)(B)(C) 修正後。全て化学 OFF (混合のみ)、3000 step 以下の切り分け run は破棄予定。

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_mix_start` | 段階起動 1 次 cfl_pseudo 0.5 (旧バイナリ) | 完走するが後述の凍結境界・偽逆流を含む | 破棄予定 |
| `run_0002`–`run_0013` | 2 次化 / cfl 2 / precond 0 / coflow のみ / 単一種 / ジェット柱なし / coupling 1,0 / 平面 の切り分け (旧バイナリ) | 多成分は全て P 70–230 kPa・T>coflow の非物理場 (`run_0003` outlet∩farfield で NaN)。coflow のみ (`run_0008`) と単一種 N₂ (`run_0009`, 柱 IC) は健全 → 「多成分固有」と誤認した (実は境界凍結 (A) が単一種では見えなかった) | 破棄予定 |
| `run_0014`–`run_0029` | 400–800 step 短期診断: 有界性 / 1 次 / sfr / 純 N₂ ジェット / precond 0,1 / relax 1 / RK3 (局所 dt・時間精度) / cfl 0.1 / 一様組成 | 化学種は有界 (Y_H2 ≤ 入口値, ΣY=1)。時間精度 RK3 (`run_0024`/`run_0025`, 6000 step 120 µs) は健全 = 空間スキーム無罪。implicit は初期 400 step の圧縮過渡は物理 (ピストン) | 破棄予定 |
| `run_0030`–`run_0035` | 単一種ピストン起動で軸ノードの追従を確認 (precond 2 / axisymMethod 1 / centroidShift 0 / precond 0) | **precond 2 では軸ノードが IC (1045 K, 4 m/s) のまま隣接 (305 K, 80 m/s) に追従せず、precond 0 (`run_0035`) と explicit (`run_0026`) は追従** → 真因 (A) を特定 | 破棄予定 (エビデンス) |
| `run_0036`–`run_0046` | (A) 修正後: 単一種/多成分, 柱 IC あり/なし, precond 3, cfl 0.1, farfield=outlet, 出口列の早期追跡 | 軸は追従するが**出口列が全体で一様に逆流** (−50 m/s @100 step, P 単調上昇, 平面・precond 0 でも同一 `run_0045`/`run_0046`) → 真因 (B) を特定 (`run_0044` 出口列時系列) | 破棄予定 (エビデンス) |
| `run_0047_fix2_outlet_early` | (B) 修正後の出口列 200 step | 出口列 3.5 m/s / 101.33 kPa を保持 | 破棄予定 (エビデンス) |
| `run_0048_fix2_single_col` | (A)(B) 後, 単一種 N₂, 柱 IC, cfl 0.5, 3000 step | **健全**: P 101.0–102.0 kPa, T 305–1045 K, 軸追従 | ref (単一種基準) |
| `run_0049`–`run_0057` | (A)(B) 後の多成分: coupling 0/1/2, 一様組成 9 種, 平面, relax 1, control 0 | 異組成ジェットは**軸上ジェット核で P 148 kPa・T 1150 K に暴走** (1000 step)。一様組成 (`run_0053`) は健全、coupling 依存なし、relax 1 で悪化 → 多成分固有の第 3 の問題 | 破棄予定 (エビデンス) |
| `run_0058`–`run_0061` | dual-time (`unsteady 1, dualTime 1`, dt 2e-7 / 1e-6) と precond 0 定常 | **dual-time (`run_0059`/`run_0061`, 2 ms) と precond 0 (`run_0060`) は健全** → precond × 分離化学種の擬似時間不整合 (C) を特定 | ref (dual-time 代替の実証) |
| `run_0062_spdt1` / `run_0063_spdt2` / `run_0064_spdt0_bitcheck` | (C) `speciesPrecondDt` 1 / 2 / 0 の A/B (多成分, 柱 IC, cfl 0.5, 1500 step) | **mode 1・2 は P 100.2–102.6 kPa, T≤1045 K で健全**; mode 0 は 91–148 kPa / 1144 K で暴走 → mode 1 を既定に | ref (C の検証) |
| `run_0065_mix_col_cfl1_6k` | 全修正後の混合 (化学 OFF): 2 次, cfl_pseudo 1, relax 0.5, 柱 IC, 6000 step | 有界 (P 99.3–102.2 kPa, T 304–1045 K, NaN なし)、ジェットは混合中 (軸 Y_H2 x/d 20–50 で 0.0234→0.0143, T 655 K)。`check_convergence`: NOT CONVERGED (発達途上: rms_roUx/roe falling, rms_ro/roK/roOmega plateau) | ref (混合の起点) |
| `run_0066_mix_col_cfl2_cont` | run_0065 res_6000 から cfl_pseudo 2 で継続 | **step 932 で NaN 発散** (cfl_pseudo 2 は不可; 位置は `forge_run.log`) | 破棄予定 |
| `run_0067_mix_col_cfl1_cont` | run_0065 res_6000 から cfl_pseudo 1 で 20000 step 継続 (計 26000 step) | 有界 (P 101.2–102.1 kPa, T 304–1045 K, NaN なし)。軸 Y_H2: x/d 10 0.0178 / 20 0.0076 / 30 0.0046 / 50 0.0044 (下流はまだ発達中)。`check_convergence`: NOT CONVERGED (rms_ro/roe/roK/roOmega falling, rms_roUx/roUy plateau 4–6e-6)。次: さらに継続 → 反応 ON | active (混合場の最新) |
| `run_0068_react_cfl1` | **反応 ON** (jac 2, ji 5, coupling 2, tci 0), run_0067 res_20000 から cfl_pseudo 1, 10000 step | 有界。step ~6000 で x/d≈20 に自着火 (OH 3e-3)、火炎基部が x/d 10.7 (8000) → 7.7 (10000) と上流へ移動中、Tmax 1930 K。`check_convergence` NOT CONVERGED (過渡) | ref (着火の起点) |
| `run_0069_react_cont` | run_0068 res_10000 から 20000 step 継続 (tci 0 = 層流化学) | 火炎基部が x/d 4.6 (4000) → 3.4 (12000) → **0 (16000 以降, リップ付着火炎)**。Tmax 1800–2020 K, 軸 T max 1593 K (z/d≈30)。`check_convergence` NOT CONVERGED (全列 plateau)。TCI なし RANS の典型 (平均温度で層流反応率を評価し着火過大) → PaSR へ。図 `compare_cabra_20000.png` | ref (tci 0 の結論) |
| `run_0070_react_pasr` | 同条件で **PaSR (`tci 1`, tauChem 1, C_mix 1)**, run_0067 混合場から 30000 step | 着火は同じく x/d≈21 (5000)、火炎基部の上流伝播が遅くなるだけで 7.8 (10000) → 2.9 (20000) → 1.2 (30000) と**ほぼ付着**。軸 T は tci 0 と同様 (T>600 K が z/d 12.7, 実験 14)。`check_convergence` NOT CONVERGED。図 `compare_cabra_30000.png` | ref (PaSR C_mix 1) |
| `run_0071_react_pasr_cmix4` | PaSR C_mix 4 (τ_mix 4 倍 → κ 小) の感度, 30000 step | 推移は C_mix 1 とほぼ同一 (21.2 → 7.9 → 4.4 → 1.9 → 1.4 → 1.2)。**PaSR の κ は基部の上流伝播を止めない** (近傍せん断層の平均場が既に可燃・高温)。NOT CONVERGED。図 `compare_cabra_30000.png` | ref (PaSR 感度の結論) |
| `run_0072_react_li` | **機構 A/B (plan §5.1 P0-3)**: run_0068/0069 と同一条件・同一 IC (run_0067 `res_20000`)、機構だけ Li 2004 (`tools/mechanisms/h2co_li2004_cantera.yaml`) に交換, tci 0, 30000 step | **付着推移は Jachimowski とほぼ同一**: 着火 x/d (Y_OH>2e-4) 20.8 (6k) → 10.5 (8k) → 7.2 (10k) → 3.3 (16k) → **0.0 (20k 以降リップ付着で安定)** (Jach は累積 26k で付着)。**機構交換では付着は変わらない = 支配因子は TCI 欠如と確定**。NOT CONVERGED (過渡)、NaN なし | ref (**機構 A/B の結論**) |
| `run_0073_react_li_cont` | run_0072 `res_30000` (付着済) から Li 機構のまま +30000 step 継続。目的: 付着後の統計定常確認 (`check_quasisteady --quantity ignx,tmax,exit_massflux,exit_hflux,exit_y_H2O`)。run_0072 末尾は出口質量流束が 0.023→0.045 kg/s と回復途中 (環状逆流 y 54–78 mm) | (実行中) | active |

## 機構着火遅れ比較 (2026-09-05, plan §5.1 P0-3)

`ign_delay_mech_compare.py` (Cantera, Cabra 混合線・断熱定圧 1 atm・τ = max dT/dt) による Jachimowski 9sp20r の 1000–1100 K 妥当性検証。成果物: `ign_delay_mech_compare.{csv,png,log}`。

| 機構 | τ_min @T_c 1015 K | @1045 K | @1075 K |
| --- | --- | --- | --- |
| Jachimowski 9sp20r (forge 現行) | 4.41 ms | **1.32 ms** | 0.70 ms |
| Li 2004 | 6.82 ms | 1.54 ms | 0.83 ms |
| Burke 2012 | 28.99 ms | 2.73 ms | 1.21 ms |
| GRI3.0 H₂ subset (既知不良) | 94.0 ms | 14.0 ms | 2.44 ms |

ξ_MR は 0.025–0.045 (文献 ≈0.05 と整合)。**Jachimowski は Cabra 検証実績のある Li 2004 より 15–35 % 早く、Burke 2012 の 1/2 (1045 K)〜1/6.6 (1015 K)** — 検証済み機構に対し系統的に早着火側で、リップ付着・BK 早着火を悪化させる向き (支配因子は TCI 欠如)。GRI subset の 10 倍遅れは Benim 2020 の「使用不可」報告と整合し、比較系の健全性確認になっている。Burke/Li の Cantera YAML は `solver_density_cuda/tools/mechanisms/h2_burke2012_cantera.yaml` / `h2co_li2004_cantera.yaml` (forge 用変換は未)。次: Li で反応 ON A/B。

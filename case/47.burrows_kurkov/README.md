# case/47.burrows_kurkov — Burrows–Kurkov 超音速 H₂ 燃焼 (有限速度化学 Phase 3 検証)

NASA TM X-2828 (Burrows & Kurkov 1973, [`papers/combustion/`](../../papers/combustion/)): 加熱空気 (vitiated air) M2.44・1270 K・1 atm
(質量分率 O₂ 0.258 / N₂ 0.486 / H₂O 0.256) の 2D 試験部に、後向きステップ面のスロット (高さ 0.4 cm, リップ 0.076 cm) から H₂ を
M1・254 K・1 atm で壁平行に噴射。x=35.6 cm の出口で組成 (モル分率)・全温・ピトー圧・M のプロファイル。
上流の断面は 5.1×8.9 cm (ステップ上)、ステップ後は上壁が 9.38→10.5 cm へ直線拡大 (境界層補償)。
計画: [`plans/active/chemistry-finite-rate-h2.md`](../../plans/active/chemistry-finite-rate-h2.md)、仕様: [`methods/chemistry.md`](../../methods/chemistry.md)。

## セットアップ

- `mesh/make_bk_mesh.py` (gmsh API): 平面構造 quad (node 用)。上流ダクト 5 cm + スロット 2 cm + 燃焼器 35.6 cm、68809 ノード、
  下壁第 1 間隔 0.05 mm (y⁺≈30 目標で SST 壁関数)。**品質 PASS** (AR max 38, skew max 0.02; 修正前は transfinite の進行方向不一致で skew 0.98 FAIL)。
- `setup_bk_case.py run_dir [--chem 0|1] [--tci 0|1] [--cfl] [--conv 0|1] [--relax] [--nstep] [--restart bk.h5]`:
  Cantera で入口状態 (空気 ρ 0.2422, U 1784 m/s; H₂ ρ 0.0967, U 1217 m/s) を計算し BC/config/メッシュ変換 (node)/IC を作る。
  BC: `inlet_uniformVelocity` ×2 (組成付き), `outlet_statPress`, `wall_isothermal` (**キーは `Ts`**, 350 K, 壁関数)。
  IC: 一様主流 + y<4 mm の帯を H₂ ジェット状態、壁ノードは u=0・T_wall で整合 (流速込み roe のままだと T=0 → NaN)。
- `restart_conserved.py SRC_res.h5 DST_input.h5`: TP 多成分の同一メッシュ restart (保存量コピー)。
- 参照: `exp_exit_profiles.csv` (Fig.6 組成), `exp_exit_ttot.csv` (Fig.13 全温, T_ref 2380 K) — 図の目視読み取り (±0.03)。
  `compare_exit.py run_dir step out.png` で出口プロファイルを重ねる。

## 計算 run 一覧

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_mixing` | 化学 OFF, 2 次, cfl_pseudo 2, 一様 IC から直接 | 初回: 壁温キー `T` (正しくは `Ts`) で壁 0 K → step 1 NaN。修正後も step 3 でリップ後縁の H₂/空気せん断層で NaN (IC の不連続) | 破棄予定 (起動失敗記録) |
| `run_0002_mixing_start` | 段階起動: 1 次風上, cfl_pseudo 0.5, y<4 mm 帯を H₂ IC, 3000 step (全壁 no-slip 等温) | 完走するが残差が step 1000 以降指数増大。**入口の壁角ノード** (x=−0.02/−0.05 の H₂・空気入口と壁の共有ノード) で P 90 bar: `inlet_uniformVelocity` の u 指定と壁 Dirichlet (u=0, T 固定) が同一ノードで衝突し質量が溜まる (node モードの入口角問題, **未解決**) | 破棄予定 (切り分け記録) |
| `run_0003_mixing_2nd` | run_0002 restart → 2 次, cfl 2 | step 20 で入口上壁角から NaN | 破棄予定 |
| `run_0004_mixing_start_slipup` / `run_0005_mixing_2nd` | 上流ダクト・スロット壁を slip (physID 5) に分離、燃焼器上下壁 no-slip | 起動は安定化したが x=0 の slip→no-slip 切替点 (上壁) が前縁特異点 (P 7–9 bar, T 2200 K) になり 2 次で NaN | 破棄予定 |
| `run_0006_mixing_start_topslip` | **上壁も slip** (no-slip 等温壁はステップ面と燃焼器下壁のみ), 1 次, cfl 0.5, 3000 step | 安定 (残差単調減少) | ref (起動段) |
| `run_0007_mixing_2nd` | run_0006 restart → 2 次, cfl 2, 20000 step (**混合基準**) | 完走・NaN なし。`check_convergence`: `rms_ro` 2e-6〜2e-5 で振動 (リップ後流の非定常, NOT CONVERGED/振動)。出口: H₂ 層が実験より薄い (X_H₂=0.5 が y≈1.0 cm, 実験 1.5 cm), 主流 M 2.55 (上壁 slip で拡大部が過膨張; 実験 2.2–2.44)。`exit_mixing.png` | ref (混合基準) |
| `run_0008_react_jac2` | **反応 ON** (Jachimowski 9 種, jacobianMode 2, coupling 2, cfl_pseudo 1), run_0007 res_20000 から restart, 20000 step | 完走・NaN なし。**自己着火を再現**するが着火点は step 5000→20000 で x 8.9→5.3 cm と上流へ **DRIFTING** (実験の発光観察は 18–25 cm: 早すぎる)。出口 (x=35.6 cm, step 20000): X_H₂O ピーク 0.51 (実験 0.50) だが位置 y=1.48 cm (実験 2.0 cm)、全温ピーク T_tot/2380 K ≈ 1.10 (実験 1.18)。差の主因は混合不足 (上流壁 slip で BL/初期せん断層厚ゼロ) と考えられる。残差は振動プラトー (`rms_ro` 3e-5)。`exit_react_20000.png` | ref (Phase 3 一次結果) |
| `run_0009_start_profile` – `run_0012_profile2_2nd_cfl05` | 入口 BL プロファイル (`inletProfile: 1`, 1/7 乗則 + Crocco–Busemann ρ) + 全壁 no-slip の系列 | 1 次は安定だが 2 次 (cfl 0.5 でも) で入口上壁角 (x=−0.05, y=0.0938) が P 11–15 bar → NaN。**node モードの入口×壁角ノードは未解決の制約** (別 plan 候補: 角ノードで入口状態を優先 / 壁ピンを接線成分のみに) | 破棄予定 (切り分け記録) |
| `run_0014_start_cw` → `run_0015_cw_2nd_cfl05` → `run_0016_mixing_cw` | **入口∩壁角の根治後** (`mesh.nodeInletCornerWall: 1` で変換, 全壁 no-slip 等温 350 K, 入口 BL プロファイル 1/7 乗則 δ=6 mm + Crocco–Busemann ρ): 1 次 cfl 0.5 → 2 次 cfl 0.5 → 2 次 cfl 2 (20000 step, 混合基準) | **全段安定・NaN なし** (角ノード P≈1 atm)。混合: X_H₂=0.5 が y=1.11 cm (旧 1.0, 実験 1.5)、主流 M 2.57 (実験 2.2: 壁関数 BL が薄く拡大部の blockage 不足)。`rms_ro` 2e-6 falling / `rms_roe` プラトー (NOT CONVERGED) | ref (混合基準 v2) |
| `run_0017_react_cw` | run_0016 restart → **反応 ON** (jac 2, coupling 2, cfl_pseudo 1, 20000 step) | **実験に大きく近づいた**: 出口組成 (H₂/N₂/O₂) がほぼ実験に重なり、X_H₂O ピーク 0.51 @ y=1.78 cm (実験 0.50 @ 2.0)、全温ピーク T_tot/2380 K = 1.06 @ 1.70 cm (実験 1.18 @ 1.9)。着火点 x≈9 cm (5000→20000 step で 9.8→8.9 cm、ほぼ定常; 実験の発光 18–25 cm より上流)。`check_convergence` は `rms_ro` 3.6e-6 / `rms_roe` 10 で末尾微増 (振動プラトー, NOT CONVERGED)。`exit_react_20000.png` | ref (**Phase 3 現行ベスト**) |
| `run_0018_react_pasr` | run_0017 restart + **PaSR** (`tci: 1`, C_mix 1, Kolmogorov τ_mix, τ_c=1/max\|J_ss\|) | **消炎** (OH max 2e-7, κ min 0.03): 最速ラジカル時間ベースの τ_c が短すぎ κ が全域で小さい → `tciTauChem: 1` (燃料/酸化剤消費時間) を追加 | ref (PaSR 感度・失敗記録) |
| `run_0019_react_pasr_tauchem1` | PaSR, `tciTauChem: 1` (H₂/O₂ 消費時間), C_mix 1 | κ (火炎内平均) 0.96 → ほぼ No-TCI と同じ (X_H₂O 0.50 @ 1.78 cm, T_tot 1.06, 着火 7.5 cm) | ref (PaSR 感度) |
| `run_0020_react_pasr_cmix001` | PaSR, `tciTauChem: 0`, C_mix 0.01 | κ 平均 0.70: 発熱が下がる (T max 2163 K, X_H₂O 0.47, T_tot 0.91) が**着火位置は 8 cm で動かない** → 早い着火は TCI でなくリップ/壁近傍のモデル (壁関数・等温壁・入口 k/ω) 起因と考えられる | ref (PaSR 感度) |
| `run_0021_react_delta10` | run_0017 restart, 入口 BL δ 6→**10 mm** (報告書「初期境界層厚さ約 1 cm」) | X_H₂=0.5 @ 1.48 cm (実験 1.5), X_H₂O ピーク 0.51 @ 1.82 cm, 着火 7.9 cm, 主流 M 2.51 (変わらず) | ref (δ 感度) |
| `run_0022_react_kx10` | 同 + 入口 k ×10 (ω 固定 → μt/μ ×10), δ 10 mm | 混合が強まり X_H₂=0.5 @ 1.62 cm, X_H₂O ピーク 0.51 @ 2.21 cm (実験 2.0: やや過混合), 全温ピーク 1.07 (幅は実験に近い), **主流 M 2.40** (BL 肥厚で blockage 増), 着火 6.7 cm | ref (入口乱流感度) |
| `run_0023_v2_react_cfl05` → `run_0024_v2_react` | **メッシュ v2** (`make_bk_mesh_v2.py`: 上壁側にもクラスタ, 93k ノード, 品質 PASS), δ 10 mm, run_0017 から cross-mesh 補間 (`interp_field.py`), cfl 0.5 → 1.0 | 壁 y⁺ が全壁で 29–43 (v1 は上壁 ~1300)。出口: X_H₂O ピーク 0.51 @ 1.79 cm, X_H₂=0.5 @ 1.47 cm, 着火 7.8 cm。**主流 M は 2.50 のまま** → 実験の出口 M 2.2 は幅 5.1 cm の側壁 BL による blockage (2D では表現不能) と判断し追わない | ref (v2 基準) |
| `run_0025_v2_delta10_kx3` | **現行ベスト**: メッシュ v2 + δ 10 mm + 入口 k ×3 (run_0022 から補間), 20000 step | 出口 X_H₂O ピーク **0.505 @ 1.96 cm** (実験 0.50 @ 2.0), X_H₂=0.5 @ 1.56 cm (実験 1.5), 全温ピーク 1.08 @ 1.85 cm (実験 1.18 @ 1.9), 着火 5.6 cm (実験 18–25 cm)。`rms_ro` 9.3e-6 falling (NOT CONVERGED)。`exit_react_20000.png` | ref (**Phase 3 ベスト**) |
| `run_0026_regr_oldbin` / `run_0026_regr_newbin` / `run_0027_regr_newbin_rerun` | 2026-09-05 の出口 γ 修正 (plan outlet §2.12) の A/B: run_0025 収束場から 300 step を修正前/修正後バイナリで再実行 + 修正後の再現性 | 修正前後で残差 max 0.57 % 差 (再現ノイズ 1.5e-5) = 壁近傍の亜音速出口列が正しく変わった挙動変更。出口プロファイル (compare_exit) の再確認は未了 | 破棄予定 (回帰記録) |

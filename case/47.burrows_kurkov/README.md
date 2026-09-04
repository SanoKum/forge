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
| `run_0002_mixing_start` | **段階起動**: 1 次風上, cfl_pseudo 0.5, implicitRelax 0.5, y<4 mm 帯を H₂ IC, 3000 step | (記入待ち) | active |
| `run_0003_mixing_2nd` | run_0002 restart → 2 次, cfl_pseudo 2, 20000 step (混合基準) | (記入待ち) | active |

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

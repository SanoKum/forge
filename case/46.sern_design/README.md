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
| `run_0001_smoke_euler_cell` | S4(a) スモーク: S1 テスト 4 の設計 (M_in 2.5, ランプ 15°, カウル 5°/1H, M_c 3.9, f 0.45, p_ext/p_in 0.05, M∞ 6) を cell Euler (slip) で評価し MOC の C_T=0.9666 と照合。段階起動 soft 3000 + 本段 6000 (cfl 4) | (実行中) `metrics.json`, `residual_history.png`, `MESH_QUALITY.txt` (PASS: AR max 120, skew 0.29) | active |

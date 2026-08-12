# case/41.wind_tunnel_design — ①軸対称風洞 (モード F) の設計・帰還検証

[Phase 3 子 plan](../../plans/active/tooling-nozzle-phase3-windtunnel.md) の検証ケース。
条件: Md=4 / Pt=1 MPa / Tt=800 K / rt=10 mm / 円弧 R=2 (Ru=Rd 単一) / 目標 = 中心線
マッハ Bézier (自由 CP 3 + 出口 C2 (Md,0,0))。評価は
`design/.venv-opt/bin/python -m forge_design.evaluate.runner_wt problem_m4.yaml run_NNNN_<slug>`。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_v1_euler` | **v1 チェーン初回 E2E** (逆 MOC 設計壁 [非粘性] → 321×65 メッシュ → forge Euler cell/slip, soft3000→12000 step) | 品質 PASS。**出口軸 M=4.024 (+0.6%)**・ΔM_rms=0.047・ΔM_max=0.119 (3.0%Md)。誤差構造: ①接合部波の軸着地スパイク (x≈1.3cm, +0.12 — 円弧/設計壁の C1 接合の曲率跳び起因)、②減速部の滑らかな欠損 (x≈10–12cm, −0.11 — v2 帰還の主対象)。`achieved_vs_target.csv` / `target_axis_M.csv` / `wall_design.csv` | active (**v1 基準**) |

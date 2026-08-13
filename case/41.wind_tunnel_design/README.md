# case/41.wind_tunnel_design — ①軸対称風洞 (モード F) の設計・帰還検証

[Phase 3 子 plan](../../plans/active/tooling-nozzle-phase3-windtunnel.md) の検証ケース。
条件: Md=4 / Pt=1 MPa / Tt=800 K / rt=10 mm / 円弧 R=2 (Ru=Rd 単一) / 目標 = 中心線
マッハ Bézier (自由 CP 3 + 出口 C2 (Md,0,0))。評価は
`design/.venv-opt/bin/python -m forge_design.evaluate.runner_wt problem_m4.yaml run_NNNN_<slug>`。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_v1_euler` | **v1 チェーン初回 E2E** (逆 MOC 設計壁 [非粘性] → 321×65 メッシュ → forge Euler cell/slip, soft3000→12000 step) | 品質 PASS。**出口軸 M=4.024 (+0.6%)**・ΔM_rms=0.047・ΔM_max=0.119 (3.0%Md)。誤差構造: ①接合部波の軸着地スパイク (x≈1.3cm, +0.12 — 円弧/設計壁の C1 接合の曲率跳び起因)、②減速部の滑らかな欠損 (x≈10–12cm, −0.11 — v2 帰還の主対象)。`achieved_vs_target.csv` / `target_axis_M.csv` / `wall_design.csv` | active (**v1 基準**) |
| `run_0002_v2_loop` | v2 帰還ループ初版 (固定 ω=0.4, 8パス) | 中域は収束するが高 M 末端でリンギング (0.12↔0.16) — 過ゲイン | 破棄予定 (過ゲインの記録) |
| `run_0003_v2_tr` | v2 帰還 trust-region 版 (ω0=0.3) | 13 パスで masked ‖ΔM‖∞ = 0.45% Md。**ただし出口 ε は旧・非軸対称重みで測っており確定値でない** ([plan §9](../../plans/active/tooling-nozzle-phase3-windtunnel.md)) | active (要再測定) |
| `run_0004_v3_cfdmap` | v3 (CFD 場マップ) のつもりだったが、**マップ生成が毎パス失敗し凍結マップにフォールバック**。かつ曲率ブレンド 1.2 が有効 | **無効** — 実験条件が二重に汚染 | 破棄予定 |
| `run_0005_noblend` | fill 修正の切り分け (ブレンドなし) | 切り分けには寄与したが、ε 測定式が未修正 | 破棄予定 |
| `run_0006_v2_clean` / `run_0007_v3_clean` | 「クリーン A/B」のつもりが (i) 呼び出し側 fallback により**ブレンド 1.2 が有効のまま**、(ii) **同名 run の作り直し** (AGENTS ルール違反)、(iii) ε 測定式・収束ゲート未修正 | **無効** — 途中で停止。以後は必ず新連番 (run_0008 以降) | 破棄予定 |

> **注記 (2026-08-14)**: 上記 0004–0007 は実験条件が汚染されているため、いかなる結論の根拠にも使わない。
> 測定器 ($\varepsilon_M$ の軸対称質量流束重み化) と各パスの CFD 健全性ゲートを直してから、
> **新しい連番で 2×2 の A/B** (現行 / B1 / B2 / B1+B2) を投入する ([plan §9.2](../../plans/active/tooling-nozzle-phase3-windtunnel.md))。

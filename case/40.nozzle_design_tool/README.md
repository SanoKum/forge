# case/40.nozzle_design_tool — ノズル設計ツール (forge_design) の検証ケース

`design/forge_design` (ノズル設計ツール — [親計画](../../plans/active/tooling-nozzle-design-tool.md) /
[Phase 0 子 plan](../../plans/active/tooling-nozzle-phase0-foundation.md)) の評価パイプライン検証用ケース。
問題定義 YAML (`problem_bell_smoke.yaml`) から「区分構成壁 → 構造化 quad msh4.1 →
`convertGmshToForge` → 品質ゲート → 準 1D 等エントロピー IC → forge (SLAU+SST+陰解法) →
収束判定 → 推力メトリクス」を 1 コマンドで実行する:

```bash
PYTHONPATH=design python3 -m forge_design.evaluate.runner \
    case/40.nozzle_design_tool/problem_bell_smoke.yaml \
    case/40.nozzle_design_tool/run_NNNN_<slug> [--steps N] [--prepare-only]
```

条件はいずれも ③ベル (Pt 4 MPa / Tt 1500 K / 背圧 20 kPa / rt 10 mm / ε 9 /
θa 20° / L/rt 7)、メッシュ 221×65 (14,080 cells)、品質 VERDICT PASS (AR max 19 / skew max 0.42)。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_bell_smoke` | E2E 初回スモーク。SST の roK/roOmega を IC 未設定 (既定 0) | step 2–6 で rms_roK/roOmega に NaN → EOS 床洗浄後プラトー (既知指紋の再現)。パイプライン自体は完走し η=0.981 | 破棄予定 (NaN 指紋の記録) |
| `run_0002_bell_smoke_sstic` | IC に roK/roOmega (k=1, ω=18000) を追加、壁第一セル y+~70 (壁関数) | **NaN なし**。残差 1.1–1.7 桁低下後プラトー (`CONVERGENCE_VERDICT.txt`: NOT CONVERGED/stalled — cell モード atomicAdd ノイズ床水準)。**推力は 4000→12000 step で 0.002% ドリフト** (F=1961 N, η=0.9790, ṁ=1.299 kg/s, `metrics.json`)。`residual_history.png` あり | active (Phase 0 E2E 基準) |

注記: 残差プラトーの深掘り (真の定常収束品質) は Phase 2 (Rao 照合) で扱う。η の妥当性
(0.98 前後はベルノズルとして自然なオーダー) も Phase 2 の照合が正式判定。

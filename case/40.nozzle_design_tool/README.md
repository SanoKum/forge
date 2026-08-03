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
| `run_0003_bell_L6` | dv 応答確認: L/rt=7→6 (`problem_bell_L6.yaml`) | η=0.9708, F=1944.7 N (短縮 → 発散損失増で η 低下 ✓)。品質 PASS・NaN なし・プラトー同傾向 | active (dv 感度の記録) |
| `run_0004_bell_L9` | dv 応答確認: L/rt=7→9 (`problem_bell_L9.yaml`) | η=0.9837, F=1970.5 N (延長 → η 単調増 ✓)。ṁ は 3 run で 0.03% 一定 (スロート同一の整合) | active (dv 感度の記録) |
| `run_0005_bell_smoke_node` | node 再計算の初回 (冷間 IC + implicit + 細壁 1e-3) | step 2 から全列 NaN。後の切り分けで真因 = 壁第一セル過細 + レシピ不一致と確定 → [調査 plan](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) §2.5 | 破棄予定 (発散の記録) |
| `run_0006_bell_node_warm` | **node 初収束**: 実証レシピ (壁 5e-3・warm start from run_0002・explicit RK3 cfl0.1・1次・katoLaunder) | **VERDICT: PASS / ALL PASS** (全列 3.0–4.6 桁低下)。η=0.9809, F=1964.8 N (cell run_0002: η=0.9790 と 0.2% 差) | active (node 基準) |
| `run_0007_bell_node_2nd_imp` → `_2nd_expl` | 2 次化・陰解法化の試行 (run_0006 収束場から) | **いずれも step 4–10 で発散** — node×SST の陰解法と 2 次精度は未解決のソルバ課題 (case/29 の実績も 1 次 explicit のみ)。plan §2.5 参照 | 破棄予定 (発散の記録) |
| `run_0008_bell_node_runner` | runner 実装後の E2E 一気通貫確認 (`--warm-from`, 24000 step) | **VERDICT: PASS / ALL PASS** (3.4–5.1 桁・全列 falling)。η=0.9817, F=1966.4 N — node レシピが runner 既定で完結 | active (node runner 基準) |

**dv 感度まとめ (η–L トレードオフの応答確認)**: L/rt = 6 / 7 / 9 → η_CF = 0.9708 / 0.9790 /
0.9837。物理的に正しい単調応答で、パイプラインが最適化の評価器として機能することを確認
(Phase 2 のパレート軸そのもの)。

**node モードの状態 (2026-08-03 更新)**: 実証レシピ (1 次・explicit cfl0.1・warm start・壁
5e-3) を runner の node 既定に実装し **node で収束可能** (run_0006 ALL PASS)。ただし
2 次精度・陰解法は未解決 ([plan](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)
の残課題) のため、速度・精度が要る評価では cell を併用する。run_0001–0004 は壁 1e-3/2e-3
世代のメッシュ (現 problem yaml は 5e-3) である点に注意。

注記: 残差プラトーの深掘り (真の定常収束品質) は Phase 2 (Rao 照合) で扱う。η の妥当性
(0.98 前後はベルノズルとして自然なオーダー) も Phase 2 の照合が正式判定。

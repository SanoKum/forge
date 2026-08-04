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
| `run_0008_bell_node_runner` | runner 実装後の E2E 一気通貫確認 (`--warm-from`, 24000 step) | **VERDICT: PASS / ALL PASS** (3.4–5.1 桁・全列 falling)。η=0.9817, F=1966.4 N — node レシピが runner 既定で完結 | active (旧 node 1次 explicit 基準) |
| `run_0009_node_sst_imp_repro` | 課題1再現試行: 1次+陰解法 cfl_pseudo 2 (bndFirstOrder あり), 300 step | 「step 4 発散」は再現せず (300 step 健全)。旧発散記録は 2 次との交絡と判明 | 破棄予定 (切り分け記録) |
| `run_0010_node_2nd_repro` | 課題2再現試行: 2次+explicit cfl0.15 (bndFirstOrder あり), 300 step | 300 step 健全。反実仮想 (bndFirstOrder 除去, scratchpad) は step 8–10 で発散 = **課題2の因子は bndFirstOrder 欠落と確定** | 破棄予定 (切り分け記録) |
| `run_0011_node_2nd_imp` | 2次+陰解法 cfl2 の短尺確認 (300 step) | 300 step 健全 (全列 1.4–2.0 桁低下) — ただし長尺は下記 0012 で発散 | 破棄予定 (切り分け記録) |
| `run_0012_node_2nd_imp_full` | 2次+陰解法 cfl2 の 12000 step (軸修正**前**) | **step 10612 で発散** (roK がベル部近軸で e-fold ~2000 step の指数成長 → roOmega inf)。課題1の実態 = 遅発性 | 破棄予定 (発散の記録) |
| `run_0013_node_imp_long` | 1次+陰解法 cfl2 の 12000 step (軸修正**前**) | **step 7875 で発散** (同モード)。= 課題1は convMethod 非依存、軸行真空化×SST k シートが真因 ([plan §2.6](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | 破棄予定 (発散の記録) |
| `run_0014_node_2nd_expl_long` | 2次+explicit cfl0.15 の 12000 step (bndFirstOrder あり) | **完走・全列 2.5–2.8 桁 falling** = 課題2は bndFirstOrder で解決済みの確証 | active (課題2解消の証拠) |
| `run_0015_node_imp_axisdir` | **nodeAxisDirichlet=1** + 1次+陰解法 cfl2, 12000 step | 完走。rms_ro 1.4e-9 / roK 2.5e-7 (floor)、roK 成長モード消滅。軸床ピン 0 ノード。η=0.9816 (1次 explicit 0.9817 と一致) | active (軸修正の効き) |
| `run_0016_node_2nd_imp_axisdir` | **最終目標: 2次+陰解法 cfl2 + 軸 Dirichlet** (runner E2E, 12000 step) | **VERDICT: ALL PASS**・quasisteady ALL STEADY。η=0.9907, F=1984.4 N, ṁ=1.2923。**12000 step ≈ 20 秒** (explicit レシピ ~7 分の ~20 倍) | active (**node 生産基準**) |
| `run_0017_cell_kato` | cell 対照: run_0002 + katoLaunder (warm from run_0002) | η=0.9790 = run_0002 と同一 → **node−cell の η 差 +1.2% は katoLaunder では説明されない** (2次化の node/cell 離散差) | active (η 帰属の対照) |
| `run_0018_bell_L6_node` | L スイープ node 取直し: L/rt=6 (壁 5e-3 修正版 yaml, warm from run_0016) | 完走 (4.0–6.1 桁低下, roOmega は床で振動 2.9e-3↔5.6e-3)。**η=0.9835** | active (dv 感度 node) |
| `run_0019_bell_L9_node` | 同 L/rt=9。cross-geometry warm は cfl4 直投入で発散 → **段階起動** (soft explicit 3000 step → 陰解法 12000 step) | **VERDICT: PASS**。**η=0.9957** | active (dv 感度 node) |
| `run_0020_node_runner_cfl4` | 2次+陰解法 **cfl_pseudo 4** + bndFirstOrder の E2E (旧既定) | **VERDICT: PASS / ALL PASS** (3.2–4.6 桁)。η=0.9907 (cfl2 の run_0016 と一致) | active (bndFirstOrder あり世代の記録) |
| `run_0021_node_2nd_imp_nobfo` | **bndFirstOrder 撤去の主検証**: 全域 2次+陰解法 cfl4+軸 Dirichlet | 12000 step NaN なし。rms_ro ~1e-7 プラトー (入口/収縮部の微小リミットサイクル, ro 相対 ~1e-5)・SST 4.1–4.2 桁。quasisteady **ALL STEADY**。η=0.9896 | active (全域 2 次の主検証) |
| `run_0022_node_runner_nobfo` | **現 runner 既定** (全域 2次+陰解法 cfl4, bndFirstOrder なし) の E2E 基準 | run_0021 と同等 (η=0.9896, ALL STEADY, プラトー同性状) | active (**runner 既定基準**) |
| `run_0023_bell_L6_nobfo` | L/rt=6 を現既定で取り直し (warm from run_0021) | **VERDICT: PASS (converged)**。η=0.9822 | active (dv 感度 node 現行) |
| `run_0024_bell_L9_nobfo` | L/rt=9 を現既定で取り直し (段階起動: soft explicit 3000→陰解法 12000) | rms_ro 9.2e-9 (実質床, flat 判定)・ALL STEADY。η=0.9942 | active (dv 感度 node 現行) |
| `run_0025_node_wallsmooth` | 壁列エントロピー市松の種切り分け (warm start 壁 T 平滑化) | res_0 は市松 5K に清浄化されるが **12000 step 後に run_0022 とビット同一の 206K 市松が再生成** = 市松は定常解固有 ([plan §2.8](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)) | active (壁市松の証拠) |
| `run_0026_node_su2axis` | **SU2 流軸対称 (`axisymMethod: 1`)** の一次検証 (全域2次+陰解法 cfl4, nodeAxisDirichlet **なし**, warm from run_0021 収束場) | **軸行が真に健全** (床0・kシート消滅・軸線滑らか)。12000 step 完走・η=0.9904・ALL STEADY。ただし implicit は喉部近軸 limit cycle で rms_ro ~9e-5 プラトー (explicit は 3e-7 到達 = 空間健全) → [SU2化 plan](../../plans/active/axisymmetric-su2-source-formulation.md) | active (axisymMethod 1 基準) |
| `run_0029_node_adiabfix` | **断熱壁熱流束リーク修正の反実仮想** (現既定構成, 修正のみ, warm from run_0022) | **壁 T 市松 mean 206→43.7K (−79%)・壁温 1631→1195K** (旧値は Tt 超えの非物理、SU2 は 1428K)。η=0.9884 (−0.12%)・ALL STEADY。以後の node 断熱壁 run はこの物理が既定 | active (**断熱修正基準**) |
| `run_0028_su2_sst` | **SU2 v8.5 クロスチェック** (同一メッシュを gmsh で .su2 化, 同一 BC, RANS-SST/ROE/MUSCL/Euler implicit+FGMRES, 手順書テンプレ) | 15000 iter 完走 + 固定 CFL10 継続 5000 iter。**SU2 も rms[Rho] 10^-5.1〜-5.6 でプラトー** (CFL_ADAPT 共振ではない=固定 CFL でも同水準)。軸は完全健全 (T zigzag 0.5K)。**壁 T 市松なし** (mean ~4K vs forge node 206K, 出口リップ 1 点 244K のみ) = 壁市松は forge node 固有の確証。η_CF=0.9573 (F=1917.6N, ṁ=1.2811) — forge より 2-3% 低いが乱流入口 BC・壁 y+~19-27 low-Re が未整合の粗比較 (真値照合は Phase 2) | active (SU2 対照) |
| `run_0027_node_rfloor` | **r 床 (`axisRFloor: 3.0e-4`, ユーザ提案)** の検証 (method 0, nodeAxisDirichlet **なし**, 全域2次+陰解法 cfl4) | 軸健全 (床0)・η=0.9906・rms_ro ~1e-7 プラトー。床帯縁に有界の k 帯 (~230; scratchpad soak +12000 step でビット不変の平衡)。素朴な面床は閉性破れで発散 → 離散閉性面積 (A_closure) で解決 | active (axisRFloor 基準) |

**軸処理 3 方式の比較 (全域 2 次+陰解法 cfl4, bndFirstOrder なし)**:

| 方式 | 軸ノード | row1 k 帯 | 収束 (rms_ro) | η_CF |
| --- | --- | --- | --- | --- |
| `nodeAxisDirichlet: 1` (現 runner 既定) | 解かずミラー | 3.4 | 1e-9 床 (PASS 可) | 0.9896 |
| `axisymMethod: 1` (SU2 流) | 真に解く | **1.35** | implicit ~1e-4 プラトー (explicit 3e-7) | 0.9904 |
| `axisRFloor: 3e-4` (r 床+閉性補正) | 解く (環状近似) | 230 (有界) | ~1e-7 プラトー | 0.9906 |


**dv 感度まとめ (η–L トレードオフの応答確認)**:
L/rt = 6 / 7 / 9 → η_CF は cell (2026-08-03 前半, 壁 1e-3 メッシュ) 0.9708 / 0.9790 / 0.9837、
**node 全域 2次+陰解法 (現既定, run_0023 / 0021 / 0024, 壁 5e-3)** 0.9822 / 0.9896 / 0.9942
(bndFirstOrder あり世代 run_0018/0016/0019 は 0.9835 / 0.9907 / 0.9957)。いずれも物理的に
正しい単調応答で、node−cell はほぼ一定の +1.1〜1.2% オフセット (トレンド保存)。オフセットの
帰属: katoLaunder ではない (run_0017)・軸修正でもない (修正前 run_0012 も 0.991) → node 側
2 次化の離散差。真値の判定は Phase 2 (Rao 照合・格子収束) で行う。

**node モードの状態 (2026-08-04 更新)**: 残っていた 2 課題は解決した
([plan §2.6–2.7](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md))。
2 次発散・陰解法の遅発性発散は**いずれも node 軸対称の軸行真空化が種** (ベル部で軸 CV が
radial 圧力平衡からデカップルし EOS 床まで過膨張 → 偽せん断 → SST k シート) — ソルバ新機能
**`nodeAxisDirichlet: 1`** (軸ノードを radial 隣接からの対称 Dirichlet に置換) で根治。
当初レシピに入れていた `bndFirstOrder: 1` (境界隣接 1 次化) は軸病変のマスキングであり、
壁境界層の精度を損なうため**撤去** (ユーザ指示 2026-08-04)。**現 runner 既定 =
全域 2 次 (convMethod 1) + 陰解法 (blockDPLUR, cfl_pseudo 4) + nodeAxisDirichlet**。
1 評価 12000 step ≈ 20 秒 (旧 explicit レシピ比 ~20 倍)。残差は形状により入口/収縮部の
微小リミットサイクルで ~1e-7 プラトーになるため (cell 全域 2 次と同格)、**計量は
`check_quasisteady.py` STEADY を確認して使う**。warm start は引き続き必須で、
**同メッシュ node 場からが最良**。異形状/cell 場からの warm は cfl4 直投入で初手発散する
ことがあり、その場合は段階起動 (soft explicit 3000 step → 陰解法、run_0019/0024 の手順)。
run_0001–0004 は壁 1e-3/2e-3 世代のメッシュ (現 problem yaml は全て 5e-3 に統一済み)。

注記: 残差プラトーの深掘り (真の定常収束品質) は Phase 2 (Rao 照合) で扱う。η の妥当性
(0.98 前後はベルノズルとして自然なオーダー) も Phase 2 の照合が正式判定。

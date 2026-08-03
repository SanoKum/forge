# セッションプロンプト: node ノズルの陰解法・2 次精度発散の修正

以下を新セッションの最初のプロンプトとして使う。

---

node (median-dual) のノズル計算で残っているソルバ課題 2 点を修正してほしい。
基準文書は [`plans/active/boundary-node-nozzle-wall-outlet-stability.md`](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)
(調査マトリクス・実証レシピ・残課題の全記録) — まずこれを読むこと。

## 背景 (確定済み事実 — 再調査不要)

- node ノズル (case/40 ③ベル, 傾斜壁+超音速流出) は**実証レシピでは完全収束する**:
  壁第一セル `wall_first_frac: 5e-3`・収束場からの warm start・explicit RK3 cfl0.1・
  convMethod 0 (1 次)・`bndFirstOrder/nodeWallDirichlet/axisCentroidShift: 1`・katoLaunder。
  成功例 = `case/40.nozzle_design_tool/run_0006_bell_node_warm` / `run_0008_bell_node_runner`
  (いずれも VERDICT PASS / ALL PASS、η_CF ≈ 0.981 で cell と 0.2–0.3% 整合)。
- リグレッションではない (case/29 `run_0038_node_sst_wc` verbatim が現行バイナリで健全)。
- **課題 1: node × SST × 陰解法 (timeIntegration 11, blockDPLUR)** — 収束済み node 場からの
  引き継ぎでも step 4 で roOmega から発散 (cfl_pseudo 2)。case/29 でも implicit 成功は
  層流のみ (run_0035 cfl2–5)。発散種は壁隣接第一内点。
- **課題 2: node × 2 次精度 (convMethod 1)** — 同条件で step 7–10 に発散 (cfl 0.3/0.15 とも、
  explicit でも)。case/29 の node SST も全て 1 次だった。発散種は同じく壁隣接。
- 動機: explicit cfl0.1 は 1 評価 ~7 分と遅く、設計ツール (plans/active/tooling-nozzle-*) の
  最適化ループには陰解法+2 次が必須。ユーザは node 主体の方針。

## 再現手順 (数分/run)

```bash
# ベース (PASS するレシピ, warm start 必須):
PYTHONPATH=design python3 -m forge_design.evaluate.runner \
  case/40.nozzle_design_tool/problem_bell_smoke.yaml <新run_dir> \
  --warm-from case/40.nozzle_design_tool/run_0006_bell_node_warm/res_12000.h5 --prepare-only
# ↑で生成された solverConfig.yaml を編集して切り分け:
#   課題1: timeIntegration 3→11, blockDPLUR 0→1, cfl 0.1→2.0  → step~4 で roOmega NaN
#   課題2: convMethod 0→1 (explicit のまま cfl 0.15)          → step~7-10 で NaN
# 実行は必ず tools/run_case.sh <絶対パス> (hook 準拠)。detectNaN: 1 は生成済み。
# NaN 位置の特定: res_nan_*.h5 の VALUE を nozzle.h5 の /CELLS/centCoords と
#   /BCONDS/<id>/iCells で引く (本セッションで使った手法は plan §2 参照)。
```

## 修正の手がかり (plan §3 より)

- 課題 1: blockDPLUR の k/ω 行 — 壁 Dirichlet/壁関数ノード・出口ノードの Jacobian 整合。
  関連: `plans/accepted/turbulence-node-sst-wallfunction.md`,
  `time_integration-scalar-dplur` 系, `discretization-node-wall-implicit-dirichlet.md`
  (DeleteValsRowi)。
- 課題 2: 壁近傍ノードの勾配再構成 (Green-Gauss / gradLSQ) と `bndFirstOrder` の適用範囲。
  case/29 `run_0040_node_sst_lsq` (LSQ, 1 次) は健全 → 再構成 (MUSCL/2nd) 側の壁隣接パスが濃厚。

## ルール (AGENTS 準拠の要点)

- ソルバ変更前に plan の `status`/方針を更新 (§3 を実装計画に具体化)。コード変更後は
  **full rebuild** ([stale-build-struct-layout-trap] — 差分ビルドは CUDA obj 取りこぼしあり)。
- 検証: ①課題別の再現 run が発散しなくなること、②実証レシピ run (run_0006 相当) が
  ビット劣化しないこと、③最終 = node SST **2 次+陰解法** 12000 step で
  `check_convergence.py` PASS、η_CF が cell (run_0002: 0.9790) と ~1% 一致。
  ④既存 node 検証ケースの非退行 (case/26 平板 node 500 step 健全 — 手順は plan §2.5)。
- run は `case/40.nozzle_design_tool/run_NNNN_<slug>` (連番は 0009 から)、README の
  run 一覧表を同期。「収束した」報告には VERDICT 貼付必須。
- 完了時: 修正が入ったら forge_design runner の node 既定を 2 次+陰解法に引き上げ、
  L スイープ (run_0003/0004 相当) を node で取り直して cell との整合を README に記録。

# node-centered の `node*` オプション削減 (整合セットを既定化し、旧規約系フラグを撤去)

## メタ

- **area**: `architecture`
- **status**: `in-progress` (ユーザ指示 2026-08-16「node なんちゃらオプションが多すぎる、追って削減」)
- **related_docs**: [`methods/discretization.md`](../../methods/discretization.md), [`methods/axisymmetric/implementation.md`](../../methods/axisymmetric/implementation.md)
- **related_plans**: [`architecture-node-centroid-value-position.md`](architecture-node-centroid-value-position.md) (整合セットの根拠 = case/43)
- **created**: `2026-08-16`

## 1. 現状の `mesh.node*` / 関連フラグと処遇案

| フラグ | 現既定 (node) | 意味 | 処遇 |
|---|---|---|---|
| `nodeValueAtNode` | **1** | 値の位置 = ノード座標 (r̄ は `rEff` で分離) | 既定固定 → **フラグ撤去** (converter が最初からノード座標+双対重心別データセットを書く形へ = 親 plan Step 2/3) |
| `nodeReconEdgeMidpoint` | **1** | 2 次再構成の目標 = エッジ中点 (SU2 流) | 既定固定 → **撤去** (旧「双対面重心」目標は削除) |
| `nodeAxisUrDirichlet` | **1** (軸対称) | 軸 u_r=0 の状態+残差+Jacobian 三点セット | 既定固定 → **撤去** (壁 no-slip と同様に「node × 軸対称なら常に」) |
| `gradLSQ` | **2** (node) | 勾配 LSQ (内部隣接のみ) | node は 2 固定、cell は従来 0 → 「node では GG 禁止」を明文化。フラグは cell 用に残す |
| `nodeAxisDirichlet` | 0 | 軸ノードを第一内点コピーで置換 (対症) | **撤去** (保存を破る・軸を 1 次化する。三点セットで不要) |
| `axisCentroidShift` (converter) | 1 | centCoords ← 双対重心 | converter が「centCoords=ノード座標 + `/CELLS/dualCentroid`」を書くようにし **撤去** (solver 側スワップも不要に) |
| `nodeMidpointFx` | 0 | fx を 0.5 に固定 | 値=ノードでは幾何 fx が自動的に中点相当 → **撤去** |
| `nodeWallDirichlet` | 1 | 壁 no-slip の三点セット | 既定固定 → **撤去** (常に ON) |
| `axisRFloor` | 0 | r 床 (軸行真空化の対症) | 軸行真空化は値位置整合で消える → **撤去候補** (case/40 で再確認後) |
| `hoopAreaFromClosure` | 0 | hoop 面積 = 離散閉性 | 生産で無害・自由流厳密 → **既定 1 → 固定 → 撤去** |
| `bndFirstOrder` | 0 | 境界 1 次化 | 既に使用禁止 → **撤去** (別 plan) |
| `nodeWallStressEdgeKernel`, `nodeWallViscGradFlux`, `nodeOmegaWfDirichlet` 等 | — | 壁処理の実験系 | 本 plan では触らず、壁 plan で整理 |

到達形: `mesh` ブロックの node 関連は `discretization: node` と `isAxisymmetric` だけ。

## 2. 手順

1. **統合前 (今)**: 既定を整合セットに切替済み (solverConfig.cpp)。旧規約は明示 0 で再現可能 (回帰用)。case/41 生産チェーン (Euler axis-Mach / NS) で既定の妥当性を確認 → OK なら `feature/median-dual-3d` へ統合。
2. **統合後 1**: 生産 run 系列 (case/41, 40, 29 等) を既定で再走し、旧 Dirichlet 系列との差を run 表に記録。
3. **統合後 2**: converter を「centCoords=ノード座標 + dualCentroid 別データセット」に変更 (h5 フォーマット +1 データセット、旧 h5 は solver 側スワップで互換)。
4. **統合後 3**: 上表「撤去」フラグを削除、docs/plans/README を同期。cell モードはビット不変を全段で確認 (case/08, 24, 36 cell)。

## 3. 変更ログ

- `2026-08-16` — 起案。node 既定を整合セット (値=ノード / エッジ中点 / 軸 u_r 三点 / LSQ) に切替 (`61e09a3b`〜)。
- `2026-08-16` (夜) — **手順 4 の第 1 弾を実装**: `nodeAxisDirichlet` (enforceAxisMirror / zeroAxisAllResiduals / axis_rep /
  5 行 decouple ごと)・`nodeMidpointFx` を削除、`nodeValueAtNode` / `nodeReconEdgeMidpoint` / `nodeAxisUrDirichlet` を
  常時 ON にしてキー撤去、node の `gradLSQ` は 2 固定 (他値はエラー)。旧キーを書いた config は起動時に明示エラー。
  runner から旧キー参照を除去。残: `axisCentroidShift` (converter 書式変更 = 手順 3)、`axisRFloor`、`hoopAreaFromClosure`
  の既定化、`nodeWallDirichlet` の固定化。


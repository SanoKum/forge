# node 壁サーフェス出力の修正 (primal 面トポロジ退避 + 壁関数 twall モデル値化)

## メタ

- **area**: `architecture` (出力 / 壁診断)
- **status**: `done`
- **related_docs**:
  - `docs/diffusion/implementation.md` 「node-centered の壁摩擦応力 twall」節
- **related_plans**:
  - `design/accepted/diffusion-node-wall-viscous-distance.md` (node 壁 twall 集約)
  - `design/accepted/architecture-node-centroid-value-position.md` (node viz は primal セル + Center='Node')
- **created**: `2026-06-27`
- **owner**: `CFD Dev`

## 1. 目的

node モードの壁サーフェス出力 `res_wall_*` が ParaView で**何も表示されない**問題と、壁せん断 `twall`
出力が**壁関数メッシュで物理値の ~十数倍**になる不整合を直す。前者は XDMF トポロジが壊れていたため、
後者は `twall` が解像 μ_total·∂_n u で `utau`/`y+` (対数則) と不整合だったため。

## 2. スコープ

- **やる**: (a) 壁トポロジを XDMF 妥当に (cell 2D 線含む)、(b) node 壁を primal 境界面で接続し
  Center='Node' で描く (体積の `vizCONNE` と同思想、一度退避し h5 永続化)、(c) 壁関数 active 時に `twall`
  を modeled τ_w=ρu_τ² に再スケール、(d) 出力専用カーネルの改名。
- **やらない**: cell 壁関数 (wt=1) の twall 再スケール (AddTauWall は node 専用。別途)。

## 3. 関連 docs と前提

- 体積の node viz は `mesh.vizCONNE` (primal セル接続) + Center='Node' (architecture-node-centroid-value-position)。
  壁はその境界版が無かった。
- `Tau_Wall=ρu_τ²` は `ransWallFunction` (computeWallFrictionSST) が viscousFlux より前に算出し、AddTauWall が
  運動量に課す値 (diffusion-node-wall-viscous-distance)。再スケールはこれを流用 (新規計算なし)。

## 4. 設計方針

### 4.1 壁トポロジ (XDMF Mixed) の修正
旧 `output.cpp` は要素タイプを nn=2/3/4 しか分岐せず、**node の壁半割面 (1 ノード, nn=1) で未初期化型
(ゴミ 24747) を書き ParaView 不可視**、2D 線 (nn=2) も Polyline をノード数カウント無しで書き不正だった。
可変ノード型 Polyvertex(1)/Polyline(2) は `[type, nNodes, node...]` とカウントを入れ、固定型
Triangle(4)/Quad(5) は `[type, node...]` とする。

### 4.2 node 壁を primal 境界面で接続 (一度だけ退避)
`replacePrimalWithDual` で境界も半割面 (1 ノード) に置換され面を作れないため、置換前に各 bcond の
**primal 境界面の global ノードid列を `bcond.vizBfaceNodes` に退避**し、h5 `/BCONDS/<id>/vizBface{Nodes,Sizes}`
に永続化 (readMesh で復元)。出力時 node モードはこれで壁ノードを線/面に接続し、値はノード単位に並べ替え
**Center='Node'** で書く。cell モードは従来どおり面単位 (Center='Cell')。

### 4.3 壁関数時の twall モデル値化
出力専用カーネルで twall を書いた後、`Tau_Wall[W] > 0` (=壁関数 active のゲート) のとき向きは解像 traction・
大きさを τ_w=ρu_τ² に再スケール。`Tau_Wall<0` (壁解像/層流, nullptr 渡し) では解像値=真の τ_w なので無補正。
→ どちらのモードでも twall が物理的な壁せん断で `utau`/`y+` と整合し、正しい Cf が出る。

### 4.4 改名
出力専用 (res/y+ 不触) を明確にするため `viscousFlux_wall_node_d` → `wallStressForOutput_node_d`。

## 5. 実装ステップ

1. `mesh/mesh.hpp`: bcond に `vizBfaceNodes`。
2. `mesh/gmshReader.hpp`: `replacePrimalWithDual` で退避、`writeInputH5` で h5 保存。
3. `mesh/mesh.cpp`: `readMesh` で復元。
4. `output/output.cpp`: 壁 CONNE/値/xmf を node 面・Center='Node' 対応、XDMF エンコーダ修正。
5. `cuda_forge/viscousFlux_d.cu`: カーネル改名 + `Tau_Wall` 引数 + 再スケール (ガードは AddTauWall と同形)。
6. `docs/diffusion/implementation.md` 同節を更新。

## 6. 検証

- 再変換した case/36 node メッシュ: `res_wall_4` が 1224 本の Polyline・Center='Node'・1226 ノード値・全
  局所 index < npoints (VALID)、Ps/twall が空間的に整合。
- twall 再スケール: case/36 node SST で `twall/(ρu_τ²)` = 1.00–1.02 (修正前 17–20)。

## 7. 影響範囲

- h5 スキーマ追加 (`/BCONDS/<id>/vizBface*`) → **既存 node メッシュ h5 は `convertGmshToForge` 再変換が必要**
  (無くても cell は不変、node 壁は旧挙動)。cell モード・cell の res_wall は挙動不変 (XDMF は spec 準拠化)。
- `docs/diffusion/implementation.md` 更新。

## 8. 完了条件

- [x] 壁トポロジ XDMF 妥当化 / node primal 面接続 / twall モデル値化 / 改名
- [x] case/36 で検証
- [x] docs 更新 / `status: done` / accepted 配置 / `design/README.md` 同期

## 9. 変更ログ

- `2026-06-27` — 壁トポロジ + node primal 面 viz + h5 スキーマ (commit dcfb3e1)、twall 再スケール + 改名
  (commit 1e4fc46)。**既存 node メッシュは再変換が必要**。

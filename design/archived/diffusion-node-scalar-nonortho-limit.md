# node-centered flux の dcc に node 座標を使う (双対 CV 重心由来の見かけ非直交を排し SST omega 爆発を根治)

## メタ

- **area**: `diffusion`
- **status**: `superseded` <!-- → architecture-node-centroid-value-position.md (centCoords=node + ghostless で根治。ncx 別配列は廃止) -->
- **related_docs**:
  - `docs/diffusion/implementation.md`
  - `docs/discretization/implementation.md` (node 近壁・axisCentroidShift)
- **related_plans**:
  - 同根: [`diffusion-node-wall-viscous-distance.md`](../accepted/diffusion-node-wall-viscous-distance.md) (壁面の dcc 退化。本 plan は内部面)
  - 関連: [`discretization-node-wall-implicit-dirichlet.md`](../accepted/discretization-node-wall-implicit-dirichlet.md)
- **created**: `2026-06-22`
- **owner**: (未割当)

## 1. 目的

node-centered で SST RANS が omega 爆発で起動できない (case/36: roOmega が step3 で 1e22→inf)。**真因はリミッタ不足でなく、flux の `dcc` (セル間距離) に「双対 CV 重心」(`axisCentroidShift` でノードからシフト) を使っていること**。双対面はノードエッジに直交するので、dcc に **node 座標**を使えば直交 (over-relaxed 係数 `delta=area/cosθ` の `1/cosθ`=1)。ところが双対重心を使うと近壁の半割 CV で重心がノードからずれ、見かけの非直交 `1/cosθ` が**実測 max 1.46e6** に発散し、μ_t×近壁 omega 勾配を増幅して omega 拡散 flux を爆発させる。本 plan は **node mode の flux dcc を node 座標に切替**え、これを根治する。

> ユーザー指摘 (2026-06-22): 「変にリミット入れる前に、そんな非直交面が存在しないはずなのに存在しているのが問題」「cell では存在せず node で生じている」。検証で正しさを確認 (node 座標で `1/cosθ`=1.00、双対重心で 1.46e6)。当初の over-relaxed 係数リミッタ案は band-aid として棄却。

## 2. スコープ

- **やる**: node mode の flux 幾何で **dcc に node 座標**を使う。`ncx/ncy/ncz` device 配列 (node mode=ノード座標, cell mode=centCoords) を追加し、scalar 拡散 + 運動量粘性 (内部面) の dcc に適用。
- **やらない**:
  - 壁面 (half-face/ghost) の dcc 退化 (別 plan [`diffusion-node-wall-viscous-distance.md`](../accepted/diffusion-node-wall-viscous-distance.md))。
  - 対流 flux の fx/dcc、勾配 (calcGradient) の node 座標化 (残発散の追加対策として後続検討)。
  - `centCoords` (双対重心) 自体は volume/軸対称 r 重み/source 用に**温存** (axisCentroidShift の意図は保持)。
  - cell モード (ncx=ccx でビット不変)。

## 3. 関連 docs と前提

- 内部面拡散/粘性 flux は over-relaxed `delta = |dcc|·|S|²/(dcc·S) = |S|/cosθ`。直交 (dcc∥S) で `delta=|S|`、非直交で発散。**median-dual の双対面はノードエッジに直交**するので、dcc=ノード間ベクトルなら厳密直交。
- `axisCentroidShift=1` (既定) は軸上 R=0 CV のゼロ回転体積回避のため `centCoords`=双対重心にする。これは volume/軸対称には必要だが、**flux の dcc に使うとノードずれ (近壁 max 0.0075) が見かけ非直交を生む**。`axisCentroidShift=0` 全体フリップは volume/勾配と不整合で step0 NaN (棄却)。→ flux dcc だけ node 座標に分離するのが正解。
- 切り分け: `scalarDiffusion=0` で omega 爆発停止 (拡散項が犯人)、近壁 `delta/area` max=1.46e6 (centCoords) vs 1.00 (node 座標)。

## 4. 設計方針

`ncx/ncy/ncz` (cell 配列, geom) を追加: node mode の内部 CV (`ic<nCells`, `iNodes={ic}`) は `nodes[iNodes[0]].coords`、ghost と cell mode は `centCoords`。flux kernel の dcc 計算で `ccx`→`ncx` に差し替える (kernel 内ロジック不変、wrapper で渡す配列を変えるだけ)。cell mode は `ncx=ccx` で完全ビット不変。

## 5. 実装ステップ

1. `variables.hpp` の cell 配列リストに `ncx/ncy/ncz` 追加 (allocVariables が device 確保)。
2. `variables.cpp` `setStructuralVariables_d`: `ncx` を充填 (node mode=ノード座標 / それ以外=centCoords)、device へ upload、free。
3. `cuda_forge/scalarTransport_d.cu` wrapper: scalar 拡散 kernel に `ccx`→`ncx` を渡す。
4. `cuda_forge/viscousFlux_d.cu` wrapper: 内部面 `viscousFlux_d` に `ccx`→`ncx`。

## 6. 検証

- **ビルド**: native (struct 不変、`variables.hpp` 変更で関係 obj 再コンパイル)。
- **回帰 (cell)**: case/36 cell SST 50 step で baseline (git stash) と比較、差は atomicAdd 自己揺れと同オーダー (0.0271 vs 0.0272) = 無影響を確認済。
- **主検証 (case/36 node SST, implicit, nodeWallDirichlet=1)**: omega 爆発が消え収束に向かう。
- **判定基準**: ① cell ビット不変、② node SST が omega 爆発しない、③ omega/mean flow が低下。

## 7. 影響範囲

- 触るファイル: `variables.{hpp,cpp}`, `cuda_forge/scalarTransport_d.cu`, `cuda_forge/viscousFlux_d.cu`。
- ドキュメント: `docs/diffusion/implementation.md`。

## 8. 完了条件

- [x] `docs/diffusion/implementation.md` に node 座標 dcc を追記
- [ ] 残発散 (case/36 step1657) の追加対策 (対流 fx/dcc・勾配の node 座標化、または near-convergence 安定化) を別 plan/§9 で判断
- [ ] [`design/README.md`](../README.md) を更新
- [ ] `status: done` 化 (残発散の整理後)

## 9. 変更ログ

- `2026-06-22` — 初稿は over-relaxed 係数リミッタ案だったが、ユーザー指摘「非直交面は人工物」を受け再検証。**node 座標で `1/cosθ`=1.00、双対重心で 1.46e6** を確認し、真因=flux dcc が双対重心を使うこと、と判明。リミッタ案を棄却し **node 座標 dcc** に方針変更。
- `2026-06-22` — **実装 (scalar 拡散 + 運動量粘性) + 検証**。`ncx/ncy/ncz` 追加、両 flux の dcc を node 座標へ。
  - **効果 (case/36 node SST, implicit, cfl_pseudo=1.0)**: 無修正 step3 爆発 → scalar のみ step87 → **scalar+viscous で step1657 まで完走**。omega **2.1 桁単調低下** (6.7e5→5.7e3, ~step1500)、mean flow 収束 (rms_roUy 3.9→0.30, rms_ro 0.022→0.0016)。
  - **cell 回帰**: baseline 比 atomicAdd 自己揺れ同オーダー (0.0271/0.0272) = ビット不変相当。
  - **残**: step1657 で最終発散 (1656 step 健全収束の後)。対流 flux の fx/dcc・勾配 (calcGradient) も centCoords を使うため同系の near-convergence 残不安定の可能性。wall viscous half-face (ghost) の dcc 退化は別 plan。追加対策は次段。

# node-centered の centCoords を「値の位置 (ノード座標)」に統一し、双対重心/軸半径を分離する

## メタ

- **area**: `architecture`
- **status**: `in-progress`
- **related_docs**:
  - `methods/discretization.md`
  - `methods/diffusion.md`
- **related_plans**:
  - 置換: [`diffusion-node-scalar-nonortho-limit.md`](../archived/diffusion-node-scalar-nonortho-limit.md) (ncx 別配列でなく本 plan の centCoords 統一で解決)
  - 依存: [`diffusion-node-wall-viscous-distance.md`](../accepted/diffusion-node-wall-viscous-distance.md) (壁境界の退化 dcc を弱形式化。本 plan の前提)
  - 関連: [`discretization-node-wall-implicit-dirichlet.md`](../accepted/discretization-node-wall-implicit-dirichlet.md), [`architecture-axisym-axis-singularity.md`](../accepted/architecture-axisym-axis-singularity.md)
- **created**: `2026-06-22`
- **owner**: (未割当)

## 1. 目的

node-centered で未知量の位置は**ノード**だが、現状 `CELLS/centCoords` には `axisCentroidShift=1` で**双対 CV 面積加重重心**(ノードからずれた点)が入る。このため対流再構成・勾配・拡散/粘性 dcc・fx が「値が重心にある」前提で動き、近壁の半割双対 CV で重心ズレ(max ~0.0075)が**見かけの非直交 `1/cosθ` max 1.46e6** を生み、SST omega 拡散を爆発させる (case/36 step3)。`axisCentroidShift` は名前(軸半径の話)と影響範囲(centCoords 全体・非軸対称にも波及)がズレている。

本 plan は **`centCoords` を「値の位置」= node モードではノード座標に統一**し、cell/node で同一処理にする。双対 CV 重心は別量に分離し、軸対称の r 重み/source だけが使う。これにより前出の見かけ非直交が消え (検証済: 内部面 dcc を node 座標にすると case/36 SST が step3→1657 へ改善)、`ncx` のような node 専用分岐も不要になる。

**追記 (2026-08-11)**: 双対 CV 重心/体積の**シューレース float32 桁落ち**を修正した
(`gmshReader.hpp`, A 相対座標+double 化 — y+≈1 スリバー CV で重心が ~100 μm 級に狂い CV 外へ
出て `wall_dist`/interp を破壊していた。[boundary plan §2.9](boundary-node-nozzle-wall-outlet-stability.md) 参照)。
これで「双対重心の値が正確」にはなったが、本 plan の主題 (**値の位置=ノード座標への統一**と
双対重心の分離) は未着手のまま残る。

## 2. スコープ

- **やる**:
  - `replacePrimalWithDual`: `centCoords` = ノード座標 (axisCentroidShift 分岐を撤廃)。
  - 双対 CV 重心 (面積加重) を**別配列** (`cvCentroid` / 軸半径 `r_eff`) に保存。
  - 軸対称の r 重み (`volume *= r_cell`) と軸ソースを `r_eff`(=双対重心半径 or `revolved_vol/A_planar`)に切替。`centCoords.y` を半径代用にしない。
  - **node モードの ghost を完全撤廃** (`nCells_all = nCells`)。境界 (壁/流入出/slip/axis) は**半割面の弱形式 flux** = bvar 境界値 + **内部隣接ノードから作った `∇φ·S`** で評価。ghost mirror も `pc` 距離も使わない。← centCoords=node を成立させる必須前提であり、ユーザー要求 (ghost 概念の撤廃)。
  - `axisCentroidShift` フラグの撤去/整理。
- **やらない**:
  - cell モードの挙動 (centCoords は従来どおり cell 中心・ghost あり、ビット不変)。
  - 内部スキーム自体の変更 (幾何の供給先を変えるだけ)。

## 3. 関連 docs と前提

- 現状: `axisCentroidShift=1` → `centCoords=dualCentroid`。`variables.cpp` の軸対称 r 重みは `r_cell=ccy[ic]`(=重心 y)、`r_face=pcy[ip]`(=面重心 y)。**r_face は面重心なので正しく、問題は r_cell が双対重心 y を使う点**。
- `centCoords=node` で `ccy=node.y=0`(軸上)→ `r_cell=0`→回転体積 0(=`axisCentroidShift` が避けていた問題)。よって `r_eff=revolved_vol/A_planar=双対重心 y` を別途持つ。
- **境界退化**: `centCoords=node` で壁ノードが壁面に乗り、ghost mirror `cc_ghost=cc+2((pc-cc)·n)n≈cc`→`dcc≈0`→ゼロ割 (検証: axisShift=0 の NaN 132/132 が壁近傍)。`axisShift=1` は重心を内側にずらし**たまたま**回避していた。→ 境界法線は一方向コンパクト/弱形式が必須。

## 4. 設計方針

- **値の位置 = ノード**: `centCoords[i] = nodes[i].coords` (node mode)。全 flux/勾配/再構成/補間が同一コードで node 座標を使う。内部双対面はノードエッジに直交するので over-relaxed `delta=area/cosθ` が `1/cosθ=1`(直交)になり、見かけ非直交が消える。
- **双対重心/軸半径の分離**: `cvCentroid`(または軸半径 `r_eff`)を `/CELLS` 等に保存。軸対称 r 重み・軸ソースのみ参照。`r_eff = revolved_volume / A_planar`(=面積加重半径)で定義可能。
- **ghost 撤廃・境界は内部隣接ノードの弱形式**: node モードでは境界ノードが境界上に乗るため、ghost mirror (`cc_ghost=cc+2((pc-cc)·n)n`) も `pc` 距離も法線方向に退化する (検証: axisShift=0 で壁 dcc≈0 → NaN 132/132 が壁)。**`pc` でなく内部隣接ノードへ行く**のが正しい: 境界半割面の法線項は**内部隣接ノードから構築したセル勾配 `∇φ` を用いた `∇φ·S`**(bvar が境界値で勾配を閉包) で評価し、ghost を一切作らない。
  - **隣接ノード情報は保有を確認済** (case/36: 全壁ノード 1226 個が内部双対面で内部隣接ノードを min2/mean3/max3 個持ち、内部隣接ゼロは 0 個)。`map_cell_planes`/`map_plane_cells` の内部双対面が各境界ノード→内部ノードの接続を与える。
  - 運動量は既存 `nodeWallViscGradFlux` が同型 (壁法線を `∇u·S` 化)。スカラ omega/k・流入出/slip も同方式へ。= [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) Phase 2 (全 ghost 撤廃) の完遂。

### 4.1 壁 omega BC の定式化 (2026-07-19 ユーザー承認・確定)

変更ログ 2026-06-22 の未解決項目 3 (「centCoords=node で壁ノード wall_dist=0 → SST 壁 omega
`ω=6ν/(β·y_w²)` が overflow」) の方針を確定する: **ghostless plan で実証済みの方式を踏襲**する。

- 壁 omega の実効 y は**第一オフ壁ノード距離** (`compute_wall_y_eff` の **MIN** 集約、凹コーナーは
  2-ring 探索) で評価し、壁ノードの `omega[ic]`/`roOmega[ic]` を**保存変数ごと直接 Dirichlet ピン**、
  `applySSTPointImplicit` で壁 ω 行を decouple (dω=0)、`zeroWallMomentumResidual` で res_roOmega=0。
- 根拠: この組み合わせで node SST が cell 基準と一致済み (Cf/Schlichting 0.89–0.93, peak mut/μ 202
  vs cell 199。[`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md)
  変更ログ 2026-06-22 の 3 修正 0f6d53d/2b19d1d/8b60f2c)。centCoords=node 化で壁ノードが厳密に
  壁面上 (wall_dist=0) に乗っても、この定式化は y を第一オフ壁距離から取るため破綻しない。
- **やらないこと**: 壁ノード y=0 のまま kSmall フロアで凌ぐ (overflow の温床)、k=0 Dirichlet
  (automatic WT では k はノイマンが正・過去に誤適用で悪化)、Menter 原式の y を双対重心シフトで
  ごまかす (axisCentroidShift の masking 再導入になる)。

## 5. 実装ステップ

段階導入 (各段で cell ビット不変 + node 検証を確認しつつ):

1. **境界 ghostless 化 (前提)**: node モードで ghost 生成を停止し (`nCells_all=nCells`)、全境界 (壁/流入出/slip/axis) を半割面弱形式 (bvar + 内部隣接ノードの `∇φ·S`) に。= node-boundary-ghostless Phase 2。これを先に入れないと centCoords=node が壁で NaN する。
2. **convert**: `mesh/gmshReader.hpp` `replacePrimalWithDual` で `centCoords=node 座標`、`dualCentroid` は別フィールド/データセットへ。
3. **mesh I/O**: `mesh.cpp` で `cvCentroid`/`r_eff` 読み込み・device 配列化 (`variables`)。
4. **軸対称 r 重み**: `variables.cpp` の `r_cell` を `r_eff` に。`r_face` は現状維持 (面重心)。軸ソース (`axisymmetricSource_d.cu`) の半径も `r_eff` に。
5. **後片付け**: `axisCentroidShift` フラグ削除、docs 更新。

> 注: ghostless 化 (Step 1) が最大の作業。完了すれば centCoords=node (Step 2) で内部面は自動直交化し、ncx のような分岐は不要 (cell/node 同一処理)。

## 6. 検証

- **回帰 (cell)**: 全 cell ケースがビット不変 (centCoords は cell では不変)。case/08, case/24, case/36 cell。
- **node 非軸対称 (case/36)**: SST が omega 爆発せず収束に向かう (step3 爆発の解消)。近壁 dP/dy ノイズ低減。`check_convergence.py` VERDICT 記録。
- **node 軸対称 (case/29)**: 層流・SST が回帰 (近軸固着なし、推力/中心線 Mach が従来同等)。`r_eff` で R=0 ゼロ体積が再発しないこと。
- **判定基準**: ① cell ビット不変、② case/36 node SST が爆発しない、③ case/29 node が回帰維持 (軸半径対策が効く)。

## 7. 影響範囲

- 触る: `mesh/gmshReader.hpp`, `mesh/mesh.cpp`, `variables.{hpp,cpp}`, `cuda_forge/axisymmetricSource_d.cu`, スカラ境界拡散 (`scalarTransport_d.cu`/新カーネル), `input/solverConfig.{hpp,cpp}` (axisCentroidShift 撤去)。
- ドキュメント: `methods/discretization.md`, `methods/diffusion.md`, `methods/index.md`。

## 8. 完了条件

- [ ] docs 更新 (centCoords=値の位置 / 軸半径分離 / 境界弱形式)
- [ ] 実装・検証完了 (§6)
- [ ] [`plans/README.md`](../README.md) 更新、`diffusion-node-scalar-nonortho-limit.md` を superseded 化
- [ ] `status: done`、§9 に case/36・case/29 の before/after

## 9. 変更ログ

- `2026-06-22` — 起案。ユーザー提案 (centCoords=値の位置, 軸半径は専用量, ncx 別定義を作らず cell/node 同一処理) を受け、`diffusion-node-scalar-nonortho-limit.md` の ncx 案を置換。検証: 内部面 dcc を node 座標化で case/36 SST step3→1657 改善 (ncx 実験)、`centCoords=node`(axisCentroidShift=0) は壁退化 dcc で NaN (132/132 壁) を確認 → 境界弱形式が前提と判明。
- `2026-06-22` — **Stage 1+2 実装 (commit 5c22025, 6a7a21d)**。node モードの境界スカラ (k/ω) 拡散を ghost dcc でなく **∇φ·S 弱形式** (calcGradient が bvar+内部隣接で閉包) に。Stage1=壁のみ→Stage2=全境界 (ghost 検出 `ic>=nCells`、periodic は実 partner で除外)。gated node-only で cell ビット不変。case/29 node SST 回帰 OK (rms_ro 7.95e-5→3.84e-6 収束)。隣接ノード保有確認 (全壁ノード内部隣接 2-3 個)。
- `2026-06-22` — **`centCoords=node` の連鎖ブロッカーを段階的に発見** (`axisCentroidShift=1` が複数問題を同時に masking していた)。case/36 で axisShift=0 を試すたびに次の masked 問題が露出:
  1. 内部面 dcc 非直交 (1.46e6) → centCoords=node で解消。
  2. 境界スカラ ghost dcc 退化 → Stage1+2 で解消。
  3. **(未解決) SST 壁 omega BC `ω=6ν/(β·y_w²)` が wall_dist=0 (壁ノードは壁面上) で overflow→inf** ([ransBoundary_d.cu:41,50](../../solver_density_cuda/cuda_forge/ransBoundary_d.cu#L41))。axisShift=1 は重心を内側にずらし wall_dist>0 にしてこれも隠していた。**node-centered では壁 omega を「壁ノード(y=0)」でなく「第一内部ノードの距離」で定式化する必要**があり、SST 壁処理の設計判断を要する (誤ると乱流が狂う) ため次段でユーザーと方針決定。
  4. **(未解決) エネルギー(roe)・運動量の境界拡散 dcc 退化**。層流(SST なし)でも centCoords=node で step0 roe=NaN を確認 (roUx/roUy は健全)。`nodeWallViscGradFlux` は**壁**の運動量/熱流束を ∇u·S/∇T·S 化するが、**非壁境界(inlet/outlet/slip)の運動量・エネルギー粘性**は ghost dcc のまま退化する。
  5. (未着手) 軸対称 r_eff 分離 (case/29 軸対称回帰用)、ghost 生成停止 (`nCells_all=nCells`)、convert で centCoords=node 本実装。

  **連鎖の総括**: `axisCentroidShift=1`(重心を内側にずらす)は、**全境界・全方程式(スカラ/運動量/エネルギー)の境界拡散 dcc 退化**と**SST 壁 omega の y=0 overflow**を**一括で masking**していた。centCoords=node を成立させるには、これら全部を ghostless 化(= node-boundary-ghostless Phase 2 の全方程式版)+ omega BC 再定式化 + r_eff 分離、が必要。**Stage1+2 (境界スカラ) のみ完了**。残りは大規模で omega BC は設計判断 (node-centered SST 壁の y 定義) を要するため、ユーザーと方針確定後に段階実装する。
  **現状到達点**: Stage1+2 commit 済 (5c22025, 6a7a21d)。centCoords=node 本フリップは残りの ghostless 化 + omega BC 定式化が揃ってから。
- `2026-07-19` — **omega BC 方針をユーザー承認で確定** (§4.1 追加): ghostless plan 実証済みの
  「wall_y_eff MIN + omega 直接ピン + point-implicit decouple」を踏襲。これで未解決 3 項目のうち
  設計判断待ちは解消し、残りは実装のみ (§5 の Step 1 ghostless 完遂 → Step 2–5)。
- `2026-08-16` — **試作 `mesh.nodeValueAtNode: 1` を実装 (solver 側スワップ, worktree ブランチ `feature/node-axis-dof`)**。
  変換器は触らず `mesh.cpp readMesh` で (ゴースト生成前に) `centCoords ← node 座標`、双対重心 y を `mesh::rEff`
  に退避し `variables.cpp` の r 重み `r_cell` だけが `rEff` を使う (= §4 の「値の位置=ノード、軸半径は専用量」)。
  境界半割面の退化 (ノードと鏡映ゴーストが面上) は `calcStructualVariables_d` の fx を `0/0→0.5` にガード。
  **粘性 dcc 退化は未対応** (Step 1 ghostless が前提) なので Euler 限定・既定 0。
  根拠 = [case/43.node_axis_dof](../../case/43.node_axis_dof/README.md): 離散連続式の演算子テストで
  軸半 CV の質量残差 e=+0.83 (GG)/+0.5 (fx=0.5 or LSQ) が h 不変の O(1) 不整合であること、
  `nodeValueAtNode: 1` で丸め床に落ちることを確認 (`optest/`)。ノズル Euler (case/41 生産モデル) で軸 M 欠損
  0.033→≤0.007、cell 参照との近軸差 0.0067→0.0031。**残課題**: 出口∩軸コーナー 1 ノード (M −0.11)、
  粘性 (Step 1)、implicit の軸 roUy 行 (SU2 型の状態・RHS・Jacobian 一体射影, discretization.md §7.1 の
  過去失敗を踏まえ corner を独立検証)。設計チェーンは並行して軸直読→偶関数外挿に切替
  (`forge_design/evaluate/axis_extract.py`)。
- `2026-08-16` (追補) — **コーナー欠陥の真因確定と修正 / NS 非対応の原因特定**。
  (1) 前項の「出口∩軸コーナー」は**境界ゴーストの回転体積潰れ → dt_local 崩壊 → CV 凍結**だった
  (値位置=ノードで鏡映距離 0 → ゴーストが r=0 → r 床 1e-20 → `setDT` の `dx=vol/|S|`~1e-20 → CFL 1e13)。
  `variables.cpp` の r 重みでゴーストに所有 CV の r̄ を使う修正で解消 (case/43 run_0010: 凍結 0・残差 2.7 桁)。
  **ghost 撤廃 (Step 1) が未了だと値位置=ノードは境界でこの種の退化を必ず起こす**ことの実例。
  (2) NS は本フラグ未対応。原因は**壁近傍 2 次再構成**で、軸 DOF・CFL・段階起動のいずれとも独立
  (case/43 run_0009〜0017)。y+≈1 スリバーでのみ破綻し粗い壁格子では完走する。→ Step 1 (ghostless) に加え
  **「値=点値」への統一 (勾配・リミッタの点値評価)** を明示的な作業項目として §4 に追加する必要がある。


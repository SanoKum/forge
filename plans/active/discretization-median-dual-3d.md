# 3D median-dual (M4) — 3次元双対メッシュ生成と periodic 双対面対応

## メタ

- **area**: `architecture` (discretization 横断)
- **status**: `in-progress`
- **related_docs**:
  - [`methods/discretization.md`](../../methods/discretization.md) (median-dual 理論・2D 実装解説)
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md)
- **related_plans**:
  - [`discretization-median-dual.md`](discretization-median-dual.md) (親計画。M1–M3 実装済、本計画が M4 を独立詳細化)
  - [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) (node 境界弱形式)
  - [`turbulence-iddes-sst.md`](turbulence-iddes-sst.md) (3D node DDES の最終目的)
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 1. 目的

node-centered (median-dual) モードを **3D 混合要素 (tet/prism/hex) + spanwise periodic** に対応させ、
3D node DDES/LES を可能にする。現状 [`buildMedianDual()`](../../solver_density_cuda/mesh/gmshReader.hpp#L1305)
は `if (is3D)` で即 bail し (`dualBuilt=false`)、node モードは 2D 平面のみ。periodic は 2D dual builder にも未実装
(双対面を周期境界で対応付ける処理が無い)。本計画で 3D 双対生成と periodic 双対面対応を実装する。

**ブロッカー確認 (2026-06-28)**: SU2 三者比較で「市松はスキーム (RHS 散逸) で決まる」と判明し
([notes/investigations/backstep-lowmach-checkerboard.md](../../notes/investigations/backstep-lowmach-checkerboard.md))、
3D node では roe (内在散逸) を検証軸にする方針。その前提として 3D node 自体が必要。

## 2. スコープ

- **やる**:
  - 3D primal (tet/prism/hex/pyramid) から **一意エッジ列**を構築 (現状 3D の `planes` は primal *面*でエッジ列が無い)。
  - **3D 双対面** (エッジ {A,B} ごとのポリゴン: 各 incident cell で `edge-mid M / 面重心 F1,F2 / cell 重心 G` の
    四角パッチ集約) の面ベクトル・面積・重心。
  - **3D 双対体積** (各ノード A の sub-polyhedron: A・edge-mid・面重心・cell 重心で構成する多面体、符号付き四面体分割)。
  - **3D 境界半割面** (境界面 tri/quad を median 分割しノードへ分配、マルチマーカ corner 処理は 2D と同型)。
  - **periodic 双対面対応**: 周期境界をまたぐエッジの双対面を partner 側ノード CV と接続 (spanwise 解像乱流に必須)。
  - 検証: 双対体積総和 = primal 体積、ノード closure `Σ dS + 半割面 = 0`、**free-stream 保存** (3D 一定流で残差 0)。
  - I/O: `/DUAL` の 3D 書き出し拡張、solver 接続マップ、`output.cpp` の `Center='Node'` 3D 化、`setInitial`/probe の 3D node 化。
- **やらない**:
  - 非平面四角面の厳密分割 (Newell 法で近似、歪み大は警告)。
  - polyhedral (ieleType=0) primal の双対化 (tet/prism/hex/pyramid のみ)。
  - rotation periodic の 3D 双対 (Cartesian translation periodic を先行。rotation は後続)。
  - 移動メッシュ・AMR。

## 3. 関連 docs と前提

- 2D median-dual の現在仕様: [`methods/discretization.md`](../../methods/discretization.md)。本計画着手前に同 doc へ
  「3D 双対生成」節 (エッジ列・双対面ポリゴン・多面体体積・periodic 対応) を追記する (AGENTS 開発フロー)。
- 既存 2D 実装 [`gmshReader.hpp:1305-1590`](../../solver_density_cuda/mesh/gmshReader.hpp#L1305) のデータ構造
  (`dualVolume[nN]`, `dualCentroid[nN*3]`, `dualFaceCells[nDF*2]`, `dualFaceVect[nDF*3]`, `dualFaceArea[nDF]`,
  `dualFaceCent[nDF*3]`, 境界 `dualBcond*`/`dualBnode*`) を 3D に踏襲。2D は `nDF=nPlanes` だったが 3D は
  `nDF=一意エッジ数`。
- periodic: [`mesh::setPeriodicPartner()`](../../solver_density_cuda/mesh/mesh.cpp#L488) の node 対応付け
  (Cartesian dx/dy/dz) を双対面構築に接続する。

## 4. 設計方針

### 4.1 一意エッジ列の構築
3D 要素 (tet/prism/hex/pyramid) の局所エッジ定義テーブルから、全セルのエッジを `(min(n0,n1),max(n0,n1))` で
重複除去し `edges[]` を作る。各エッジに **incident cells** と、各 cell 内で**エッジを含む 2 面**を記録
(`edge -> {cell, faceA, faceB}` リスト)。primal 面 (`planes`) は面重心算出に流用。

### 4.2 3D 双対面 (エッジ {A,B})
各 incident cell の寄与 = 四角パッチ `(M, F1, G, F2)` (M=edge midpoint, F1/F2=エッジを含む 2 面の重心, G=cell 重心)。
面ベクトルは Newell 法 (パッチ頂点の外積総和の 1/2)、A→B 向きに統一して全 incident cell 分を加算 = 双対面ベクトル。
面積=ノルム、重心=パッチ面積加重。2D の「midpoint→centroid 区間」を 3D ポリゴンに一般化したもの。

### 4.3 3D 双対体積 (ノード A)
各 incident cell で、A 周りの sub-polyhedron を `A・(A 接続エッジの midpoint)・(A を含む面の重心)・G` を頂点とする
四面体群に分割し符号付き体積を加算。全 cell 分の総和が `dualVolume[A]`。重心 `dualCentroid` は四面体重心の体積加重。

### 4.4 3D 境界半割面
境界面 (tri/quad) を median 分割: 各構成ノード N に、`N・(N 接続エッジ midpoint)・(面重心)` のサブポリゴンの
外向き面ベクトルを分配。マルチマーカ corner (壁∩出口∩periodic) は 2D と同型に各 incident bcond へ 1 枚ずつ
(`halfByOwner`)。閉性集計 `bnodeAccum` は所有非依存で全半割面。

### 4.5 periodic 双対面対応 (核心)
周期境界面上のエッジは、partner 面側に同一形状のエッジが存在する。`setPeriodicPartner` のノード対応 (dx,dy,dz 平行移動)
でエッジ A↔A', B↔B' を対応付け、**周期エッジの双対面を内部双対面として両側 incident cell を繋ぐ** (境界半割面に
しない)。これにより spanwise の解像乱流が周期方向に連続する。実装: periodic 面を通常境界として半割面化せず、
partner セルを `dualFaceCells` の他端に設定 (周期 ghost ではなく直接接続、cell モードの `periodicPartner` と同思想)。

### 4.6 検証量
- 双対体積総和 = primal 体積総和 (relErr < 1e-6)。
- ノード closure: `Σ_f (符号付き双対面ベクトル) + Σ 境界半割面 = 0` (periodic 内部面は両側で相殺)。
- **free-stream 保存**: 3D 一定流 (ρ,u,p 一様) で全ノード残差が機械精度。非直交・周期で特に重要。
- 負体積 0。

## 5. 実装ステップ

1. **docs 先行**: [`methods/discretization.md`](../../methods/discretization.md) に「3D 双対生成」節を追記 (§4 の設計)。
2. **エッジ列**: `gmshReader.hpp` に要素別ローカルエッジテーブル + 一意エッジ構築 + `edge→{cell,face,face}` マップ (§4.1)。
3. **3D 双対面/体積/重心**: `buildMedianDual()` の `is3D` 分岐を実装 (§4.2–4.3)。Newell 法・四面体体積ヘルパ追加。
4. **3D 境界半割面**: §4.4 を実装 (2D の `halfByOwner`/`hcentByOwner` を 3D 面サブポリゴンへ一般化)。
5. **periodic 双対面**: §4.5。`bconds` の periodic マーカを検出し partner ノード対応で双対面を内部接続。
   `mesh.cpp` の periodic 接続と整合。
6. **検証チェック**: §4.6 を `buildMedianDual()` 末尾に追加 (体積・closure・負体積)、free-stream は別途 run。
7. **I/O・solver**: `/DUAL` 3D 書き出し、`output.cpp` Center=Node 3D、`setInitial`/`point_probes` の 3D node 化、
   `replacePrimalWithDual` の 3D 対応確認。
8. **検証 run** (§6)。

## 6. 検証

- **単体/ビルド**: native ビルド。`convertGmshToForge` を 3D tet/hex メッシュへ適用し体積総和・closure・負体積を確認。
- **検証ケース** (易→難):
  1. **3D 直方体一定流** (tet/hex): free-stream 保存 (残差機械精度) — 双対面ベクトル整合の最小確認。
  2. **3D periodic チャネル/box**: spanwise periodic で一定流保存 + 周期方向連続性。
  3. `case/05.sod_shock_tube` の 3D hex 押し出し: shock 位置・一定流が cell/解析解と整合。
  4. **`case/18.backstep` 3D (`backstep.h5`, spanwise periodic)**: roe+RANS で cell/2D node と整合 → roe+LES →
     KEEP 3D。DDES (`turbulence-iddes-sst`) は本計画完了後。
- **判定基準**: cell モード完全不変。3D node の free-stream 保存 (機械精度)・体積/closure 合格・periodic 連続。
  `verify-node-and-cell-both`・AGENTS 収束/準定常手順に従う。

## 7. 影響範囲

- ファイル: `mesh/gmshReader.hpp` (中核), `mesh/mesh.cpp` (periodic 接続), `output/output.cpp`,
  `input/setInitial.hpp`, `probe/point_probes.cpp`, `methods/discretization.md`。
- 既存: cell モード・2D node は不変 (`is3D` 分岐の新規追加なので 2D 経路ビット不変)。
- docs: `methods/discretization.md` + `methods/index.md` 整合。

## 8. 完了条件

- [ ] `methods/discretization.md` に 3D 双対生成節を追記
- [ ] 実装・検証完了 (§6: free-stream 保存・体積/closure・periodic 連続・cell 不変)
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/` へ移動、親 plan の M4 を done 反映、[`plans/README.md`](../README.md) 同期

## 9. 変更ログ

- `2026-06-28` — 初稿 (計画)。3D node DDES 移行のため M4 を親 plan から独立詳細化。2D 実装の構造
  ([gmshReader.hpp:1305-1590](../../solver_density_cuda/mesh/gmshReader.hpp#L1305)) を踏襲し 3D 双対面/体積/
  境界半割面/periodic 対応を設計。SU2 検証で「3D node は roe を散逸軸にする」方針を前提化。
- `2026-06-28` — **ステップ 1–4・6 (検証チェック) を実装・検証** (ブランチ `feature/median-dual-3d`)。
  - docs: [`methods/discretization.md`](../../methods/discretization.md) §2.5「3D 双対生成 (M4)」を追記
    (エッジ列・双対面・多面体体積・境界半割面・periodic・検証量)。
  - 実装: [`gmshReader.hpp`](../../solver_density_cuda/mesh/gmshReader.hpp) に `buildMedianDual3D()` を追加し
    `buildMedianDual()` の `is3D` 分岐から呼ぶ。§4.1 一意エッジ列 (要素局所面テーブル `nodesOrderPlanes` の
    周回ノード対から `edge→{cell,2面}` を導出; 別エッジテーブル不要)、§4.2 Newell 法による双対面パッチ
    `(M,F1,G,F2)` を A→B 統一集約、§4.3 四面体群 `(node,M,F,G)` による双対体積・面積加重重心、§4.4 境界面
    (tri/quad) の median 分割半割面、§4.6 体積総和/closure/負体積チェックを 3D 化。コードは要素タイプ非依存
    (tet/prism/hex/pyramid 共通パス)。
  - **付随バグ修正**: primal 体積計算が tet を `eleType.name=="tet"` で判定していたが登録名は `"tetra"` で
    不一致 → **全 tet セルの primal 体積が 0** だった既存バグを修正 (cell モードの tet メッシュも影響)。
  - 検証 (native `convertGmshToForge`, node モード):
    - **hex** (`case/35` tg.msh, 32³ box): 体積 dual=248.05=primal=248.05 (relErr 3.5e-8)、closure
      normalized=2.5e-6、境界半割面積=primal 面積、負体積 0。
    - **tet** (生成した単位 box, 197 tetra): 体積 dual=1=primal=1 (relErr 3.1e-7)、closure 9.7e-8。
    - prism/pyramid は同一コードパス (tri/quad 面プリミティブ・table 駆動) で未実行 run のみ。
  - **未着手**: §4.5 periodic 双対面 (現状は周期境界も通常境界として半割面化)、free-stream 保存 run
    (検証ケース 1、solver/IO の 3D node 化=ステップ 7 が前提)、ステップ 7 (I/O・solver 接続)、検証ケース 2–4。

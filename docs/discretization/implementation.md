# 離散化レイアウト (cell / node) 実装解説

本ドキュメントは、cell-centered / node-centered (median dual) 両対応の実装方針と
ソース対応を記す。理論は [theory.md](theory.md) を参照。

## 1. 設計の核: discretization mesh 抽象

forge の GPU カーネル群 (対流・勾配・粘性・block-DPLUR) は、消費する幾何が次の 4 群に限られる:

1. **CV 中心量**: `volume[ic]`, `ccx/ccy/ccz[ic]` ([variables.cpp](../../solver_density_cuda/variables.cpp))
2. **面量**: `sx/sy/sz[ip]`, `ss[ip]`, `pcx/pcy/pcz[ip]`, `fx[ip]`, `dcc[ip]`
3. **面 → CV**: `map_plane_cells_d[2*ip]` ([mesh.cpp](../../solver_density_cuda/mesh/mesh.cpp))
4. **CV → 面 CSR**: `map_cell_planes_index_d` / `map_cell_planes_d`

これらは「CV と面」の抽象だけで、primal セル形状を参照しない。`fx`/`dcc` も CV 重心と
面重心からランタイム再計算される ([calcStructualVariables_d.cu](../../solver_density_cuda/cuda_forge/calcStructualVariables_d.cu))。

→ **median-dual を前処理で構築し、同じ `nCells/nPlanes/map_plane_cells/CSR` レイアウトに
「CV=ノード, 面=エッジ双対面」として詰めれば、内部カーネルは無変更で動く。**

`mesh` 構造体 ([mesh.hpp](../../solver_density_cuda/mesh/mesh.hpp)) は既にこの汎用抽象になっている。
よって solver はこの抽象だけを消費し、**供給側 (mesh を作る前処理) を 2 実装に分ける**:

- **cell モード**: 現行 gmshReader 経路。`/CELLS` を CV とする。**現状維持**。
- **node モード**: 新 `dualMeshBuilder` が双対量を HDF5 `/DUAL` グループに書き、
  `readMesh` がそこから `nCells := nNodes`, `nCells_ghst := 0` として読む。

切替は `solverConfig` の `discretization: cell | node` (既定 `cell`)。
**`nCells` を `nCV` にリネームしない** (差分肥大回避)。node モードでは `nCells == nNodes` と読み替える
(その旨を [mesh.hpp](../../solver_density_cuda/mesh/mesh.hpp) の doc コメントに明記)。

## 2. 双対メッシュ生成 (前処理)

配置: `convertGmshToForge` 前処理 ([gmshReader.hpp](../../solver_density_cuda/mesh/gmshReader.hpp))。
既存の primal 幾何構築 (セル重心・面重心・向き合わせ) と、**既に構築済みだが未使用の
`node.iCells` / `node.iPlanes`** ([mesh.hpp:19-20](../../solver_density_cuda/mesh/mesh.hpp)) を起点に再利用する。

生成物:
- **双対体積** $V_p$: 各 primal セルを「セル重心–面重心–エッジ中点–ノード」で細分し各ノードへ集計。
  `Σ双対体積 == Σprimal体積` を assert。
- **双対面** (primal エッジ ↔ 1 枚): `surfVect = Σ パッチ面積ベクトル`, `surfArea = |Σ|`,
  `centCoords = 面積加重重心`。`plane_cells = {n0, n1}`、法線は n0→n1 向きに統一。
- **境界半割双対面**: 各 primal 境界面を構成ノードへ分配し `bnode_face_{sx,sy,sz,ss,pc}` を格納。
  `Σ半割面積 == primal face 面積` を assert。
- **閉性チェック**: 各内部ノードで `Σ(双対面ベクトル, 符号付き)` の max ノルムを print。
- 出力: HDF5 `/DUAL` グループ (可視化用 primal は `/MESH /CELLS` に残す二層構造)。

## 3. カーネル別の対応状況

| カテゴリ | ソース | node 対応 |
| --- | --- | --- |
| 対流 (Roe/SLAU/AUSM) MUSCL | [convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) | 無変更で成立 |
| Green–Gauss 勾配 | [calcGradient_d.cu](../../solver_density_cuda/cuda_forge/calcGradient_d.cu) | 無変更 (閉性が前提) |
| 粘性 over-relaxed | [viscousFlux_d.cu](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) | 無変更 |
| block-DPLUR 陰解法 | [timeIntegration_d.cu](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) | 無変更 |
| 構造量 (fx/dcc) | [calcStructualVariables_d.cu](../../solver_density_cuda/cuda_forge/calcStructualVariables_d.cu) | 無変更 (CV 重心と面重心から再計算) |

## 3.5. M1 実装方針 (確定): 双対を primary mesh として書き出す

設計検討の結果、**M1 は solver 側を一切変更せずに実現できる**ことが判明した (当初計画の
`readMesh` dual 分岐や新 BC カーネルは不要だった)。具体的には:

- `convertGmshToForge` が node モードで `buildMedianDual()` → `replacePrimalWithDual()` を実行し、
  双対メッシュ (CV=ノード、内部面=双対面、境界面=半割双対面) を **そのまま `/MESH /PLANES /CELLS
  /BCONDS /VALUE` に primary mesh として書き出す** ([gmshReader.hpp](../../solver_density_cuda/mesh/gmshReader.hpp))。
- solver の `readMesh` はそれを通常メッシュとして読み、**境界半割面 (1 セルの境界 plane) から既存ロジックで
  自動的にゴースト CV を生成する**。よって `setMeshMap_d`・構造量カーネル・対流フラックス・既存 BC
  カーネル (`wall_d`/`inlet_*`/`outlet_*`/`slip_d`) はすべて無変更で node-centered を扱える。
- **境界はゴースト経由の弱形式**: 各境界ノードの半割面に対しゴースト CV を置き、既存 BC カーネルが
  ゴースト状態を設定 → 対流フラックスがゴースト経由で境界フラックスを与える (cell-centered と同型)。
  ゴースト位置は `pc = node + h·n_out` (h=0.5√双対体積) で非退化にする (M1 は非粘性 1 次のため位置は
  フラックスに影響しない)。**より厳密な弱形式境界カーネルは M2+ の精度改善として残す**。
- **可視化 (対応済み)**: node モードでは `replacePrimalWithDual()` が primal セル接続を `/VIZMESH/CONNE`
  に退避し、`writeInputH5`・`output.cpp` がそれを使って **primal セルトポロジ + `Center='Node'`** で
  XDMF を書く (CV index == primal node index)。cell モードは従来どおり `Center='Cell'`。
  両モードとも自己整合な XDMF/HDF5 を出力し ParaView で読める。

## 4. 書き換えが必要な箇所 (M2 以降)

- **境界条件** ([boundaryCond_d.cu](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu)):
  現状はゴーストセル法 (`bplane_cell` / `bplane_cell_ghst` に対称/反射値)。node では境界ノードが
  境界上に乗りゴーストが置けない。新 `boundaryFlux_node_d.cu` で半割双対面に弱形式フラックスを
  `atomicAdd(res_*[node])`。ディスパッチは [boundaryCond.cpp](../../solver_density_cuda/boundaryCond.cpp) でモード分岐。
- **readMesh のゴースト生成** ([mesh.cpp:401-455](../../solver_density_cuda/mesh/mesh.cpp)):
  node モードでスキップし、`nCells_ghst = 0`、境界ノード + 半割面を構築。
- **I/O** ([output.cpp](../../solver_density_cuda/output/output.cpp)): node モードで XDMF
  `Center='Cell'`→`'Node'`、`Dimensions=nCells`→`nNodes`。Topology/Geometry は可視化用 primal を流用。
- **初期条件** ([setInitial.hpp](../../solver_density_cuda/input/setInitial.hpp)): 領域判定ループを node 座標へ。
- **プローブ** ([point_probes.cpp](../../solver_density_cuda/probe/point_probes.cpp)): KD-tree を CV 中心=ノード座標で構築。
- **implicit/AMG の edge 分類** ([mesh.cpp:742-780](../../solver_density_cuda/mesh/mesh.cpp)):
  `nNormalPlanes` ベースの CV 隣接で内部 edge と境界 edge を厳密に分類。
- **軸対称** (副次, [variables.cpp](../../solver_density_cuda/variables.cpp) /
  [axisymmetricSource_d.cu](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)):
  双対体積に `r_node` 重み、`r_eff = volume/A_planar` を双対で再定義、軸上半割マスクをノード版に移植。

## 7. M3: node-centered 弱形式境界 (実装方針)

理論は [theory.md](theory.md) §6。**cell-centered のゴースト経路 ([boundaryCond_d.cu](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu))
は一切変更せず**、node 専用の弱形式を**別ファイル** `cuda_forge/boundaryNode_d.cu` (+ `.cuh`) に実装する。

**設計判断 (既存ロジックの流用)**: 各 BC 種別の境界状態 $Q_b$ は、既存 BC カーネルが**従来どおり
ゴースト CV に書き込む値**を流用する (slip 鏡像・inlet 規定・outlet 外挿の物理を再実装しない)。
弱形式カーネルはこの $Q_b$ (= ゴースト値) と半割面幾何を使い、勾配・フラックスへの寄与だけを置き換える。
ゴースト CV は残すが、**勾配・対流フラックスの境界寄与をゴースト平均から弱形式へ差し替える**のが要点。

**node モードのパイプライン差し替え** (`cfg.discretization=="node"` で分岐、cell は不変):
1. **勾配** ([calcGradient_d.cu](../../solver_density_cuda/cuda_forge/calcGradient_d.cu)):
   `calcGradient_cellgather_d` は CSR 中の**内部面のみ** (`ip < nNormalPlanes`) を集約 (node モードで境界面をスキップ)。
   その後 `boundaryGradientNode_d` が各境界半割面で $\phi_b \mathbf{S}_b$ (φ_b=ゴースト値) を raw 勾配へ atomicAdd。
   最後に `calcGradient_2_d` が体積正規化 (閉曲面が内部+境界で閉じる → §6.2)。
2. **対流フラックス** ([convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) wrapper):
   node モードは内部面 `[0,nNormalPlanes)` のみ処理し、境界は `boundaryFluxNode_d` が
   1 次 Riemann/SLAU フラックス $\hat F(Q_i,Q_b,n_b)S_b$ を残差へ atomicAdd。
3. **ロバスト化 (任意)**: 境界ノードに隣接する内部面の MUSCL 再構成を 1 次に落とすオプション
   (高マッハ壁近傍の過大再構成対策)。まず弱形式勾配で改善するか確認し、不足なら導入。

**検証**: bump hiM (Mach1.65) の node 2 次 MUSCL が**発散しなくなる**ことを `check_convergence.py` で確認
(M2 で step~400 発散 → M3 で安定収束が目標)。cell モードの全既存ケースが**従来とビット一致**(無変更) を回帰確認。

## 5. 設定

[solverConfig.hpp](../../solver_density_cuda/input/solverConfig.hpp) に `discretization` (`cell`/`node`, 既定 `cell`)。
設定の意味は [.github/forge-solver-settings.md](../../.github/forge-solver-settings.md) を参照。

## 6. 検証

[.github/plans/discretization-median-dual.md](../../.github/plans/discretization-median-dual.md) §6 を参照。
要点: cell モードで既存ケースが従来と完全一致すること (デグレ無し) を回帰確認した上で、
node モードを同一メッシュで実行し sod の shock 位置、tet ケースの DOF/コスト/精度を比較する。

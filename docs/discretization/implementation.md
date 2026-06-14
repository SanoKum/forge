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
- **軸対称の軸 (R=0) 注意点 (対応済み)**: solver は回転体積を `volume = A_planar × ccy`(ccy=セル重心 R)で
  実行時計算する。node-centered で CV 重心にノード座標 (軸上で R=0) を使うと**軸ノードの回転体積が ≈0 → 時間積分が
  除算爆発**する (case/29 conical で 272 軸 CV が体積比 1.5e20、step3 で NaN)。対策: `replacePrimalWithDual` が
  dual セル重心に**双対 CV の面積加重重心** (`dualCentroid`, 軸上ノードでも R>0) を使う。これで回転体積が正しくなり
  (体積比 2.4e4, cell の 1.2e4 と同等)、平面ケースは無回帰 (内部ノードでは重心≈ノード)。FV 的にも CV 重心が正しい中心。
  **ただし残課題**: これは軸上のゼロ体積爆発を解消するだけで、node-centered の near-axis 精度問題は別途残る。
  case/29 conical Euler の収束場で、node は軸-inlet 角に非物理 P オーバーシュート (P>chamber Pt)・軸-exit に Mach
  スパイクが出る (bulk は cell と一致)。軸 BC・inlet 角・r 重みの near-axis 処理は今後の課題 ([architecture-axisym-axis-singularity] と関連)。
  **試行 (いずれも発散→不採用)**: (1) 軸ノードで `roUy=0` 毎ステージ強制 (roe から半径 KE 除去) → 過渡で roe<0→全域 NaN。
  (2) 軸ノードで半径方向圧力ソース抑制 → ソースは軸 CV の釣り合いに load-bearing で step~51 破綻。
  **結論: 半径ソースは軸 CV にも必要で安易な対称強制は不可**。baseline (ソース維持) は収束し corner オーバーシュートのみ残る。
  infra (`mesh.axis_flag_d`, `enforceAxisSymmetry_d`) は残置・既定 off。corner の正攻法は別途要検討。
  **SU2 流の確認 (ユーザ提案)**: SU2 はノードベース (中点双対) で、軸を `MARKER_SYM=(axis)` の**対称面**として扱い、
  対称ノードで**法線 (=半径) 運動量成分を残差から射影**する。forge で同方式を試したが、forge の implicit は
  **block-DPLUR (5×5 連成)** のため:
  - (3) `res_roUy=0` 射影 (全 flux+source 後) → 単独では連成 solve が補正を漏らし**無効** (Uy 164→136 でほぼ不変)。
  - (4) implicit commit で `roUy=0` 強制 → 線形 solve と不整合で**発散級に悪化 (Mach~1000)**。
  → **SU2 流対称面を効かせるには block-DPLUR の Jacobian row を軸ノードで整合修正する必要**があり (SU2 は
  Jacobian を修正している)、commit/残差だけの hack では不可。これは delicate な implicit-solver 改修なので
  open issue とし、infra (`zeroAxisRadialResidual_d`、commit の `axis_flag` 引数) は残置・既定無効。
  baseline (無介入) は収束し corner 1 セルのオーバーシュートのみ。
  **explicit 検証 (ユーザ要請)**: conical Euler を explicit で回すと **node は step1 で即発散** (cfl=0.1 でも)。
  原因は **near-axis の explicit 剛性**: 軸 CV は r 重みで面積→0 のため spectral radius→0、explicit 更新
  `res·CFL/spectral` が発散 (cell は step3483 まで持つが非粘性ノズル startup でやはり発散)。よって本ノズルは
  **explicit では node が回らず implicit (block-DPLUR 無条件安定) が前提**。explicit 検証には near-axis の
  dt/spectral フロアが別途要る。射影 (roUy 残差=0) は今 node+axisym で既定 ON だが、explicit が回らない・implicit は
  連成で漏れるため、単独では corner を直せない。
  **入口境界の確認 (ユーザ要請)**: 収束場の入口面 (x=inlet) を r 方向に見ると **r≳0.003 で入口 BC は正常**
  (P≈4MPa=chamber Pt・軸方向流入正)。異常は **軸-inlet 角 (r<0.002) のみ**: 軸方向**逆流 Ux≈-798**・Uy≈+164・P 6.82MPa。
  → **入口 BC 自体は問題なし。軸特異性が角を多方程式的に汚染** (radial だけでなく axial 逆流・P も)。roUy 単独射影では
  直らず、near-axis は全方程式整合 (対称面 + dt/spectral フロア + 連成 Jacobian) の包括対策が要る。
  **SU2 処方の全試行 (ユーザ提案: y<EPS で軸ソース OFF + roUy/Uy=0)**: 両処理を実装 (`axisymmetricSource` の
  `axis_flag` でソース OFF、`enforceAxisSymmetry_d` で roUy=0/Uy=0/roe 整合) し、**併用かつ全配置** (applyBconds
  =更新前、block-DPLUR commit 内、commit 後 post-update) で試したが、**implicit は全て Mach~1000 に発散**。
  理由: **block-DPLUR (defect-correction) は solve の外で roUy/roe を手術されると非整合**になり増幅する。SU2 が安定なのは
  **対称条件を Jacobian (行) の中で課す**から。explicit は recipe 併用でも軸 CV が **step1 で発散** (別系統の near-axis
  剛性)。**結論: SU2 処方は原理的に正しいが、forge では「外部状態手術」では不可。block-DPLUR の軸ノード 5×5 Jacobian で
  roUy 行を decouple する実装 (= SU2 と同じ内挿) が必須**。全フック (axis_flag/enforceAxisSymmetry/zeroAxisRadialResidual)
  は実装済み・既定無効で残置。baseline (無介入・ソース ON) は収束し corner 1 セルのみオーバーシュート。
  **in-Jacobian decouple も試行 (SU2 と同じ「Jacobian 内で対称化」)**: block-DPLUR の軸ノード 5×5 で roUy 行 (index2) を
  単位行に置換し rhs[2]=0 (`timeIntegration_d.cu` の `axis_flag`) → solve の中で一貫して dq_roUy=0。結果:
  **bulk は正しい (mean Mach 1.883) が、corner 2 セル (inlet-axis, exit-axis) が Mach~1007 に発散**(P 22MPa)。
  **重要: corner の発散は roUy ではなく軸方向運動量・密度**(Ux 逆流 −798 と同系)。**= corner CV (r=0 で両境界半割面が
  r 重み→0 の極小 CV) は多方程式的に ill-posed で、roUy の対称化だけ (外部手術でも in-Jacobian でも) では直らない**。
  → corner は構造的処置が要る可能性 (例: 軸-境界角ノードを隣接へ merge/除外する、専用の角クロージャ)。
  baseline (無介入) が現状ベスト (収束・bulk 正しい・corner は P≤6.8MPa の局所オーバーシュートに留まる)。全フックは既定無効。
  **SU2 メッシュ/角処理の調査結果 (ユーザ要請, `run_su2cmp_su2_euler`)**:
  - SU2 メッシュ = forge node メッシュと**同一** (NPOIN=27472、軸-inlet 角は**同じ共有ノード id 2**, x=-0.0587/y=0、
    marker は inlet101/axis272/wall272/outlet101 で forge bcond と一致)。**トポロジは差が無く、差は BC 処理のみ**。
  - **SU2 の共有角は clean**: 角で Ux=+54.78 (正常流入, 逆流でない)・**Uy=0 厳密**・P=3.99MPa (=chamber, overshoot 無し)。
    SU2 は**軸全体で Uy=0 を厳密強制** (max|Uy|=0)。
  - forge の roUy 対称化 (外部手術 / in-Jacobian 単行 decouple) は**全て発散**。SU2 が clean なのは、SU2 の対称面が
    **残差 + Jacobian + 解を一体で扱う coherent な対称面演算子 (CSymmetryPlane 相当)** だから。forge の piecemeal な
    1 行 decouple は、軸 CV の roUy 連成が持つ陰的減衰を壊し、軸方向・密度を不安定化させる (Mach~1000)。
  - **結論**: corner を SU2 並に clean にするには、**SU2 の対称面演算子を一体移植**する必要がある:
    (i) 対称面フラックスを圧力のみ (質量/エネルギー流束 0)、(ii) 残差から法線運動量成分を除去、(iii) 同じ射影を
    block-DPLUR の 5×5 Jacobian へ整合適用 (1 行 hack でなく射影行列 P=I-nn^T で両側変換)、(iv) 解にも v_n=0。
    SU2 ソースは本リポジトリにバイナリのみで未読。bulk は完成済みのため残作業は corner=この対称面演算子の移植。
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

### 7.0 [解決] node-centered 軸対称の near-axis corner

長く open だった「軸-境界 (inlet/outlet) 角 CV の P オーバーシュート/逆流」を **2 修正の併用で解決** (cell/SU2 一致):
1. **軸 (R≤eps) で軸対称ソース OFF** (`axisymmetricSource_d` の `axis_flag`、node+axisym のみ)。ソース
   `res_roUy += P·A_planar` は**唯一 r 重みされない項**で `source/volume = P/r → ∞` (r→0) が発散源。対流項は
   r_face/r_cell がキャンセルし有限。軸でソースを切るのが本質 (SU2 の y<EPS 相当、ユーザ提案)。
2. **境界半割面の重心に真の面積加重重心 (R≥0) を使う** (`dualBnodeCent`)。旧 `node+h·n_out` 便宜は pcy≈0/<0 で
   入口/出口 BC の r 重み実効面積を ~0 にし **BC が軸近傍 corner CV に届いていなかった**。真の重心 (R>0) で入口が
   コーナーに届き P=chamber に。(`axisCentroidShift=1` で CV セル中心に面積加重重心を使うのも維持: r_cell と r_face を
   一致させ対流項の r をキャンセルさせるため。)

**検証 (case/29 conical & bell, 軸対称 Euler implicit)**: PASS。corner P=3.99MPa(=chamber, overshoot 0 セル、旧 6.82e6)、
Ux=+51.7(流入, 旧 -798 逆流)、Uy≈0 — **SU2 (P=3.99,Ux=+54.8,Uy=0) と一致**。near-axis<3% mean|ΔM|=0.016/max0.049
(旧 max0.95) と中/外帯並み。平面 (bump loM) 無回帰。**explicit は startup の exit-wall 角で依然脆く発散** (cell explicit も
同様; 軸でなく supersonic-startup ロバスト性) → implicit が実用。

### 7.1 診断結果 (重要): hiM 発散は「境界値」でなく「壁近傍 2 次再構成」の問題
弱形式実装の前に切り分け診断を実施 (`bndFirstOrder`: 境界隣接 CV の再構成勾配を 0 にし境界側を 1 次化、
[calcGradient_d.cu](../../solver_density_cuda/cuda_forge/calcGradient_d.cu) の `zeroBndNodeGradient_d`、
フラグは `mesh.bnode_flag_d`)。結果:
- **`bndFirstOrder=1` で bump hiM node の NaN 発散 (step~400) が解消** (explicit は ~1e-3 のリミットサイクル、
  **implicit + bndFirstOrder で完全収束 PASS**)。cell モードは flag 既定 0 で無影響 (ビット不変)。
- **→ 真因は近壁 2 次 MUSCL 再構成のロバスト性であって、境界値 (φ_b) ではない**。slip 壁では弱形式の φ_b は
  現状の `0.5(node+ghost)` と一致するため、**弱形式境界だけでは hiM は直らなかったはず**。`bndFirstOrder` が
  実効的な最小修正。
- **弱形式境界 (§7) は引き続き粘性 (壁せん断・熱流の正しい評価) で必要**だが、hiM 安定化とは別件として整理する。
  よって `boundaryNode_d.cu` の新規作成は粘性着手時に回す (既存 BC カーネルの Q_b 流用方針は不変)。

### 7.2 node モード: ゴースト全撤廃へ (段階導入)

理論は [theory.md](theory.md) §6.3/§6.4。node-centered は最終的に**境界ゴーストを全撤廃**し、境界ノード/半割面に
BC を直接課す。段階導入し、まず粘性壁の発散 (case/29 viscous node の exit-lip NaN) を Phase 1 で止める。
専用計画: [`discretization-node-boundary-ghostless.md`](../../.github/plans/discretization-node-boundary-ghostless.md)。

**Phase 1 — 壁ゴースト撤廃 ＋ 壁 Dirichlet ＋ 壁優先コーナー所有** (流入出/slip/axis はゴースト維持):

1. **壁優先所有・壁 plane 不出力** ([gmshReader.hpp](../../solver_density_cuda/mesh/gmshReader.hpp)
   `buildMedianDual`/`replacePrimalWithDual`): 各境界ノードを優先度 (wall>inlet>outlet>slip/axis) で 1 bcond に
   所有。**wall 所有ノードは境界 plane を出力しない** (→ ghost 生成されない) が、wall bcond のノード列 (iCells) は
   h5 に残す。コーナー (wall 所有) は出口 plane を持たない → 出口ゴーストの矛盾消滅。
2. **wall_flag_d** ([mesh.cpp](../../solver_density_cuda/mesh/mesh.cpp) `setMeshMap_d`,
   [mesh.hpp](../../solver_density_cuda/mesh/mesh.hpp)): `axis_flag_d` と同パターンで wall 種別 bcond の iCells を 1 に。
3. **壁 no-slip カーネル** ([axisymmetricSource_d.cu](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)
   に追記、軸カーネルと同型):
   - `enforceWallNoSlip_d` (一度きり state 初期化): 壁ノードで `roe -= 0.5*(roUx²+roUy²+roUz²)/ro` 後
     `roU{x,y,z}=0`, `U{x,y,z}=0` (`enforceAxisSymmetry_d` の 3 成分版)。IC 確定後 1 回。
   - `zeroWallMomentumResidual_d` (毎反復残差射影): 壁ノードで `res_roU{x,y,z}=0` (`zeroAxisRadialResidual_d` の
     3 成分版)。**残差射影なので block-DPLUR と整合** (state 直書きは Mach1000 発散実績)。
4. **ディスパッチ/呼び出し**: [boundaryCond.cpp](../../solver_density_cuda/boundaryCond.cpp) は node モードで
   wall/wall_isothermal を no-op 化。[main.cpp](../../solver_density_cuda/main.cpp) で `enforceWallNoSlip` を IC 後 1 回、
   `zeroWallMomentumResidual` を `assembleResidual` 末尾 (軸射影 `zeroAxisRadialResidual` の後) に毎反復。
   全変更を `discretization=="node"` ＋ `wall_flag_d` でゲートし cell モードは無変更。

**なぜ壁ゴーストを消してよいか** (theory §6.3): 断熱壁の壁半割面フラックスは恒等的に 0 (u·n=0、運動量は
Dirichlet、熱流束 0)、近壁せん断は壁ノード (u=0) と内部ノードを繋ぐ**内部双対面**が担うため。等温壁は壁半割面に
熱流束項を別途加える (Phase 1 拡張、case/29 は断熱で不要)。

**Phase 2 — 残り (inlet/outlet/slip/axis) のゴースト撤廃**: 半割面ジオメトリ (`dualBnodeId/Vect/Cent`,
`dualBcondOffset`) を device 転送し、§7 の `boundaryFluxNode_d`/`boundaryGradientNode_d` で半割面上に直接 flux を
評価、共有 flux ループは node モードで境界 plane をスキップ。implicit 寄与も flux/Jacobian に内在化。

## 5. 設定

[solverConfig.hpp](../../solver_density_cuda/input/solverConfig.hpp) に `discretization` (`cell`/`node`, 既定 `cell`)。
設定の意味は [.github/forge-solver-settings.md](../../.github/forge-solver-settings.md) を参照。

## 6. 検証

[.github/plans/discretization-median-dual.md](../../.github/plans/discretization-median-dual.md) §6 を参照。
要点: cell モードで既存ケースが従来と完全一致すること (デグレ無し) を回帰確認した上で、
node モードを同一メッシュで実行し sod の shock 位置、tet ケースの DOF/コスト/精度を比較する。

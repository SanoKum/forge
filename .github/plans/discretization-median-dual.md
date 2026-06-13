# cell-centered / node-centered (median-dual) 両対応化

## メタ

- **area**: `architecture`
- **status**: `in_progress`
- **related_docs**:
  - `docs/discretization/theory.md`
  - `docs/discretization/implementation.md`
  - `docs/architecture/overview.md`
- **related_plans**:
  - `architecture-axisym-axis-singularity.md` (軸対称 near-axis、本対応で改善余地あり=副次)
- **created**: `2026-06-14`
- **owner**: `CFD Dev`

## 1. 目的

現状 cell-centered 固定の forge を、`solverConfig` のフラグで node-centered (中点双対 median-dual)
にも切替できる両対応にする。完了時には、同じ実行ファイル・同じメッシュで cell / node 両方の
計算が走り、tet 非構造ケースで node の自由度削減 (≈1/5.5) によるコスト改善と精度を実測できる。

## 2. スコープ

- **やる**:
  - median-dual メッシュ生成 (前処理)。双対体積・エッジ双対面・CSR・境界半割面・閉性チェック。
  - `discretization: cell | node` フラグと `readMesh` の dual 分岐。
  - node-centered 境界フラックス (弱形式) カーネル。
  - I/O (`Center='Node'`)・初期条件・プローブの node 対応。
  - cell モードの完全な後方互換 (既存ケースが従来と一致)。
- **やらない (将来 / 別 plan)**:
  - atomicAdd を避ける gather 版への全面最適化 (M2 以降の性能改善として切り出し)。
  - 軸対称双対体積の完全対応は副次 (M3)。`architecture-axisym-axis-singularity.md` と連携。
  - node-centered 専用の高次スキーム研究。

## 3. 関連 docs と前提

- 理論: `docs/discretization/theory.md` (median-dual の CV/面/閉性/境界半割面)。
- 実装方針: `docs/discretization/implementation.md` (discretization mesh 抽象、カーネル対応表)。
- 前提: GPU カーネルは「CV + 面 + 接続」抽象だけで書かれており、内部スキームは双対上で無変更で成立
  (調査済み)。差異は前処理と境界・I/O に局在する。

## 4. 設計方針

詳細は `docs/discretization/implementation.md`。要点のみ:

- **供給側 2 実装**: cell モード = 現行 gmshReader (`/CELLS` を CV)、node モード = 新 `dualMeshBuilder`
  (`/DUAL` を CV)。solver は汎用 `mesh` 抽象 (`nCells/nPlanes/map_plane_cells/CSR`) を消費。
- **`nCells` をリネームしない**。node モードでは `nCells == nNodes`, `nCells_ghst == 0` と読み替える。
- **境界のみカーネル新規**: 半割双対面 + 弱形式フラックス。ゴースト経路は node モードで未使用化。
- **閉性 assert を前処理に必須実装** (`Σ双対面ベクトル = 0`、`Σ双対体積 = Σ領域体積`、`Σ半割面積 = 境界面積`)。

## 5. 実装ステップ

### M1 — 双対生成 + Euler 1 次が回る (最小マイルストーン)
1. **[完了]** `buildMedianDual()` を `mesh/gmshReader.hpp` に追加 (2D: triangle/quad)。双対体積・
   エッジ双対面 (= primal plane と 1:1)・境界半割面・閉性チェック・`/DUAL` HDF5 出力。
   `convertGmshToForge` が `discretization=="node"` で呼ぶ。`input/solverConfig.{hpp,cpp}` に
   `discretization` フラグ追加 (既定 `cell`)。**3D median-dual は未実装 (M4)**。
2. `mesh/mesh.cpp` `readMesh` に dual 分岐 (`nCells:=nNodes`, `nCells_ghst:=0`)、ゴースト生成スキップ。
3. 境界最小実装: slip + supersonic in/out の node 弱形式カーネル (`cuda_forge/boundaryFlux_node_d.cu`)。
   ディスパッチ `boundaryCond.cpp`。
4. 対流 1 次・粘性 OFF・RK 陽解法で **2D ケース**を実行。
   **注意**: `case/05.sod_shock_tube` の `1D_shock_tube.msh` は実際には 3D hex 押し出し
   (nBPlanes=1998=499×4+2) なので 2D builder の対象外。M1 の最初の 2D 検証は
   `case/24.laminar_channel_bl` (quad)・`case/08.bump`・`case/23.axi_nozzle` の `*_2d.msh` を使う。
   3D の sod は 3D median-dual (M4) 完成後。

### M2 — 2 次 MUSCL + 勾配 + implicit + tet 比較
5. Green–Gauss 勾配を dual 有効化、MUSCL 2 次。`calcGradient_d.cu` (無変更で動く想定、検証対象)。
6. `case/08.bump` で wall BC node 版 + block-DPLUR dual。
7. tet 非構造ケースで cell vs node を同一メッシュ比較 (DOF 数・実行時間・精度)。← 動機の本丸。

### M3 — 粘性 + 軸対称 (副次)
8. viscous dual 有効化、層流ケース (`case/15` or `24`)。
9. 双対体積の `r` 重み + 軸上マスク、`case/23.axi_nozzle` で near-axis 改善有無を比較。

### M4 — I/O・プローブ・3D 混合要素
10. `output.cpp` の `Center='Node'`、`setInitial.hpp` / `point_probes.cpp` の node 化、3D tet/prism ケース。

## 6. 検証

- **単体 / ビルド**: native (WSL) で `solver_density_cuda` ビルド成功。
  `convertGmshToForge` を tet メッシュに適用し `/DUAL` 生成・閉性 assert・体積総和一致・負体積無しを確認。
- **検証ケース**: `case/05.sod_shock_tube` (Euler/shock) → `case/08.bump` (内部流/壁/implicit) + tet 比較
  → `case/15` or `24` (粘性) → `case/23.axi_nozzle` (軸対称) → 3D 混合要素。
  各 run は新規 `run_*` に出力、`residual_history.png` 生成、NaN/発散・全残差列の収束を AGENTS.md 手順で確認。
- **判定基準**:
  - cell モード: 既存ケースの残差・物理量が従来と**完全一致** (デグレ無し)。
  - node モード sod: shock 位置・一定流保存が cell 版/解析解と整合。
  - tet ケース: 同一メッシュで DOF 数減・コスト改善・許容精度を定量比較。

## 7. 影響範囲

- 触るモジュール: `mesh/gmshReader.hpp`, `mesh/mesh.cpp`, `input/solverConfig.hpp`,
  `cuda_forge/boundaryFlux_node_d.cu` (新規), `boundaryCond.cpp`, `output/output.cpp`,
  `input/setInitial.hpp`, `probe/point_probes.cpp`, `variables.cpp` (軸対称・副次)。
- 既存ケース: cell モードは無影響 (フラグ既定 `cell`)。
- ドキュメント: `docs/discretization/{theory,implementation}.md`, `docs/index.md`,
  `docs/architecture/overview.md`, `.github/forge-solver-settings.md`。

## 8. 完了条件

- [x] 関連 `docs/discretization/theory.md` 作成済み
- [x] 関連 `docs/discretization/implementation.md` 作成済み
- [ ] M1 実装・検証完了 (§6 の sod 判定を満たす)
- [ ] M2–M4 実装・検証完了
- [ ] `.github/plans/README.md` の状態を更新
- [ ] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-06-14` — 初稿。調査 (現行は CV+面抽象で内部カーネルは双対上無変更で成立) を反映し、
  両対応方針と M1–M4 ロードマップを確定。docs theory/implementation を整備。
- `2026-06-14` — **M1-1 完了 (2D)**。`buildMedianDual()` (2D triangle/quad) + `discretization`
  フラグ + `/DUAL` HDF5 出力を実装し、`convertGmshToForge` をネイティブビルド
  (`build-dual`, RTX 3060 arch86)。`case/24.laminar_channel_bl/mesh/channel.msh` で検証:
  双対体積総和 relErr=5.7e-8、閉性 max|ΣdS| 正規化=1.5e-5 (float32 丸め)、境界半割面積=primal
  面積 (0.01/0.01/0.2 一致)、nNodes(CV)=13065・nDualFaces=25864・nBHalf=532。アルゴリズム正当性を確認。
  **副次発見**: `case/05.sod_shock_tube` は 3D hex 押し出しのため 2D 検証は `case/24` 等を使用 (§5 M1-4)。
  **次**: M1-2 (`readMesh` の dual 消費) → M1-3 (node 境界カーネル) → M1-4 (2D ケース実走)。

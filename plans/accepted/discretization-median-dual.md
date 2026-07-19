# cell-centered / node-centered (median-dual) 両対応化

## メタ

- **area**: `architecture`
- **status**: `done`  <!-- 親計画としてクローズ (2026-07-19)。現役の入口は子 plan (末尾の変更ログ参照) -->
- **related_docs**:
  - `methods/discretization.md`
  - `methods/discretization.md`
  - `methods/architecture/overview.md`
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

- 理論: `methods/discretization.md` (median-dual の CV/面/閉性/境界半割面)。
- 実装方針: `methods/discretization.md` (discretization mesh 抽象、カーネル対応表)。
- 前提: GPU カーネルは「CV + 面 + 接続」抽象だけで書かれており、内部スキームは双対上で無変更で成立
  (調査済み)。差異は前処理と境界・I/O に局在する。

## 4. 設計方針

詳細は `methods/discretization.md`。要点のみ:

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
2. **[完了・方針変更]** 当初は `readMesh` の dual 分岐を予定したが、検討の結果 **solver 無変更**で
   実現できると判明。`replacePrimalWithDual()` ([gmshReader.hpp](../../solver_density_cuda/mesh/gmshReader.hpp))
   が双対を primary mesh として `/MESH /PLANES /CELLS /BCONDS /VALUE` に書き出し、solver の `readMesh` が
   境界半割面から自動でゴースト CV を生成する。`setMeshMap_d`・カーネル・BC はすべて無変更。
3. **[完了・方針変更]** 新 BC カーネルは M1 では**不要**。境界半割面にゴースト CV を置き、既存
   `wall_d`/`inlet_*`/`outlet_*`/`slip_d` がゴースト状態を設定 → 対流フラックスが弱形式の境界
   フラックスを与える (cell-centered と同型)。厳密な弱形式 BC カーネルは M2+ の精度改善に回す。
4. **[完了]** 対流 1 次・粘性 OFF・RK3 陽解法で 2D ケースを cell/node 両モードで実走・比較。
   検証ケースは `case/24.laminar_channel_bl` (quad)。
   **注意**: `case/05.sod_shock_tube` の `1D_shock_tube.msh` は実際には 3D hex 押し出し
   (nBPlanes=1998=499×4+2) なので 2D builder の対象外。3D の sod は 3D median-dual (M4) 完成後。

### M2 — 2 次 MUSCL + 勾配 + implicit + tet 比較 **[一部完了]**
5. **[完了]** Green–Gauss 勾配 + MUSCL 2 次を dual で検証 (カーネル無変更で動作)。
6. **[完了]** `case/08.bump` (2D 版 `bump_2d.geo` を新規生成) で loM/hiM・explicit/implicit を cell vs node 検証。
7. **[完了]** 三角形 (tet 相当) メッシュ `bump_2d_tri.geo` で DOF/コスト比較。
   **結果は §9 の `2026-06-14 M2` 参照**。要点: 低マッハ 2 次は node≈cell・implicit も良好だが、
   **高マッハ 2 次 MUSCL が壁近傍で発散** (M1 のゴースト便宜由来) → **M3 弱形式境界が必要**と確定。

### M3 — 弱形式境界 (前倒し) + 粘性 + 軸対称
7.5. **[診断完了]** hiM 発散の切り分け診断 `bndFirstOrder` を実装 (境界隣接 CV の再構成勾配を 0 化)。
   **結果: hiM node の NaN 発散は「壁近傍 2 次再構成のロバスト性」が真因で「境界値」ではないと確定**。
   `bndFirstOrder=1` で発散解消 (implicit と併用で完全収束 PASS)。詳細は §9 `M3 診断`。
   → 弱形式境界 (`boundaryNode_d.cu`) は hiM の特効薬ではなく**粘性 (壁せん断/熱流) 用**として粘性着手時に作る。
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
- ドキュメント: `methods/discretization/{theory,implementation}.md`, `methods/index.md`,
  `methods/architecture/overview.md`, `procedures/solver-settings.md`。

## 8. 完了条件

- [x] 関連 `methods/discretization.md` 作成済み
- [x] 関連 `methods/discretization.md` 作成済み
- [ ] M1 実装・検証完了 (§6 の sod 判定を満たす)
- [ ] M2–M4 実装・検証完了
- [ ] `plans/README.md` の状態を更新
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
- `2026-06-14` — **M1-2/M1-3/M1-4 完了 (2D)・方針変更**。検討の結果、**solver は一切変更せず**に node-centered を
  実現 (当初の `readMesh` dual 分岐・新 BC カーネルは不要)。`replacePrimalWithDual()` で双対を primary mesh
  として書き出し、solver の `readMesh` が境界半割面からゴースト CV を自動生成 → 既存カーネル・BC が無変更で動く。
  境界はゴースト経由の弱形式 (既存 `wall_d`/`inlet_*`/`outlet_*` を再利用)。
  **検証** (`case/24.laminar_channel_bl`, 非粘性 1 次 RK3 陽解法 3000 step, cell vs node 同一メッシュ):
  両モードとも exit0・NaN/Inf=0・3000 step 完走。残差はほぼ一致 (rms_ro 2.67e-6 vs 2.75e-6、
  rms_roUx 9.66e-4 vs 9.59e-4、rms_roe 0.82 vs 0.84、残差履歴の形状もほぼ重なる)。
  最終場も一致 (ro_mean 1.1788 vs 1.1785、Ux_mean 13.60 vs 13.53、P_mean 101537 vs 101506、
  物理的に健全 Ps≤P≤Pt・ro>0・NaN無し)。run: `case/24.laminar_channel_bl/run_dual_{cell,node}_m1/`
  (比較図 `m1_cell_vs_node_residuals.png`)。**注**: rms_roe プラトーは低マッハ非粘性 duct の陽解法 3000 step
  の緩慢収束で両モード共通 (dual 起因ではない)。`output.cpp` の XDMF/CONNE は cell topology 前提のため
  node 可視化は未対応 (M4)。
  **次**: M2 (2 次 MUSCL + 勾配 + implicit を dual で有効化、tet 比較) → M3 (粘性・軸対称) → M4 (3D median-dual・I/O)。
- `2026-06-14` — **M2 検証 (2D, 非粘性)**。`case/08.bump` の 2D 版 (`mesh/bump_2d.geo` 新規・quad、
  `mesh/bump_2d_tri.geo` 新規・unstructured 三角形) を生成し cell vs node を比較。出力修正 (`/VIZMESH`) で
  両モード ParaView 可読。
  - **M2a 2 次 MUSCL+勾配**: 低マッハは node≈cell。`case/24` channel・bump loM (subsonic) とも完走・NaN無し、
    最終場一致 (bump loM: ro_mean 1.2131 vs 1.2126、Ux_mean 175.7 vs 176.0、P range 同等)。**勾配/MUSCL/limiter
    カーネルは dual 上で無変更で動作**。
  - **高マッハの限界 (重要)**: bump hiM (Mach1.65 supersonic) は **node の 2 次 MUSCL が step~400 で発散**
    (cell は安定)。CFL 0.25・limiter venkata/barth いずれでも発散 (~280-410) → CFL/limiter 強度の問題ではない。
    **1 次は node でも安定** (rms_ro 5.8e-6 収束)。発散場所を pre-NaN snapshot で特定 = **下壁 (slip) の bump
    後縁 x≈1.98,y≈0.009** (強膨張域)。根本原因は **M1 のゴースト便宜による境界ノード勾配の不正**で、強勾配の
    高マッハで 2 次再構成が増幅して破綻。低マッハは勾配が緩く顕在化しない。**→ M3 の弱形式境界で解消すべき**
    (境界 1 次化の暫定策も可)。
  - **M2b implicit (block-DPLUR)**: dual で**無変更動作**。bump loM implicit cfl_pseudo=5 で両モード完走・NaN無し、
    node は rms_ro 8.6e-7 まで収束 (cell 7.3e-4 より深い)。
  - **M2c DOF/コスト (動機本丸)**: 三角形 bump で **node は DOF 半減 (CV 6436 vs cell 12548、~1.95×減)**。
    ただし explicit の per-step 実時間はほぼ同等 (node 18.3s vs cell 19.3s/4000step)。理由: **対流フラックス
    コストは面 (エッジ) 数律速で、双対の内部面数 (~19k) は primal の面数とほぼ同じ**。DOF 半減の効果は
    **メモリと implicit 線形系サイズ**に出る (explicit per-step では出にくい)。
  run: `case/08.bump/run_dual_{loM,hiM}_{cell,node}_m2/`, `run_dual_imp_loM_{cell,node}_m2b/`,
  `run_dual_tri_loM_{cell,node}_m2c/` (比較図 `m2_bump_cell_vs_node.png`)。
  **次 (M3 を前倒し)**: 弱形式境界カーネル (境界ノード半割面に直接フラックス) で高マッハ 2 次の壁近傍を解消 → 粘性。
- `2026-06-14` — **SU2 参照検証 + 入口/軸重複の切り分け (自律検証)**。
  - **入口×軸 重複の確認 (ユーザ仮説)**: 角 CV (id=2) は inlet_Pressure と slip(軸) の両 bcond に属するが、
    `inlet_Pressure_d`/`slip_d` とも **ghost のみ書込み・内部ノード状態は上書きしない** (状態競合は無し)。
    真因は、角 (r≈0) で**両境界半割面が r 重みで面積→0** になり BC が実質効かず、軸対称の半径方向圧力ソースが
    不釣合になること。入口面を r 方向に見ると **r≳0.003 で入口 BC は正常** (P≈4MPa=Pt・流入正)、異常は軸-inlet 角のみ。
  - **cell-centered は SU2 と一致 (検証済み)**: `run_su2cmp_su2_euler` の SU2 軸対称 Euler を参照に、conical 中心線
    Mach を比較 — cell は exit Mach 3.79 vs SU2 3.71、平均 |ΔM|=0.073、Pmax 一致 (3.99MPa)、overshoot 無し。
  - **node は bulk で SU2/cell に一致・アーティファクトは near-axis に限局**: conical を半径帯で node vs cell 比較すると
    **r 30-50%: 平均|ΔM|=0.014、r 70-95%: 0.011** (max Mach も一致)、**r<3% のみ大 (max|ΔM|=0.95, exit 軸 Mach 5.18,
    inlet 角 P 7.44MPa)**。中心線プロット (`m3_conical_centerline_su2_cell_node.png`) でも cell/node が SU2 と重なり、
    node だけ inlet 角・exit 軸に数点の外れ値。bell も同様 (near-axis のみ)。
  - **結論 (確定)**: **node-centered 軸対称 Euler は bulk で正しい (SU2/cell に ~1% 一致)。欠陥は r=0 軸線の局所処理のみ**。
    explicit は軸 CV の spectral→0 で step1 発散 (implicit 前提)。near-axis 包括対策 (全方程式対称+dt/spectral フロア+
    連携 Jacobian) が唯一の残課題。run: `run_dual_eul_{conical,bell}_{cell,node}_m3/`, `run_su2cmp_su2_euler/`。
  - **explicit 追検証 (detectNaN + subsonic)**: supersonic conical explicit は **step2 で exit-wall 角 (x=exit,r=rmax) の
    outlet 近傍**が発散 (detectNaN: ro NaN @ CV(x=0.086,r=0.033)、P→2.5e11)。さらに**低圧subsonic (Pt=110k) でも explicit は発散**。
    → **node-centered 軸対称は Mach によらず explicit で回らない (軸 CV が explicit 不安定) → implicit 必須**。
    よって SU2 流残差射影 (res_roUy=0) の explicit 有効性は本系では検証不能。検証には平面 (case/24, 既に node 粘性壁 OK) か、
    軸を持たない gentle 系が要る。**結論不変: node 軸対称は implicit 前提、near-axis 包括対策が残課題。**
- `2026-06-14` — **M2 収束判定の訂正 (反省)**。上記 M2 報告の多くは **NaN チェックのみで収束を確認しておらず誤り**だった。
  新ツール `solver_density_cuda/tools/check_convergence.py` で全残差列を再判定した結果:
  - **explicit の channel M1・bump loM/hiM は全て未収束のトランジェント** (0.6〜2.5 桁しか低下せず `falling`)。
    cell と node が**ほぼ同一に推移**するのは discretization 等価性の傍証ではあるが「収束・一致」ではない。
  - **bump loM implicit (cfl_pseudo=5, 8000step)**: **node は PASS (5.1 桁低下)、cell は ~2.3 桁で停滞 (NOT CONVERGED)**。
    cfl_pseudo=20 は両モードとも NaN 発散 (memory の cfl50 安定は別 case)。
  - 教訓を `AGENTS.md` の「収束確認 (必須)」に反映: **check_convergence.py の VERDICT 提示を必須化**。
    今後は「収束/一致」と書く前に必ずツールを通す ([[convergence-check-discipline]])。
  - **収束済み比較を実施** (`case/24` channel implicit, MUSCL, cfl_pseudo=5, 12000step):
    **cell は PASS** (rms_ro/roUx/roe 3.3〜3.5 桁低下・falling、roUy は初期厳密 0 の near-zero)。
    **node は「still converging」** (主残差 3.3 桁 falling、roUy のみ 1.4 桁で収束途中=要追加 step、停滞ではない)。
    well-developed 状態の場は一致: ro_mean 1.17927 vs 1.17928、Ux_mean 35.25 vs 35.20 (0.14%)、P≈101325 一様、NaN無し。
    収束流速 Ux≈35 (M≈0.1、case 設計値) は **未収束 M1 トランジェントの Ux≈13.5 と全く異なる** = M1 比較が
    トランジェント同士だった証左。run: `run_dual_imp_{cell,node}_m2chk/`。
  - ツール改良: init=0 の near-zero 列を inactive 扱い、未収束理由を `still converging`(falling)/`stalled`(flat)/
    `diverged`(NaN/rising) で区別表示。
- `2026-06-14` — **M3 診断: hiM 発散の真因を特定**。弱形式境界の実装前に切り分け診断 `bndFirstOrder`
  (境界隣接 CV `mesh.bnode_flag_d` の再構成勾配を `zeroBndNodeGradient_d` で 0 化=境界側 1 次) を実装し bump hiM で検証:
  - **`bndFirstOrder=1` で node hiM 2 次 MUSCL の NaN 発散 (step~400) が解消**。explicit は ~1e-3 のリミットサイクル
    (min rms_ro 7.8e-4 で振動、発散ではない)、**implicit + bndFirstOrder で全列 PASS 完全収束**。
  - cell モードは flag 既定 0 で無影響 (新コードは `bndFirstOrder==1` で完全ガード)。
  - **結論**: hiM 発散は**近壁 2 次再構成のロバスト性**が真因で、**境界値 (φ_b) ではない** (slip では弱形式 φ_b=現状平均で一致)。
    弱形式境界は hiM の特効薬ではなく**粘性 (壁せん断・熱流) の正確評価用**と整理し、粘性着手時に `boundaryNode_d.cu` を作る。
  - 新フラグ `bndFirstOrder` (mesh: 既定 0) を `solverConfig` に追加。run: `run_dual_hiM_node_m3{diag,imp}/`。
  **次**: 粘性 (viscous) を dual で有効化する際に弱形式境界 (`boundaryNode_d.cu`) を実装。軸対称 (M3 後半) も。
- `2026-06-14` — **node-centered 軸対称 (case/29) + 軸 R=0 修正**。軸ノード(R=0)が CV 中心になると回転体積
  `A_planar×R` が ≈0 → 除算爆発 (case/29 conical で 272 軸 CV、体積比 1.5e20、step3 で NaN)。
  **修正**: dual セル重心に**双対 CV の面積加重重心** (`dualCentroid`, gmshReader) を使用 (軸上でも R>0)。
  → 軸 CV 0 個・体積比 2.4e4 (cell 1.2e4 と同等) に改善。平面 bump loM node は無回帰 (ro_mean/Ux_mean が M2 と一致)。
  **case/29 Euler 軸対称 (conical, viscMethod=0 真の Euler, 1 次, implicit cfl_pseudo=2)**: cell/node とも残差 PASS (rms_ro ~3e-9)。
  **注**: 初め viscMethod=1(Sutherland)で「visc=0 でも粘性が残り」発散していた→真 Euler は viscMethod=0。
  run: `case/29.bell_vs_conical/run_dual_eul_conical_{cell,node}_m3/`。
  **訂正 (収束場で再評価)**: 当初「Mach 一致で validated」と報告したが、それは **res_0=IC を比較していた誤り**
  (outStepInterval=5000>4000 step で収束場が dump されていなかった)。**収束場 (res_4000) で再比較すると、bulk は一致
  (Mach mean 1.846 vs 1.847) だが node に near-axis オーバーシュート**: P max **6.82MPa (> chamber Pt 4MPa=非物理)** が
  軸-inlet 角 (x=inlet, r=0)、Mach max **4.915 vs cell 4.015** が軸-exit。27 セルが P>4MPa、全て inlet-axis。
  → **CV 重心修正で軸ゼロ体積爆発は解消したが、node-centered の near-axis (r=0 境界角) 処理に残課題**。
  bell・2 次・viscous の前に near-axis 処理 (軸 BC・inlet 角・r 重み) の検討が要る。
- `2026-06-14` — **M3 切り分け (重要): 壁は OK、軸が blocker**。
  - **near-axis 角オーバーシュートの試行は全て発散** (軸ノード roUy=0 強制→roe<0で全域NaN、軸ソース抑制→step51破綻)。
    半径方向圧力ソースは軸 CV にも load-bearing で、安易な対称強制は不可。baseline は収束 (corner 1 セルのみ過大)。
    infra (`mesh.axis_flag_d`/`enforceAxisSymmetry_d`) は残置・既定 off。
  - **bell Euler 軸対称**: conical と同様 cell/node とも PASS、bulk 一致 (Mach mean 1.926 vs 1.926)、node に同じ
    inlet-axis corner オーバーシュート (P 6.82MPa, 1 セル)。run: `run_dual_eul_bell_{cell,node}_m3/`。
  - **node-centered no-slip 粘性壁は planar で OK と判明**: `case/24` channel viscous (平面, 非軸対称) を node で実行 →
    **発散せず cell と追従** (ro_mean 1.178 一致, NaN無し, 両者 still converging)。run: `case/24/run_dual_visc_{cell,node}_m3diag/`。
    → **既存ゴースト壁 + CV 重心で node 粘性壁は動く。弱形式 no-slip 壁 (`boundaryNode_d.cu`) は当面不要**。
    implicit 12000step で well-developed (主残差 roUx: cell 4.0dec / node 4.4dec、roUy は near-zero 横流れで両者 plateau)、
    収束場一致 (ro_mean 1.18064 vs 1.18060、Ux_mean 10.68 vs 10.93 ≈2.3%、NaN無し)。
  - **結論**: case/29 viscous (軸対称) の step5 発散は **壁ではなく軸 (near-axis)** が原因。M3 の残 blocker は
    node-centered の **near-axis (r=0) 処理**であり、これは [[architecture-axisym-axis-singularity]] と同種の難問。
    壁・平面 viscous・Euler (bulk) は node で動く。**次の本丸は near-axis 専用の腰を据えた対策**。
- `2026-06-14` — **陽解法 + 残差射影 R←R−(R·n)n の検証 (ユーザ要請)**。軸の対称フラックスは現状ゴースト/slip
  (`replacePrimalWithDual` が軸半割面を 1 セル境界 plane→ゴースト生成、`slip_d` が鏡像、対流カーネルが 1 次で処理)。
  残差射影 (`zeroAxisRadialResidual`=軸で res_roUy=0=R←R−(R·n)n) は assembleResidual に**常時有効**。
  - **explicit で射影は完璧に機能**: detectNaN ダンプで **軸ノードの Uy=0 厳密** (271/272 finite, max|Uy|=0、SU2 と同じ)。
  - **だが explicit は依然発散**: 死因は roUy でなく **exit-axis 角の ro(密度)** (subsonic: step2, CV x≈0.085/r≈0.0004 で ro NaN
    4 セル; supersonic は exit-wall 角)。**= 射影は必要だが十分でない。corner CV (r≈0 の極小 CV, outlet∩axis) が密度で
    explicit 不安定 (多方程式)**。implicit の inlet-axis corner オーバーシュートと同根。
  - **残課題は corner CV の構造/コヒーレント処理**で確定 (roUy 対称は解決済み)。候補: 角 CV を隣接へ merge して極小 CV を
    作らない / 角専用 dt フロア / SU2 CSymmetryPlane を flux+残差+Jacobian 一体移植。
- `2026-06-14` — **[解決] near-axis corner (implicit Euler)**。長く open だった inlet/outlet-axis 角の P オーバーシュート/
  逆流を **2 修正の併用で解決** (commit ec3372b)。①**軸 (R≤eps) で軸対称ソース OFF** (`axisymmetricSource_d` に
  `axis_flag`、node+axisym のみ)。ソース `res_roUy += P·A_planar` は唯一 r 重みされない項で `source/volume = P/r → ∞`
  (r→0) が発散源、対流項は r_face/r_cell キャンセルで有限 → 軸でソースを切るのが本質 (SU2 の y<EPS 相当、ユーザ提案)。
  ②**境界半割面の重心に真の面積加重重心 (R≥0) を使う** (`dualBnodeCent`)。旧 `node+h·n_out` 便宜は pcy≈0/<0 で入口/出口
  BC の r 重み実効面積を ~0 にし **BC が軸近傍 corner CV に届いていなかった** → 真の重心 (R>0) で入口がコーナーに届き
  P=chamber に。(`axisCentroidShift=1` で CV セル中心に面積加重重心を使うのも維持: r_cell と r_face を一致させ対流項の r を
  キャンセル。) **検証** (case/29 conical & bell, 軸対称 Euler implicit): PASS。corner P=3.99MPa(=chamber, overshoot 0
  セル、旧 6.82e6)、Ux=+51.7(流入, 旧 −798 逆流)、Uy≈0 — **SU2 (P=3.99,Ux=+54.8,Uy=0) と一致**。near-axis<3% 帯
  mean|ΔM|=0.016/max0.049 (旧 max0.95) と中/外帯並み。平面 (bump loM) 無回帰。run:
  `case/29.bell_vs_conical/run_dual_eul_{conical_node_srcoff_keepshift, bell_node_m3, conical_cell_m3}/`。
  **残**: explicit は startup の exit-wall 角で依然脆く発散 (cell explicit も同様、軸でなく supersonic-startup ロバスト性) →
  implicit が実用。viscous (軸対称) を次に検証。
- `2026-06-14` — **viscous (軸対称) node は exit-lip corner で発散 (node 固有・軸とは別問題)**。corner 修正後に
  node 粘性 (laminar, viscMethod=1) implicit を再検証 → **step5-8 で `ro` NaN**。detectNaN ダンプの空間特定で
  **NaN 種は exit 出口リップ角** (x∈[0.073,0.086], 出口 x=0.0862; r∈[0.029,0.033], 壁 rmax=0.0331。最下流 NaN セルは
  (0.0859,0.0328)=出口リップそのもの。NaN の 292/309 が壁近傍 r>0.9rmax、**軸近傍 r<0.05rmax は 0**)。
  → **軸 corner 修正は無関係に効いており (軸では発散していない)、これは別の「出口リップ角 CV」問題**。
  - **node 固有と確定**: 同一メッシュ・同一設定の **cell 粘性は完走・NaN 無し** (`run_dual_visc_conical_cell_m3`,
    step7999 まで rms_ro=2.1e-5)。node だけ exit-lip で爆発。
  - **cfl では治らない (構造問題)**: cfl 2.0→0.3 で発散が step5→112 に遅延するだけで最終的に NaN。
    → startup transient の stiffness でなく、出口リップの極小 CV (vol~2.5e-9) で no-slip 壁ゴースト (Ux=Uy=0) と
    supersonic outlet ゴースト (外挿) が同一 CV に同居し、粘性応力 (超音速→0 の du/dy を極小 CV 幅で評価) が爆発する
    構造的問題。inlet-axis corner と同類の「境界∩境界の極小 corner CV」issue (軸でなく壁∩出口版)。
  - **残課題 (確定)**: 出口リップ等の boundary∩boundary 角 CV のコヒーレント処理。候補 (inlet-axis と共通):
    角 CV を隣接へ merge して極小 CV を作らない / 角専用の粘性フラックス制限・dt フロア / wall BC と outlet BC の
    ゴースト整合 (角では wall 優先で no-slip を満たす)。run: `run_dual_visc_conical_node_{cornerfix(o2),cf_o1,cf_cfl03,
    cf_nan,cf_cfl03_nan}/` (全て発散、NaN 種=exit-lip)、cell 比較 `run_dual_visc_conical_cell_m3/` (収束)。
- `2026-06-14` — **node 内部双対面の面補間係数を中点 (fx=0.5) に固定するフラグ `nodeMidpointFx` (既定 0=OFF)**
  (レビュー指摘)。`dualFaceCent` がエッジ中点でなく面積加重重心のため `fx` が歪み面で 0.5 からずれていた問題を、
  標準的な median-dual エッジ補間 $\phi_f=\tfrac12(\phi_A+\phi_B)$ に揃える。`calcStructualVariables_d` に `nodeMode`
  (`discretization=="node" && nodeMidpointFx==1`) を渡し内部面 (`ip<nNormalPlanes`) で `fx=0.5` 固定。
  cell モード・境界半割面は不変。methods/discretization.md §7.4。
  **検証** (node 層流 40k, `run_divufix_node_after` 幾何fx vs `run_fx_node` fx0.5、build-lsq vs build-lsqfx、divu 込み):
  近壁 `dUxdy` checkerboard roughness 99pct **12.84→8.23 (−36%)** 低減。**SU2 (axisym laminar 同条件 `run_su2cmp_su2_lam`)
  壁圧比較で fx ON/OFF は区別不能 (<0.5%)、両者 SU2 に平均 3.2% 一致** → fx=0.5 は SU2 一致を悪化させず安全。
  局所 Ux は出口リップ近傍で 1376→430 と変わるが圧力場に伝播せず、SU2/forge とも未収束のため near-wall 速度細部の
  厳密な是非は未確定 → 既定 OFF で opt-in 残置。node SST は node-viscous 脆弱性 (§7.2)+SST 剛性で baseline build から
  発散し検証不可。
- `2026-07-19` — **親計画としてクローズ (status=done, accepted へ移動, ユーザー承認)**。目的の
  cell/node 両対応化は達成: 同一バイナリ・同一メッシュで cell/node が走り、Euler/粘性/軸対称/SST/
  periodic まで検証済み (node は SU2/cell と bulk 一致、tet で DOF ~1/2)。**現役の入口は子 plan**:
  3D/periodic = [`discretization-median-dual-3d`](discretization-median-dual-3d.md)、境界 ghostless =
  [`discretization-node-boundary-ghostless`](discretization-node-boundary-ghostless.md)、centCoords 統一 =
  [`architecture-node-centroid-value-position`](architecture-node-centroid-value-position.md)、LSQ 勾配 =
  [`discretization-lsq-gradient`](discretization-lsq-gradient.md)。**未消化の検証 1 件 — case/29 node
  viscous (軸対称) の exit-lip 発散が ghostless 5a/5b/5e で解消したかの再検証 run — は ghostless plan
  側へ残課題として移管** (同 plan 変更ログ 2026-07-19 追補)。本文の M2 コスト比較・near-axis/corner
  解決の経緯は歴史記録として保持。

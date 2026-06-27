# node-centered 境界のゴースト撤廃 (段階導入: まず壁, 次に流入出/slip)

## メタ

- **area**: `boundary` / `architecture`
- **status**: `in_progress`
- **related_docs**:
  - `methods/discretization/theory.md` (§6.3 壁 Dirichlet, §6.4 コーナー所有優先)
  - `methods/discretization/implementation.md` (§7.2 ゴースト全撤廃へ)
- **related_plans**:
  - 親: [`discretization-median-dual.md`](discretization-median-dual.md) (node 化本計画)
  - 関連: [`architecture-axisym-axis-singularity.md`](../accepted/architecture-axisym-axis-singularity.md) (軸 corner, 残差射影の先例)
- **created**: `2026-06-14`
- **owner**: `CFD Dev`

## 1. 目的

node-centered (median-dual) の境界はcell-centered 用の**ミラーゴースト方式**を流用しており、node では
幾何学的に不整合 (境界ノードが境界上に乗るのにゴースト鏡像を置く)。これが ① 壁ノード速度が厳密 0 にならない、
② 壁∩出口コーナーで矛盾ゴースト 2 個 → exit-lip 発散 (case/29 viscous node の step5-8 NaN) を生む。
**最終的に node モードの境界ゴーストを全撤廃**し、境界ノード/半割面に BC を直接課す。完了時、node viscous
(case/29) が cell/SU2 と整合して収束し、壁ノード速度が厳密 0 になる。

## 2. スコープ

- **やる (Phase 1)**: 壁ゴースト撤廃 ＋ 壁ノード Dirichlet (u=0) ＋ 壁優先コーナー所有。流入出/slip/axis は
  当面ゴースト維持。断熱壁 (`wall`) のみ (等温壁の熱流束項は拡張)。
- **やる (Phase 2)**: 残り (inlet/outlet/slip/axis) のゴーストも半割面弱形式へ置換し node を完全 ghost-free に。
- **やらない**: cell モードの変更 (全工程で無変更・ビット一致を保証)。3D 混合要素の境界 (別途)。等温壁熱流束
  (case/29 は断熱で不要、必要時 Phase 1 拡張)。

## 3. 関連 docs と前提

- 理論: [theory.md](../../methods/discretization/theory.md) §6 (ゴースト vs 弱形式, §6.3 壁 Dirichlet, §6.4 コーナー所有)。
- 実装方針: [implementation.md](../../methods/discretization/implementation.md) §7 / §7.2。
- 先例: 軸対称 near-axis の**残差射影** (`zeroAxisRadialResidual_d`) が block-DPLUR と整合した
  ([architecture-axisym-axis-singularity.md](../accepted/architecture-axisym-axis-singularity.md))。**state 直書きは Mach1000 発散**
  実績があるため、壁 no-slip も残差射影で実装する。
- 前提: 境界ノードを優先度で 1 bcond に所有させる。種別は `readBcondConfig` が前処理前に `bcond.bcondKind` を
  埋める ([convertGmshToForge.cpp:35,45]) ので利用可。

## 4. 設計方針

**なぜ壁ゴーストを消してよいか** (theory §6.3): 壁ノードを u=0 (Dirichlet) にすると、断熱壁の壁半割面フラックスは
恒等的に 0 (u·n=0 で質量・エネルギー 0、運動量は Dirichlet が吸収、断熱で熱流束 0)。近壁の内部ノードが感じる
壁せん断は、壁ノード (u=0) と内部ノードを繋ぐ**内部双対面**の粘性フラックスが担う。よって壁ゴースト/壁半割面の
flux 寄与は不要。

**コーナー処理 (マルチマーカ emit, 2026-06-22 改訂)** (theory §6.4): 当初は優先度 (wall>inlet>outlet>slip/axis) で
コーナーを 1 bcond に所有させる方式だったが、これは**非所有側の半割面 flux が欠落しコーナー CV が未閉**になる
(壁所有コーナーで出口流出が消える→出口で壁面速度が 0 にならない等の症状)。改めて **SU2 流のマルチマーカ emit
(`ow=ib`)** を採用: 各 bcond は自分の境界面に属するノードの半割面を emit し、コーナーは incident する全 bcond から
半割面 ghost を 1 枚ずつ得る (壁側 mirror ghost + 出口側 ghost)。両半割面ベクトル総和でコーナー CV が閉じ、出口/
入口 BC がコーナーに正しく届く。壁 no-slip はコーナー含め `wall_flag` (壁 iCells) + `enforceWallNoSlip`/
`zeroWallMomentumResidual` の Dirichlet が u=0 を厳密化する (壁優先は所有でなく Dirichlet の上書きで実現)。

**実装上の Dirichlet** = state を一度 0 初期化 ＋ 毎反復で運動量残差を射影してゼロ (陰解法整合)。

## 5. 実装ステップ

1. **マルチマーカ半割面 emit** — [gmshReader.hpp](../../solver_density_cuda/mesh/gmshReader.hpp)
   (`buildMedianDual`/`replacePrimalWithDual`)。`ow=ib` で各 bcond は自分の境界面に属するノードの半割面を emit する
   (全 bcond が ghost を持つ。壁も mirror ghost を持つ)。コーナーは incident する全 bcond の iCells/iPlanes に重複出現し、
   それぞれの半割面 ghost を得る。wall bcond の iCells (コーナー含む) は h5 に残し solver の wall_flag 構築に使う。
   〔旧方式 (優先度所有・壁 plane 不出力) はコーナー CV 未閉のため廃止。〕
2. **wall_flag_d** — [mesh.cpp](../../solver_density_cuda/mesh/mesh.cpp) `setMeshMap_d` /
   [mesh.hpp](../../solver_density_cuda/mesh/mesh.hpp)。`axis_flag_d` と同パターンで wall 種別 bcond の iCells を 1 に。
   空 plane bcond でも readMesh が落ちないことを確認。
3. **壁 no-slip カーネル** — [axisymmetricSource_d.cu](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)/`.cuh`
   に `enforceWallNoSlip_d` (一度きり state 初期化: KE 除去後 roU/U=0) と `zeroWallMomentumResidual_d`
   (毎反復 res_roU{x,y,z}=0)。両 wrapper を `discretization=="node" && wall_flag_d!=nullptr` でゲート。
4. **ディスパッチ/呼び出し** — [boundaryCond.cpp](../../solver_density_cuda/boundaryCond.cpp) で node 時 wall/
   wall_isothermal を no-op 化。[main.cpp](../../solver_density_cuda/main.cpp) で `enforceWallNoSlip` を IC 後 1 回、
   `zeroWallMomentumResidual` を `assembleResidual` 末尾 (軸射影の後) に毎反復。
5. **(Phase 2 — inlet/outlet/slip 弱形式化 + 壁粘性 ghostless 化, 2026-06-22 設計)**

   **背景 (バグ実証)**: node モードの出口は依然 ghost+主ループ処理で、境界ノードが境界上に乗る退化幾何
   (`d_along_n=0`、節点=半割面重心) のため、急勾配域 (壁近傍) で ghost flux が破綻する。`case/26.flat_plate_sst`
   層流 node (`run_node_lam_cont`) で実証: BL は x≤0.94 では Blasius どおりだが **x≥0.98 で壁せん断が消え、出口列
   (x=1.0) で BL 完全崩壊** (第1オフ壁 Ux=32→68=自由流、dUxdy が 9e3→7.8e6 に爆発)。出口自由流圧 (97250) は
   正しく届くのでバルクは健全、壊れるのは退化幾何×急勾配の近壁のみ。`ow=ib` のコーナー修正とは独立の既存問題。

   **方針**: 壁と同様に inlet/outlet/slip も ghostless 弱形式にし、node 主ループを内部面のみに限定する。

   - **5a. 主対流ループを内部(+periodic)のみに**: 現行 node は `convPlaneBound = nNormal_halo_Planes - nWallHaloPlanes`
     で壁 plane だけ除外し、inlet/outlet/slip の ghost plane は主ループに残っている (これが崩壊の原因)。非 periodic
     境界 plane 総数 `nBoundaryHaloPlanes` を `setMeshMap_d` で数え (plane 並びは [内部][periodic][非壁境界][壁]
     なので末尾連続)、node では `convPlaneBound = nNormal_halo_Planes - nBoundaryHaloPlanes` (=内部+periodic) とする。
     viscous 主ループは既に `ip<nNormalPlanes` で内部限定なので変更不要。
   - **5b. 全境界の弱形式対流**: [convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu)
     `convectiveFlux_boundary_d` (bvar を R 状態とする単純 FVS、pRef 控除済) を node では**壁だけでなく全非 periodic
     bcond**で起動する (現状 `nodeMode && !isWallKind` で非壁を除外している分岐を撤去)。bvar は既に各 BC kernel が
     設定済み (inlet=固定状態, outlet=`Uxb`/`Psb`, slip=`Uxb=0`+`Psb=P[ic]`→ mdot=0 pressure-only)。退化 ghost を
     一切使わずコーナーも各 incident bcond の半割面ごとに弱形式で閉じる。
   - **5c. 壁→隣接内部ノード隣接マップ + node 版 viscousFlux_wall_d** (ユーザ設計):
     - **マップ構築** (`setMeshMap_d`): 内部双対面 (`ip<nNormalPlanes`) を走査し、片側が壁ノード (`wall_flag`)・他方が
       非壁ノードの面を「壁-流体面」として CSR (`wallAdj_offset[壁ノード]`, `wallAdj_face[]`, `wallAdj_other[]`,
       面の向き符号) に収集する。**1 壁ノードが複数の内部隣接を持つ**ので CSR で集約する。
     - **新カーネル** `viscousFlux_wall_node_d`: 壁ノードごとにこのマップを回し、各隣接内部ノード M との
       **コンパクト法線差分** `μ (u_t,M - 0)/d` で壁せん断を計算。**粘性 flux は流体側ノード M にのみ加え
       (`-flux`)、壁ノード N 自身の運動量 flux は 0 にする** (N は Dirichlet=運動量方程式を解かないのでペア保存は
       不要・後段 `zeroWallMomentumResidual` とも整合)。τ_w は面積加重で**壁ノードに集約** (utau/y⁺ も)。退化 ghost も
       壁ノードの平滑化勾配も使わない (これが ④ の u_τ 過小=∇u·S 平滑化を解消する)。粘性発熱・isothermal 熱流束も
       流体側に扱う。
     - **二重計上の回避**: 壁-流体面 (N↔M) は内部 viscous 主ループも処理しているため、**内部 viscous ループは
       壁-流体面 (片側が壁ノードの内部面) をスキップ**し、当該面の流体側 M への粘性寄与は `viscousFlux_wall_node_d`
       が一手に担う (壁面を壁カーネルが単独所有)。これで M は壁せん断を 1 回だけ受け、N は運動量 flux 0 で確定する。
     - 旧 `viscousFlux_wall_d` (ghost+`bplane_cell_ghst` 前提) は node 経路から外し、cell モード専用に残す。
   - **5d. 残差順序の確認**: 既に `assembleResidual` は updateVariablesInner→`enforceWallNoSlip`(Dirichlet state)→
     境界 flux で res 設定→`zeroAxisRadialResidual`→`zeroWallMomentumResidual`(運動量のみ)→(blockDPLURSolve+
     applyBlockImplicitCorrection=update) の順。ユーザ提案の順序と一致しており変更不要 (確認のみ)。
   - **5e. block-DPLUR Jacobian の境界半割面スキップ (実装で判明した真因・必須)** —
     [timeIntegration_d.cu](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) `implicit_defect_correction_block_d`。
     5a+5b だけでは出口 BL 崩壊が直らなかった。調査で **出口境界ノードが implicit で凍結 (dq≈0)** しており、それは
     block-DPLUR 対角組立が `cell_planes` (境界半割面を含む) を回り、境界半割面 (`has_nbr=false`) でも粘性対角
     `viscous_diag = 2ν·delta/dcc` (`delta = dcc·ss²/|dcc·S|`) を課すため。**node は境界ノードが境界上に乗り
     `dcc∥境界・dcc·S≈0` で delta が爆発→対角巨大→dq≈0→境界ノード凍結** (= 出口列 BL 崩壊 = 残差プラトーの正体)。
     これは残差側の 5a+5b とは別経路なので 5a+5b だけでは効かない。**node では境界半割面の viscous_diag を
     スキップ** (`int isNode` 追加、`!(isNode && !has_nbr)` でゲート)。対流 A⁺S 対角は弱形式流出 Jacobian の近似
     として残す。cell は境界ゴーストが法線方向に非退化なので従来どおり (isNode=0、ビット不変)。

   **段階規律 (codex 指摘)**: **5a+5b は「境界対流の責務変更」なので独立に小検証して一度 commit で区切る**。5c の壁粘性
   まで一気に入れると「出口 BL 崩壊の改善」と「u_τ 補正」の効果が混ざり原因切り分けが困難になる。順序は
   **5a+5b → 検証/commit → 5c → 検証/commit**。

   **5b の注意 (codex 指摘)**: 既存 `convectiveFlux_boundary_d` は「壁 pressure-only」前提の符号・面法線・左右状態に
   なっていないかを丁寧に確認する。入口で bvar を R 状態にしたとき SLAU/ROE の流入 flux と既存 ghost 方式が
   どの程度一致するかを、まず **1 次精度 (convMethod=0) の短い run** で確認してから本検証に進む。

## 6. 検証

- **ビルド**: native `solver_density_cuda/build-dual` (RTX3060 arch86)。cell 回帰 (`case/08.bump`) がビット一致。
- **検証ケース**: [case/29.bell_vs_conical](../../case/29.bell_vs_conical) の node viscous conical/bell を**新規 run_***
  (`run_dual_visc_conical_node_wallfix/` 等) で。`tools/check_convergence.py` VERDICT=PASS。
- **判定基準**: step5-8 で落ちず完走、全残差列 (rms_ro/roUx/roUy/roUz/roe) 低下、NaN/Inf=0、`res_*.h5` で
  P≤Pt・ro>0・T>0、**壁ノードで Ux=Uy=Uz=0 を数値確認**。cell viscous (`run_dual_visc_conical_cell_m3`) と場比較。
  Euler conical/bell node (corner 修正済み) と平面 bump loM node が無回帰。`residual_history.png` 生成。

### Phase 2 専用の検証 (codex 追加 3 項 + 段階)

- **検証ケース**: `case/26.flat_plate_sst` の **新規 run** (層流 `run_node_lam_*`、SST `run_node_sst_*`)。`ow=ib` 修正済み
  メッシュを fresh convert し、収束場から restart。cell モードはビット不変 (全変更を node ゲート、`case/08.bump` で確認)。
- **(検証A, 5a+5b 後) 主ループに境界 ghost plane が残っていないことを数える**: 起動ログに `node conv main planes =
  internal + periodic` と `boundary weak planes = Σ(inlet/outlet/slip/wall)` の個数を出し、主ループ面数 = nNormalPlanes
  (+periodic) になっていること、弱形式面数 = 全境界半割面数になっていることを確認する。
- **(検証B, 5b 後) 境界ごとの mass flux / pressure flux 積分を出す**: 特に **outlet で `Psb=97250` が ghost/bvar だけ
  でなく実際の flux に効いている**こと (出口圧が場に伝播し自由流が加速する) を確認。まず convMethod=0 の短い run で
  既存 ghost 方式との一致度を見る。
- **(検証C, 5a+5b の主目的) 出口 BL 崩壊の解消**: δ99(x) が x で単調増加 (x≥0.98 で崩れない)、第1オフ壁 Ux と dUxdy が
  出口列で爆発しない (層流 run で `run_node_lam_cont` の dUxdy 9e3→7.8e6 爆発が消えること)。
- **(検証D, 5c 前後) 壁せん断精度を cell と比較**: `u_τ`, `Cf`, `y⁺`, wall shear 積分を cell baseline と並べる。
  既知の `u_τ node 1.24 vs cell 1.97` がどこまで近づいたかを数値で示す。**5c は 5a+5b の commit 後に単独で入れ、
  効果 (BL 崩壊解消 vs u_τ 補正) を分離する**。

## 7. 影響範囲

- ファイル: `mesh/gmshReader.hpp`, `mesh/mesh.cpp`, `mesh/mesh.hpp`,
  `cuda_forge/axisymmetricSource_d.cu`/`.cuh`, `boundaryCond.cpp`, `main.cpp`。(Phase 2 で convertGmshToForge/
  境界 flux 経路一式)。
- 既存ケース: cell モードは無変更・ビット一致。node モードの壁挙動が変わる (改善)。
- docs: `methods/discretization/theory.md` §6.3/§6.4, `implementation.md` §7.2 (更新済み)、`methods/index.md` は項目既存。

## 8. 完了条件

- [x] 関連 `methods/discretization/theory.md` 更新済み
- [x] 関連 `methods/discretization/implementation.md` 更新済み
- [ ] 実装・検証完了 (§6 を満たす)
- [ ] `plans/README.md` の状態を `done` に更新
- [ ] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-06-14` — 初稿。node 境界ゴースト全撤廃の段階計画 (Phase 1: 壁 Dirichlet＋壁優先コーナー, Phase 2: 残り撤廃)。
  親 plan discretization-median-dual から分離。
- `2026-06-14` — **方式確定・viscous node が回る (explicit)**。試行錯誤の結論:
  - 残差射影 Dirichlet は壁圧力寄与を落とし発散 (不正、ユーザ指摘)。壁 ghost 完全撤廃は勾配閉性を壊し悪化。
    → 採用: **壁優先コーナー所有 (commit 93ef041) + 壁の弱形式 pressure-only 対流 (bvar 速度=0, mdot=0→圧力のみ)**。
    壁 ghost は勾配/粘性 mirror 用に残す。inlet/outlet/axis は実績ある主ループ+ghost のまま (全境界弱形式は
    `convectiveFlux_boundary_d` の inflow flux が SLAU 非等価で inlet-axis corner を 6.1MPa に悪化させたため壁限定)。
  - 実装: `wall_d` bvar=no-slip(0), `normal_halo_planes_d` 末尾に壁 plane (`nWallHaloPlanes`), 主対流ループは
    node で壁除外, `convectiveFlux_boundary_d` を壁 bcond のみ起動。
  - **viscous は explicit**: 近壁極小双対 CV の粘性が block-DPLUR 近似対角で implicit 化不足 (cfl2→0.1 で発散
    step13→610=構造的)。explicit(RK3) は setDT 粘性スペクトル半径で局所 dt が縮み安定。cfl0.5 発散・**cfl0.1 安定**。
  - **検証**: node Euler conical 無回帰 (Pmax=4.475MPa=committed 一致, inlet-axis overshoot なし)。cell viscous conical
    無回帰 (rms_ro 2.02e-5 vs baseline 2.10e-5)。**node viscous conical explicit cfl0.1: 30k step 完走・NaN なし・残差
    ~3 桁低下 (still converging)・場 物理的** (P≤Pt, ro>0, T 359-2231K)。run:
    `case/29.bell_vs_conical/run_expl_test/` (residual_history.png), 長 run `run_dual_visc_conical_node_expl_long/`。
  - **残**: 陰解法 viscous node (近壁粘性 Jacobian 強化=gpu-implicit-plan 後続)、Phase 2 (inlet/outlet/slip/axis +
    スカラ輸送のゴースト撤廃)、bell viscous explicit, より深い収束。
- `2026-06-22` — **コーナー所有を優先度方式からマルチマーカ emit (`ow=ib`) に変更**。優先度所有は壁所有コーナーで
  出口半割面 flux が欠落しコーナー CV が未閉だった (出口で壁面速度が 0 にならない症状)。`gmshReader buildMedianDual` で
  各 bcond が自境界面ノードの半割面を emit し、コーナーは全 incident bcond から半割面 ghost を 1 枚ずつ得る形に。
  dead code (`ownerBc`/`bcondPriority`/`isWallKind` 優先度ロジック) を撤去、`gmshReader.hpp`/`mesh.cpp` の旧コメント
  も是正。検証 `case/26.flat_plate_sst/run_node_corner_verify/` (fresh convert→restart 500 step): closure 3e-6、
  コーナー node 2 (x=1,y=0) Ux=Uy=0、全壁ノード \|U\|=0、NaN なし、コーナーが wall+outlet 両 iCells に出現。
  〔別問題として残差プラトー rms_roUx≈0.23 (checkerboard 疑い) は未解決〕。
- `2026-06-22` — **Phase 2 着手 (設計確定)**。出口 ghost+主ループの退化幾何で近壁 BL が崩壊するバグを
  `case/26.flat_plate_sst/run_node_lam_cont` で実証 (x≥0.98 で壁せん断消失→出口列 BL 崩壊、dUxdy 9e3→7.8e6)。
  §5 Phase 2 に設計を記載: (5a) node 主対流ループを内部+periodic のみに (`nBoundaryHaloPlanes` 導入)、(5b) 全非
  periodic 境界を `convectiveFlux_boundary_d` の bvar 弱形式に、(5c) **壁ノード→隣接内部ノード CSR 隣接マップ +
  node 版 `viscousFlux_wall_node_d`** (コンパクト法線差分で τ_w 集約、流体側 M のみに flux・壁ノード運動量 flux=0、
  内部 viscous ループは壁-流体面をスキップ)、(5d) 残差順序は現行で要件充足 (確認のみ)。勾配は既に ghostless
  (`calcGradient_cellgather_d` の境界スキップ + `calcGradient_b_d`/LSQ bvar) なので Phase 2 対象外。実装はこれから。
- `2026-06-22` — **5a+5b+5e 実装・検証完了 (出口 BL 崩壊 + 残差プラトーを解消)**。5a+5b だけでは効かず、真因が
  **block-DPLUR Jacobian の境界半割面・退化粘性対角による境界ノード凍結** と判明 → **5e** (node で境界半割面の
  `viscous_diag` をスキップ) を追加。検証 (`case/26.flat_plate_sst`):
  - 層流 `run_node_lam_5e_long` (40000 step, 旧崩壊場から restart): 出口 node242 Ux 32.27→0.058、dUxdy 7.8e6→1.4e4、
    **δ99(x) 単調で Blasius 一致** (x=1.0 δ99=0.0114, 旧 1e-5 崩壊)。**rms_roUx 0.214→3.08e-5 (3.8桁)**, rms_ro 1.04e-7。NaN なし。
  - SST `run_node_sst_5e` (10000 step): rms_roUx 0.597→2.14e-3 (2.4桁 falling, 旧 0.23 STALLED)、mean flow 全列 falling。
    roK 微増は壁せん断精度 (5c) 課題。
  - cell: 全変更 node ゲート (convPlaneBound cell 分岐不変・診断 log node 限定・isNode=0 で viscous_diag 従来通り)、
    `case/08.bump` cell SLAU が NaN なし正常。検証A: `conv main planes=42691(internal 42691+periodic 0), boundary weak
    planes=659(wall 200+non-wall 459)` を起動ログ確認。
  - **残**: 5c (壁せん断 u_τ 精度)、SST roK 収束、検証B (境界別 flux 積分)。次フェーズ。
- `2026-06-22` — **検証D (u_τ) で 5c の前提が崩れる重要発見**。5a+5b+5e 後の SST (`run_node_sst_5e_long`, 40000 step)
  で u_τ を cell baseline (run_0001, **u_τ≈1.97**) と比較: node u_τ≈**1.33–1.43** (旧 1.24 から微改善も依然 ~30% 過小)。
  原因を切り分けたところ **5c (viscous flux 離散化) ではなく SST 壁での乱流粘性の非減衰**だった:
  - node 壁ノードで **k=66.8 (本来 0)**, omega=4.2e5 (`wall_y_eff`=6.1e-5 由来で有限), **vis_turb/μ≈10** (cell は近壁で 0)。
  - `rans_wall_scalar_boundary_d` (wall-resolved) は `kb=0` と**ゴースト** `k[ig]=2kb−k[ic]` を設定するが、node モード
    では壁ノードが CV 中心 `ic` で、**`k[ic]=0` をピン留めする機構が無い** (速度は `enforceWallNoSlip` で u=0 を
    ピン留めするが k には無い)。→ 壁ノード k が 67 に発達し near-wall vis_turb が過大 → BL が過拡散し u_τ 過小。
  - **結論: 5c は u_τ を直さない。正しい修正は「node モードの SST k=0 壁 Dirichlet」** (enforceWallNoSlip の k 版:
    壁ノードで k 状態=0 ＋ roK 残差=0、必要なら omega も壁ノードでピン)。handoff の ④ (∇u·S 平滑化) 仮説は外れ。
  - 検証データ: `case/26.flat_plate_sst/run_node_sst_5e_long/` (res_40000), cell 比較 `run_0001_slau_rans_implicit/`。
- `2026-06-22` — **上記 u_τ 診断を再訂正 (有効な乱流基準で確認)**。`run_0001` は Ps=99303 の**低速ほぼ層流** (U=34,
  peak mut/μ≈1) で基準にならなかった。同条件の有効基準 = **`run_0009_ewt_fine_mode1`** (cell, Ps=97250, U≈67,
  automatic WT, peak mut/μ=199, **Cf≈Schlichting (0.94–1.0)**, k_wall≈0)。これと node を突き合わせると:
  - **node はむしろ過剰乱流**: peak mut/μ=375 (cell 199), **Cf ≈ 3× Schlichting**, δ99 ≈ 1.5–2× 過厚。
  - 壁ノード k=67 は**過生成された内部 k をノイマンで追随しただけの症状**で、根本ではない (ユーザ指摘: k は
    automatic WT でノイマン=処理不要、ω こそ壁で値固定+residual0+陰解 decouple)。
  - **試した壁 BC 修正は無効/悪化 → revert 済み**: (wallTreatmentSST=1 + ω point-implicit decouple dω=0 +
    `wall_y_eff` を MEAN→MIN で ω 値是正)。過渡では壁 vis_turb≈0 になったが収束させると過剰乱流が再発・悪化
    (peak mut/μ=909, Cf 3.5×)。→ 根本は**壁 BC でなく BL 全体の k 過生成** (node 固有, おそらく歪み速度/生成項 or
    convMethod=0)。別途の SST 生成項調査が必要 (本 plan の範囲外, 次フェーズ)。
  - **確定済み (commit 1ac6a57/50867a5) の 5a+5b+5e=mean-flow の出口 BL 崩壊+残差プラトー解消はこの乱流件と独立**。
- `2026-06-22` — **過剰乱流を解消 (ユーザ+レビューAI の指摘が的中)**。「k 過生成」と片付けたのは誤りで、真因は
  **node モードで壁 ω が Dirichlet ピンされていない**ことだった (ghost BC `omega[ig]=2ω_w-omega[ic]` は壁ノードが
  CV 中心・dcc≈0 退化で omega[ic] を ω_w に固定できない → ω 過小 → k 消散 Dk=β*ρkω 不足 → 過生成 runaway。壁 k=67 は
  その症状で k=0 Dirichlet は誤り)。修正 (node 限定, commit 2b19d1d / massflux は 0f6d53d):
  - **出口で massflux=0 バグ** (5a/5b で `convectiveFlux_boundary_d` が massflux 未書込→スカラが境界で移流せず k が
    出口で蓄積) を修正: 境界 mdot を massflux に書き戻す。
  - `compute_wall_y_eff` を MEAN→**MIN** (最近接オフ壁ノード距離=正しい Δy。MEAN は Δy 過大評価で ω_w ~200x 過小)。
  - `rans_wall_scalar_boundary_d` wall-resolved 分岐で node は **omega[ic]=ω_w・roOmega[ic]=ρω_w を直接ピン** (保存変数)。
  - `applySSTPointImplicit` で壁 ω 行を **decouple (dω=0)**、k はノイマンで通常更新。
  - `zeroWallMomentumResidual` で **res_roOmega も壁で 0** に (Dirichlet 残差で rms 汚染を除去)。
  - 結果 (`run_node_sst_omegapin2`, wallTreatmentSST=0): **Cf/Schlichting 1.03/1.21/1.54 (旧 ~3.0–3.4)**、壁 mut/μ≈0
    (旧 10)、peak mut/μ 367 (旧 375–470, cell 199)、δ99 上流で cell 近接。rms_roOmega 1.1e6→0.14。x=0.89 が依然 ~1.5x
    (出口域・peak mut/μ がやや高い) は follow-up。
  - 〔反省: 「壁 BC でない/k 過生成」と早合点したのは、保存変数 roOmega の Dirichlet 整合 (ρ 変化後の再ピン) と
    omega 値 (MIN) を同時に満たさず単独で試したため。レビューAI の "omega は保存変数 Dirichlet" 指摘が決定打。〕
- `2026-06-22` — **残課題 (x=0.89 過剰乱流) の真因 = wall_dist バグを修正、SST node が cell 基準と一致 (commit 8b60f2c)**。
  ユーザ指摘「wdist はちゃんと入ってるか」が的中。`calcWallDistance_kdtree` が壁点集合に**壁半割面の重心**を使って
  おり、node では重心が壁ノードから x 方向に ~dx/8 ずれている (off-wall ノードは壁ノードの直上に整列)。近壁ノードの
  最近接壁点距離が法線距離 y でなく x ずれ (≈dx/8, 下流でメッシュ粗化に伴い増大) に支配され、**wall_dist≈1e-4·x
  (壁ノード自身も 0 でなく 8.9e-5、出口で 5e-3)** と大誤り → ω_w=60ν/(β·wall_dist²) 過小 (下流ほど) → 過剰乱流。
  修正: node では壁点に**壁ノード座標 (bc.iCells の CV 中心, 壁面上で off-wall ノードの直下)** を使う。cell は plane
  重心のまま (セル中心の直下に整列, 不変)。検証 (`run_node_sst_wdist`, 再 convert): wall_dist=0 (壁), =y (第1オフ壁)
  が全 x で成立。**SST node が cell 基準 (run_0009) と一致: Cf/Schlichting 0.89/0.91/0.93 (旧 3×), peak mut/μ 202
  (cell 199), k_wall=0, δ99 が cell 一致**。massflux(0f6d53d)+omega pin(2b19d1d)+wall_dist(8b60f2c) の 3 点で完結。

## 変更ログ (追補 2026-06-26)

- node スカラー/species 境界の **ghost 読み残り** を ghostless 化 (BC は ghost を書くが node consumer は
  境界ノード値を使う設計との不整合を解消):
  - **k/ω advection** (`scalar_advection_first_order_d`): inflow で `phi[ic1]`(ghost) でなく境界ノード
    `phi[ic0]` を使う (`isNode` ガード)。入口 Dirichlet ピンと整合し移流が真の入口値を運ぶ。
    cell ビット不変、node SST は入口移流が正しくなり ~1.4e-2 改善 (回帰でなく correctness)。
  - **species 拡散** (`species_diffusion_d`): 境界半割面で ghost mirror の `dcc≈0` が 0/0 退化するのを、
    境界ノード species 勾配の弱形式 `J_s=ρD(∇Y_s·S)` に置換 (k/ω 拡散と同型)。`speciesGradient` を
    node+粘性多成分でも計算するよう gating 拡張。
  - **species faceY advection** (`species_advection_faceY_d`): node 境界半割面は主ループ除外で `Yface[ip]`
    が stale → 境界ノード組成を使用。
  - flow(ρ,ρu,ρe) は既に bvar 弱形式で閉じ済み・全 BC (slip/wall/inlet/outlet) が bvar 5 量を書込済みで
    node 正しいことを確認 (本監査の coverage)。
  - 注意: species + node は検証ケース未整備のため上記 species 修正は **未検証** (構造は k/ω の実証済み
    パターンに準拠)。node 単成分 SST と cell は回帰確認済み。

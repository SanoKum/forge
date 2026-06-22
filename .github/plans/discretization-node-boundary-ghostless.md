# node-centered 境界のゴースト撤廃 (段階導入: まず壁, 次に流入出/slip)

## メタ

- **area**: `boundary` / `architecture`
- **status**: `in_progress`
- **related_docs**:
  - `docs/discretization/theory.md` (§6.3 壁 Dirichlet, §6.4 コーナー所有優先)
  - `docs/discretization/implementation.md` (§7.2 ゴースト全撤廃へ)
- **related_plans**:
  - 親: [`discretization-median-dual.md`](discretization-median-dual.md) (node 化本計画)
  - 関連: [`architecture-axisym-axis-singularity.md`](architecture-axisym-axis-singularity.md) (軸 corner, 残差射影の先例)
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

- 理論: [theory.md](../../docs/discretization/theory.md) §6 (ゴースト vs 弱形式, §6.3 壁 Dirichlet, §6.4 コーナー所有)。
- 実装方針: [implementation.md](../../docs/discretization/implementation.md) §7 / §7.2。
- 先例: 軸対称 near-axis の**残差射影** (`zeroAxisRadialResidual_d`) が block-DPLUR と整合した
  ([architecture-axisym-axis-singularity.md](architecture-axisym-axis-singularity.md))。**state 直書きは Mach1000 発散**
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
5. **(Phase 2)** 半割面ジオメトリの device 転送＋`boundaryFluxNode_d`/`boundaryGradientNode_d`＋共有 flux ループの
   境界 plane スキップ。

## 6. 検証

- **ビルド**: native `solver_density_cuda/build-dual` (RTX3060 arch86)。cell 回帰 (`case/08.bump`) がビット一致。
- **検証ケース**: [case/29.bell_vs_conical](../../case/29.bell_vs_conical) の node viscous conical/bell を**新規 run_***
  (`run_dual_visc_conical_node_wallfix/` 等) で。`tools/check_convergence.py` VERDICT=PASS。
- **判定基準**: step5-8 で落ちず完走、全残差列 (rms_ro/roUx/roUy/roUz/roe) 低下、NaN/Inf=0、`res_*.h5` で
  P≤Pt・ro>0・T>0、**壁ノードで Ux=Uy=Uz=0 を数値確認**。cell viscous (`run_dual_visc_conical_cell_m3`) と場比較。
  Euler conical/bell node (corner 修正済み) と平面 bump loM node が無回帰。`residual_history.png` 生成。

## 7. 影響範囲

- ファイル: `mesh/gmshReader.hpp`, `mesh/mesh.cpp`, `mesh/mesh.hpp`,
  `cuda_forge/axisymmetricSource_d.cu`/`.cuh`, `boundaryCond.cpp`, `main.cpp`。(Phase 2 で convertGmshToForge/
  境界 flux 経路一式)。
- 既存ケース: cell モードは無変更・ビット一致。node モードの壁挙動が変わる (改善)。
- docs: `docs/discretization/theory.md` §6.3/§6.4, `implementation.md` §7.2 (更新済み)、`docs/index.md` は項目既存。

## 8. 完了条件

- [x] 関連 `docs/discretization/theory.md` 更新済み
- [x] 関連 `docs/discretization/implementation.md` 更新済み
- [ ] 実装・検証完了 (§6 を満たす)
- [ ] `.github/plans/README.md` の状態を `done` に更新
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

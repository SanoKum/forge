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

**コーナー所有優先** (theory §6.4): 壁∩出口のリップは物理的に壁点。優先度 wall>inlet>outlet>slip/axis で壁が所有し、
当該ノードを壁 Dirichlet として扱い出口 bcond から除外 → 矛盾ゴースト消滅。

**実装上の Dirichlet** = state を一度 0 初期化 ＋ 毎反復で運動量残差を射影してゼロ (陰解法整合)。

## 5. 実装ステップ

1. **壁優先所有・壁 plane 不出力** — [gmshReader.hpp](../../solver_density_cuda/mesh/gmshReader.hpp)
   (`buildMedianDual`/`replacePrimalWithDual`)。各境界ノードを優先度で 1 bcond に所有。wall 所有ノードは境界
   plane を出力しない (ghost 生成されない)。wall bcond の iCells (ノード列) は h5 に残す。非 wall 所有ノードは
   所有 bcond にのみ半割面 plane 出力 (コーナーは出口 plane を持たない)。空 plane bcond を h5 writer が許容するか確認。
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

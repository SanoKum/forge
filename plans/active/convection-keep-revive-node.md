# KEEP スキーム復活 (modern API・cell/node) + node WALE — LES/ILES 構成

## メタ

- **area**: `convection / turbulence`
- **status**: `in-progress`
- **related_docs**: [`methods/convection/implementation.md`](../../methods/convection/implementation.md), [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md)
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 1. 目的

node-centered (median-dual) で **KEEP 対流スキーム + WALE SGS** による LES/ILES を回せるようにする。
KEEP は運動エネルギー・エントロピー保存の低散逸中心流束で、LES の乱流エネルギースペクトルを
損なわない。散逸は SGS (WALE) が物理的に担う、という構成を狙う。

## 2. 出発点 (調査結果)

- **KEEP**: `cfg.solver` は `SLAU/SLAU2/HLLE/ROE` のみ有効。KEEP カーネルは
  `cuda_forge/convection/legacy/convectiveFlux_ausm_keep_d.inc.cuh` の `KEEP_FVS_d` に
  **旧フラット引数 API・`ip<nNormalPlanes` 限定**で退避され dispatch から到達不能。
  さらに**真の KEEP 中心流束 (Ctilde/Mtilde/Ktilde/Itilde) は `if(false)` で無効**、
  単純中心平均 + Roe 行列散逸のみ稼働していた。
- **WALE**: `WALE_d` は現役だが `dimGrid_cell` 起動で **cell 専用**。node では起動されず
  `vis_turb=0` (= LES OFF) になっていた。
- 対流 dispatch (`convectiveFlux_d_wrapper`) は `FaceGeom`+`convPlaneBound` で cell/node を
  吸収するので、KEEP を modern API で復活させれば node でも自動で効く。

## 3. 設計・実装

### 3.1 KEEP_d (新規 `convectiveFlux_keep_d.inc.cuh`)
- `KEEP_FVS_d` を modern bundled API (`FaceGeom/PrimState/ResidualOut/LimiterFields/GradFields`)
  へ移植。シグネチャは `ROE_d` に倣う (`conv_scheme, limit_scheme, ga, cnd, geom, st, reso, lim, grd`)。
- **KEEP 中心流束を有効化** (legacy の `if(false)` ブロックを常時実行に変更)。
  - 中心流束は隣接対 `(ic0,ic1)` の**生値**で構成 (KE/エントロピー保存)。
  - 散逸は MUSCL 再構成 L/R 状態の Roe 行列 $R|\Lambda|L$ を L/R で評価し $0.5(|A_L|Q_L-|A_R|Q_R)$。
- 周回面 = `geom.nLoopPlanes` (= `convPlaneBound`)。境界 ghost (`ic1>=nCells`) は `conv_scheme=-1` で 1 次化。
- `massflux[ip]` に散逸込み総質量流束を書く (スカラー輸送整合)。
- dispatch: `convectiveFlux_d.cu` に `else if (cfg.solver=="KEEP")` を追加。cell モードの
  `skipBoundaryFluxKernel` に KEEP を追加 (主ループが境界 ghost を処理するので専用境界カーネルは skip)。

### 3.2 node WALE
- `turbulent_viscosity_d_wrapper` の `WALE_d` 起動グリッドを `dimGrid_cell` → `dimGrid_normalcell`
  に変更 (SST と同じ)。WALE_d 本体は DOF ごとの勾配/体積/wall_dist を読むだけで discretization 非依存。

## 4. スコープ

- **やる**: KEEP_d (cell/node, KEEP 中心流束 + Roe 散逸) の復活と dispatch、node WALE 有効化、
  backstep での node+KEEP+WALE LES 計算。
- **やらない (おいおい対処)**: Roe 散逸の Ducros スケーリング/低マッハ補正による LES 用の散逸低減、
  KEEP_SLAU / AUSM 系の復活、KEEP の陰解法 Jacobian (現状 explicit RK3 のみ)。

## 5. 検証

- backstep 3D node + KEEP + WALE (LESorRANS=1, LESmodel=1) を非定常で回し、発散しないこと・
  剥離せん断層/再付着/spanwise 乱れが発達することを確認する。SLAU ILES (run_node3d_iles_unsteady) と
  剥離構造・再付着長を比較。
- KEEP 中心流束の KE 保存性は一様流 free-stream 保持 + 等方性乱流減衰の挙動で確認 (段階)。

## 6. 影響範囲

- 新規: `cuda_forge/convection/convectiveFlux_keep_d.inc.cuh`。
- 変更: `convectiveFlux_d.cu` (include + dispatch + skip リスト)、`turbulent_viscosity_d.cu` (node WALE グリッド)。
- 既存挙動: `solver != "KEEP"` かつ既存 LES/RANS は**不変** (KEEP は新規分岐、WALE は cell でも normalcell で同一網羅)。

## 7. 変更ログ

- 2026-06-28: KEEP_d 移植 + dispatch、node WALE 有効化、methods 更新。ビルド・backstep 計算で検証中。
- 2026-06-28: **KEEP_d を純粋 KEEP に簡素化** (user 依頼)。Roe 行列散逸・MUSCL 再構成・リミタ・Ducros・
  `keepDissipation` 切替を撤去し、引数を `KEEP_d(ga, geom, st, reso)` に縮約。`space.keepDissipation`
  config も廃止 (残っていても無視)。検証: Taylor-Green M0.4 で cell・node とも運動量~1e-7・KE0.4%・
  エントロピー~1e-5 保存 (`case/09`)、回帰 smoke PASS。**注意**: 低マッハ checkerboard 抑制の散逸が無くなったため、
  homogeneous 方向の市松は `lowMachPrecond`/SLAU 併用で対処する (旧 keepDissipation=1 経路は廃止)。

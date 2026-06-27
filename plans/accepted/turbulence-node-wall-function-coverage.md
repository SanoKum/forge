# node-centered SST 壁関数の近壁カバレッジ修正 (第一内層ノードへ適用)

## メタ

- **area**: `boundary`
- **status**: `done`
- **related_docs**:
  - `methods/turbulence/theory.md` (§6.5 automatic wall treatment)
  - `methods/turbulence/implementation.md` (§3.7 automatic wall treatment 責務分担)
- **related_plans**:
  - 関連: [`discretization-node-wall-implicit-dirichlet.md`](../accepted/discretization-node-wall-implicit-dirichlet.md) (壁運動量3行 decouple の先例)
  - 関連: [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md)
- **created**: `2026-06-28`
- **owner**: (未割当)

## 1. 目的

node-centered (median-dual) SST の automatic wall treatment が、壁関数処理 (P_k 置換 + ω ピン + ω decouple) を **壁ノードにしか適用せず、第一内層ノードを取りこぼす**ため、近壁で k が暴走し過剰乱流になる (`case/18.backstep` 段差凹コーナーで vis_t/l=6789 vs cell 77)。本計画は **node の壁関数処理の適用先を「壁ノード」から「第一内層ノード (Normal_Neighbor)」へ移し**、cell (=壁隣接第一セルに適用、SU2 と x_R 一致で検証済み) と整合させる。完了時、node の近壁 k/ω/vis_turb が cell と整合し、SU2 とも整合する。cell の挙動は不変。

## 2. スコープ

- **やる**:
  - node モードで `wf_pk` (P_k 置換値) の書き込み先を 壁ノード `ic` → 第一内層ノード `irep` (既存の Normal_Neighbor) に変更。
  - node モードで ω ピン (`omega`/`roOmega`) の適用先を `ic` → `irep` に変更 (y は `irep` の `wall_dist` = 壁→第一内層の法線距離)。
  - ω 行 decouple (`applySSTPointImplicit_d` の `dω=0`) の判定を `wall_flag[ic]` → `wf_pk[ic]>=0` に変更し、ピン先 (=irep) に追従させる。
- **やらない**:
  - **cell モードの変更** (現行 = 壁隣接第一セルに適用で SU2 一致・検証済み。ビット不変を保つ)。
  - k を Dirichlet 化する SU2 全面置換 (route 2)。k は forge 設計どおり zero-gradient + P_k 置換で解いたまま (ω のみ Dirichlet)。
  - 運動量壁せん断 (`viscousFlux_wall_d` / `Tau_Wall`) は壁ノードのまま (no-slip 速度・壁応力は壁面の話で本件と別。`Tau_Wall` は壁ノードに保持)。
  - 3D median-dual (M4)・periodic は対象外。

## 3. 関連 docs と前提

- 現在仕様: [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5 (a)-(d)、[`implementation.md`](../../methods/turbulence/implementation.md) §3.7。
- 真因の確定根拠: [`notes/investigations/turbulence-node-wall-function-coverage.md`](../../notes/investigations/turbulence-node-wall-function-coverage.md) (k 収支実測: cell は壁隣接第一セル列が全て wf_pk≥0、node は壁ノード1点のみ・第一内層 wf_pk=-1 で標準生産 Pk=166-228 が暴走)。
- SU2 整合根拠: SU2 `CTurbSSTSolver::SetTurbVars_WF` は Normal_Neighbor (第一内層) に k,ω を設定。forge は ω のみ Dirichlet (k は P_k 置換で解く) だが、**適用先が第一内層**である点を SU2 と揃える。
- 既存インフラ: `compute_wall_friction_sst_d` ([ransWallFunction_d.cu](../../solver_density_cuda/cuda_forge/ransWallFunction_d.cu)) は node で代表内部点 `irep` (SU2 Normal_Neighbor 流, 壁内向き法線 cos 最大) を既に算出し u_τ 計算に使用。本計画はこの `irep` を wf_pk/ω ピンの**書き込み先**にも使う。

## 4. 設計方針

> **注 (2026-06-28 実装で改訂)**: 本節の当初案 (wf_pk と ω ピンの両方を irep へ移設) は race と ω 崩壊で**失敗**した。
> 実際に採用した設計は **§9 変更ログ**参照 = 「P_k 置換のみ壁+第一内層の両ノードへ、ω ピン/decouple は壁ノード据え置き」。

現状 node: `bplane_cell[ib]=ic`=壁ノード。wf_pk[ic]・ω[ic] ピン・decouple(wall_flag) が全て壁ノードに当たり、力学的に重要な第一内層ノードが標準 SST のまま → k 暴走。

変更: node の壁関数処理を `irep` (第一内層ノード) に集約する。

1. **`compute_wall_friction_sst_d` (node)**: 算出済み `irep` を per-bplane に保存 (`irep_b[ib]`、新 bvar)。`wf_pk` の書き込みを `wf_pk[ic]` → `wf_pk[irep]` に (cell は `ic` のまま=第一セル)。`Tau_Wall[ic]` は壁ノードのまま (運動量再スケール用、別件)。
2. **`rans_wall_scalar_boundary_d` (node)**: ω ピン先を `irep_b[ib]` に。y は `wall_dist[irep]` (第一内層ノードの壁距離=Δy、`compute_wall_y_eff` の近似と一致)。cell は `ic` (第一セル) のまま。
3. **`ransSource_d`**: `res_roOmega`/`src_jac_omega` の 0 化は既に `wf_pk>=0` 判定なので、wf_pk が irep に移れば**自動追従**。P_k 置換も同様。変更不要。
4. **`applySSTPointImplicit_d`**: ω decouple 判定を `wall_flag[ic]==1` → `wf_pk[ic]>=0` に変更 (ピン先=irep に追従)。`wf_pk` を引数に追加。cell は従来 decouple OFF + ransSource 0 化で dω=0 を実現しており不変 (cell の `wf_pk>=0` セルでも ransSource 0 化が効くため挙動同じ; 明示 decouple を cell にも効かせるかは実装時にビット不変を確認して決める)。
5. **壁ノード自体**: ω ピンが外れ通常 solved DOF に戻る。第一内層ノードのピンが近傍を抑える想定。挙動を §6 で実測確認 (必要なら壁ノードにも軽い zero-gradient 担保)。

cell は §1-5 で一切変更されない (`discretization=="cell"` 分岐は `ic`=第一セル維持) ことをビット不変試験で担保する。

## 5. 実装ステップ

1. **bvar 追加**: `irep_b` (per wall bplane の Normal_Neighbor node index) を `boundaryCond.hpp` の wall bvar 初期化と `mesh.hpp` の `bplaneValNames` に追加 (片方だと device 未確保で IMA)。型は index だが既存 bvar が flow_float なので float 格納→int キャスト、または int 配列を別途用意 (実装時に既存規約へ合わせる)。
2. **`ransWallFunction_d.cu`**: `compute_wall_friction_sst_d` で node の `irep` を `irep_b[ib]` に保存、`wf_pk` 書き込みを node では `irep` 先に変更。
3. **`ransBoundary_d.cu`**: `rans_wall_scalar_boundary_d` に `irep_b` を渡し、node の ω/roOmega ピンを `irep` 先・y=`wall_dist[irep]` に変更。`compute_wall_y_eff` は不要化 (irep の wall_dist を直接使用) — 残すかは実装時判断。
4. **`update_d.cu`**: `applySSTPointImplicit_d` の ω decouple 判定を `wf_pk>=0` に変更 (`wf_pk` 引数追加、wrapper 配線)。
5. ビルド (build-native, full rebuild [[stale-build-struct-layout-trap]])。

## 6. 検証

- **ビルド**: `build-native` full rebuild。
- **回帰 (cell ビット不変)**: `case/18.backstep` cell (run_0057 相当) を本変更前後で比較し **ビット不変** (cell path 不変の担保)。`case/08`/`case/24` の cell SST も差分なし。
- **主検証 (node の効き)**: `case/18.backstep` node autowall (run_0059 設定 = nodeWallDirichlet+autowall+MUSCL+nStepInner5) を再走:
  - 段差凹コーナー (x≈2, y≈0.06-0.25) の **vis_t/l が 6789 → cell 同等 (~77 オーダー) に低下**。
  - 第一内層ノードが `wf_pk>=0` (壁関数処理されている) ことを `wf_pk` 出力で確認。
  - 近壁 k が cell と同オーダー (k≈3-5、現状 52 から低下)。
- **クロスチェック維持**: node の x_R が cell (7.95)・SU2 (7.89) と整合方向 (現状 node 6.71 が cell 寄りに動くか、少なくとも乱流場が cell と整合)。`check_convergence.py` VERDICT と `check_quasisteady.py` を併記。
- **判定基準**: ① cell ビット不変、② node コーナー vis_t/l が cell オーダー、③ node 近壁 k 暴走解消、④ node 残差が発散しない (定常到達)。

## 7. 影響範囲

- 触るファイル: `cuda_forge/ransWallFunction_d.cu`, `cuda_forge/ransBoundary_d.cu`, `cuda_forge/update_d.cu` (+ wrapper)、`boundaryCond.hpp`/`mesh.hpp` (irep_b bvar)。
- cell・wall-resolved (wallTreatmentSST=0)・非 node・explicit は無影響 (分岐ガード)。
- docs: `methods/turbulence/theory.md` §6.5 と `implementation.md` §3.7 に node の適用先 (第一内層ノード) を追記。`methods/index.md` 整合確認。

## 8. 完了条件

- [ ] `methods/turbulence/theory.md` §6.5 / `implementation.md` §3.7 に node 第一内層適用を追記
- [ ] 実装・検証完了 (§6 の判定基準を満たす)
- [ ] cell ビット不変を確認
- [ ] `status` を `done` に更新し §9 に変更ログ (before/after の vis_t/l・k を数値で)
- [ ] `plans/active/` → `plans/accepted/` へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期

## 9. 変更ログ

- `2026-06-28` — 初稿。`case/18.backstep` node 過剰乱流の真因 (壁関数処理が壁ノードのみで第一内層を取りこぼす) を k 収支実測で確定 ([notes/investigations/turbulence-node-wall-function-coverage.md](../../notes/investigations/turbulence-node-wall-function-coverage.md)) を受けて起案。
- `2026-06-28` — **実装・検証完了 (done)**。当初案 (wf_pk と ω ピンの両方を第一内層 irep へ移設) は **2 つの理由で失敗**: ① 凹コーナーで底壁/段差面が同じ irep を共有し ω ピン値が race、② 壁ノードの ω アンカが外れ ω 崩壊 → k がさらに暴走 (52→246-350、run_0062/0064 で実証)。**採用した正しい設計** = 「**生産 P_k 置換 (`wf_pk`) のみ壁ノード+第一内層ノードの両方に書き、ω ピン/残差ゼロ化/decouple は壁ノードのまま据え置く**」。
  - 実装: `ransWallFunction_d.cu` で node のとき `wf_pk[irep]` も書く (両ノード被覆)。`ransSource_d.cu` の ω 残差ゼロ化は node では `wall_flag==1` 限定に変更 (wf_pk は第一内層にも付くため; `wall_flag`/`isNode` 引数追加)。`rans_wall_scalar_boundary_d` / `applySSTPointImplicit_d` (ω ピン・decouple) は **baseline のまま不変**。`build-native` full rebuild。
  - **効果 (`case/18.backstep` node autowall, run_0065/0067)**: 段差コーナー vis_t/l **6789→77**・k **52→3** (cell 同等)。場平均 vis_t/l **424→207** (cell 198)。x_R **6.71→7.63** (cell 7.95・SU2 7.89 に接近)。machmax/pmax `STEADY`。
  - **cell ビット不変** (run_0066): vis_t/l mean 197.71・max 994 が baseline と一致 (cell 経路は `wall_flag=nullptr`/`isNode=0` で `wf_pk>=0` 判定=従来同一)。
  - **残課題 (別途)**: node 壁ノードは cell 第一セルより壁から遠く ω ピンが ~1/4 低いため、**再付着せん断層の近壁で μ_t ピークが残る** (vis_t/l ~5800 vs cell 990、局所 ~27 node)。場平均・x_R は整合。診断出力 `wf_pk`/`Pk_diag`/`roK_wf` は今後の調査用に残置。
- `2026-06-28` — **追加オプション: node k Dirichlet (`nodeKwfDirichlet`, 既定 0)**。上記再付着 μ_t ピーク対策として SU2 `SetTurbVars_WF` 流の k ハード Dirichlet を実装 (第一内層ノードで `roK=ρ·k_wf`, `k_wf=ω_w·μ_t,wall/ρ=ω_w·ν·(1/g-1)`; `applySSTPointImplicit` で dk=0, `ransSource` で res_roK 0 化)。新フィールド `roK_wf` (-1=inactive)、新フラグ `turbulence.nodeKwfDirichlet`。
  - **効果 (run_0068/0069_on)**: 再付着 μ_t ピーク除去 (vis_t/l max 5769→893 = cell 994 並、再付着 vt/l 5834→5、場平均 207→180)。
  - **副作用**: 非平衡再付着で u_τ→0 → k_wf→0 と再付着乱流を抑えすぎ **x_R 7.63→8.67 と伸びる** (cell 7.95/SU2 7.89 より長い)。「場の清浄さ」と「x_R 精度」のトレードオフ。
  - **判断**: **既定 OFF** (既定=k を解く prod-fix で x_R が cell/SU2 整合)。再付着 μ_t ピークを嫌う用途のみ opt-in。`nodeKwfDirichlet=0` で **prod-fix と一致** (run_0069_off が run_0067 と mean 206.8/max 5769/x_R 7.63 完全一致, roK_wf-active=0 で確認)。ω ピンは Dirichlet 版でも壁ノード据え置き (第一内層へ移すと race で崩壊するため)。

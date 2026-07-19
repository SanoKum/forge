# node-centered 壁 no-slip の implicit Jacobian Dirichlet 化 (SU2 DeleteValsRowi 相当)

## メタ

- **area**: `boundary / time_integration`
- **status**: `done`
- **related_docs**:
  - `methods/discretization.md` (§6.3 壁 Dirichlet)
  - `methods/discretization.md` (§ 軸 MARKER_SYM の in-Jacobian decouple 先例)
- **related_plans**:
  - 親: [`discretization-node-boundary-ghostless.md`](../active/discretization-node-boundary-ghostless.md) (Phase 1 で残差射影 Dirichlet を実装。本 plan はその implicit 整合化)
  - 先例: [`architecture-axisym-axis-singularity.md`](architecture-axisym-axis-singularity.md) (軸 roUy 行の in-Jacobian decouple)
  - 関連: [`gpu-implicit-plan.md`](../active/gpu-implicit-plan.md) (block-DPLUR 本体)
- **created**: `2026-06-22`
- **owner**: (未割当)

## 1. 目的

node-centered の壁 no-slip Dirichlet (`nodeWallDirichlet=1`) は現在 **残差射影 (`zeroWallMomentumResidual`) + 状態射影 (`enforceWallNoSlip`)** だけで実装され、**block-DPLUR の Jacobian は壁ノードの運動量を近傍と連成したまま**である。このため implicit では線形 solve が残差 0 でも `Δ(ρu)≠0` を返し、壁速度が drift する (実測: case/36 で壁 `|roUx|` mean 5.6→34 に増大、毎ステージ KE mean ~4000 を剥ぎ取る持続的エネルギーシンク)。本 plan は **block-DPLUR の壁ノード 5×5 Jacobian で運動量3行を decouple** (SU2 `DeleteValsRowi` 相当) し、`Δ(ρu)=0` を solve から一貫して得る。完了時、壁速度 drift・spurious KE 除去・近壁の蹴りが消え、壁 no-slip が implicit でも厳密に保たれる。

## 2. スコープ

- **やる**:
  - `implicit_defect_correction_block_d` の既存 in-Jacobian decouple 機構を壁に拡張: `wall_flag[ic]==1` で運動量3行 (index 1=roUx, 2=roUy, 3=roUz) を単位行化 + `rhs=0`。
  - 起動側で `nodeWallDirichlet=1 && discretization=="node"` のとき `msh.wall_flag_d` を渡す (現状 `nullptr`)。
  - 既存の状態射影 (`enforceWallNoSlip`) / 残差射影 (`zeroWallMomentumResidual`) は残す (explicit 経路・IC・保険)。Jacobian decouple 後は drift が消えるので両者は実質 no-op になる想定。
- **やらない**:
  - explicit 経路の変更 (explicit は残差射影だけで `Δ(ρu)=0` が成立済み)。
  - 壁圧力寄与を運動量バランスへ正しく入れる弱形式 half-face flux (= node-boundary-ghostless Phase 2)。本 plan は「速度 Dirichlet の implicit 整合」のみ。
  - 近壁勾配の退化 (dcc/dn) 対策 (別 plan [`diffusion-node-wall-viscous-distance.md`](diffusion-node-wall-viscous-distance.md))。
  - 軸 MARKER_SYM の corner 問題 (別系統、本 plan の壁とは別)。

## 3. 関連 docs と前提

- **理論的位置づけ**: 壁 no-slip は「壁ノードの運動量方程式を `u=0` (Dirichlet) に置換、連続・エネルギーは保存式のまま」。SU2 (vertex-centered median-dual) も同じで、implicit では `Jacobian.DeleteValsRowi(momentum rows)` により運動量行を単位行化する。forge は残差のみ 0 にして Jacobian を放置しているのが差分。
- **先例 (重要)**: 軸対称 MARKER_SYM で「残差射影だけでは implicit/block-DPLUR に非整合 (Mach~1000 発散)、in-Jacobian で roUy 行 (index2) を decouple すれば整合」が確立済み ([methods/discretization.md](../../methods/discretization.md#実装) の SU2 比較節)。本 plan は同じ機構を**壁の運動量3行**へ適用する。
- **既存インフラ**: `implicit_defect_correction_block_d` は既に `axis_flag` 引数と decouple コード ([timeIntegration_d.cu:837](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L837)) を持ち、`wall_flag_d` も node モードで構築済み ([mesh.cpp:803](../../solver_density_cuda/mesh/mesh.cpp#L803))。新規データ構造は不要。
- 理論変更は軽微 (既存 SU2 比較節に壁を追記) のため、実装前に methods/discretization/{theory,implementation}.md に壁 in-Jacobian decouple の小節を追記する。

## 4. 設計方針

block-DPLUR は点 implicit (対角 5×5 ブロック + 近傍は defect-correction の `dq_old` スイープ)。壁ノードの運動量行を単位行化すると、その点ブロック solve は `Δ(ρu)=Δ(ρv)=Δ(ρw)=0` を返し、近傍方程式は壁の運動量変化 0 を off-diagonal で受ける (= 壁運動量が固定 DOF として扱われる)。完全 sparse Jacobian を組まない点が SU2 と異なるが、DPLUR スイープ収束で SU2 `DeleteValsRowi` と同等の効果になる (軸 roUy decouple と同じ近似レベル)。

既存 axis decouple (1 行) と同じ形を運動量3行へ:

```cpp
// 既存 (axis, roUy のみ):
if (axis_flag != nullptr && axis_flag[ic] == 1) {
    for (int jj = 0; jj < 5; ++jj) diag_block[2][jj] = 0;
    diag_block[2][2] = 1; rhs[2] = 0;
}
// 追加 (wall, 運動量3行 = SU2 DeleteValsRowi 相当):
if (wall_flag != nullptr && wall_flag[ic] == 1) {
    for (int row = 1; row <= 3; ++row) {
        for (int jj = 0; jj < 5; ++jj) diag_block[row][jj] = 0;
        diag_block[row][row] = 1; rhs[row] = 0;
    }
}
```

連続 (行0)・エネルギー (行4) はそのまま → 壁 ρ, ρe_total は保存式で発展。u=0 なので `roe = ρe_int`(運動学、ガスによらず成立)で、圧力は solver の EOS が復元する: **CPG (`thermalMethod=0`) は `P=(γ-1)roe`、TP (`thermalMethod=2`) は NASA-9 で `e_int→T` 反転 → `P=ρR_mix·T`**。decouple は運動量3行だけを触り **thermo 行 (エネルギー) を変えないので CPG/TP 共通に動く (thermo-agnostic)**。コーナー (壁∩出口) は `wall_flag` が壁優先所有していれば運動量 decouple される (所有規約は親 plan Phase 1)。

## 5. 実装ステップ

1. **カーネル拡張** ([cuda_forge/timeIntegration_d.cu](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)): `implicit_defect_correction_block_d` テンプレートに `geom_int* wall_flag` 引数を追加 (axis_flag の隣)、§4 の運動量3行 decouple を axis decouple の直後に挿入。`lowMachPrecond=2` 別カーネル (`implicit_defect_correction_block_lowmach_d` 系) にも同処理を反映 (使うなら)。
2. **起動側配線** ([timeIntegration_d.cu:1221](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L1221) 付近の `FORGE_BDPLUR_ARGS`): 末尾 `nullptr /* axis_flag */` の隣に `wall_flag` 実引数を追加。`(cfg.discretization=="node" && cfg.nodeWallDirichlet==1) ? msh.wall_flag_d : nullptr` を渡す。`implicit_defect_correction_d` (非ブロック経路) は対象外 (block-DPLUR 限定)。
3. **既存フックの整理**: `enforceWallNoSlip` / `zeroWallMomentumResidual` はそのまま残置 (explicit・IC・保険)。decouple 後に壁 `|roUx|`≈0・KE 除去≈0 になることを §6 で確認し、確認できれば docs にその旨記載 (将来の簡素化候補)。

## 6. 検証

- **ビルド**: native (`build-native`)。`solverConfig.hpp` struct 不変 (新フラグ無し) だが、安全のため full rebuild ([[stale-build-struct-layout-trap]])。
- **回帰 (ビット不変)**: `discretization=cell` および `nodeWallDirichlet=0` の node は `wall_flag=nullptr` 経路で**従来とビット不変**であること (case/08, case/24 で確認)。
- **主検証 (本 plan の効き)**: `case/36.passive_pseudoshock_control` を node + implicit で再走 (現 `run_node_solid_dir_lam` の設定):
  - 壁ノード `|roUx|` が step 進行で**増大しない** (現状 mean 5.6→34 → 目標 ~0 を維持)。`enforceWallNoSlip` が剥ぐ KE が ~0 に落ちる。
  - 近壁 dP/dy の振動 (現状の症状) が低減するか定量比較 (GG 隣接偏差指標、`node_lam_GG_vs_LSQ_profile.png` と同手法)。
  - 全エネルギーの壁シンクが消える (保存性改善)。
- **回帰 (収束維持)**: `case/29.bell_vs_conical` の node 層流 (run_0032/0033 系) と node SST (run_0038 系) が**従来同等以上に収束** (デグレ無し)。`check_convergence.py` VERDICT を併記。
- **SST 再挑戦 (副次)**: case/36 node SST が omega 爆発のまま改善しないか確認 (omega 近壁は別要因 A だが、壁速度 drift が omega 生成 ∇k·∇ω を汚していた可能性の切り分けになる)。
- **判定基準**: ① cell/`nodeWallDirichlet=0` ビット不変、② 壁 `|roUx|` 非増大、③ case/29 デグレ無し (VERDICT 同等)、④ 近壁 dP/dy 振動が悪化しない。

## 7. 影響範囲

- 触るファイル: `cuda_forge/timeIntegration_d.cu` (カーネル引数 + decouple 数行 + 起動配線)。`mesh`/`solverConfig` は不変 (`wall_flag_d` 既存)。
- 既存ケース: cell・`nodeWallDirichlet=0`・explicit は無影響 (`wall_flag=nullptr` か非 block 経路)。
- ドキュメント: `methods/discretization.md` §6.3 (壁 Dirichlet の implicit 整合)、`methods/discretization.md` (SU2 比較節に壁 in-Jacobian decouple を追記)、`methods/index.md` の整合確認。

## 8. 完了条件

- [ ] `methods/discretization.md` に壁 in-Jacobian Dirichlet を追記
- [ ] `methods/discretization.md` の SU2 比較節に壁を追記
- [ ] 実装・検証完了 (§6 の判定基準を満たす)
- [ ] [`plans/README.md`](../README.md) の状態を `done` に更新
- [ ] 本 plan の `status` を `done` に変更し、§9 に変更ログ (壁 drift・KE シンク・近壁 dP/dy の before/after を数値で) を記載

## 9. 変更ログ

- `2026-06-22` — 初稿。症状切り分け (case/36 node 層流 implicit で壁 `|roUx|` が drift し KE 除去がエネルギーシンク化、原因 = block-DPLUR Jacobian が壁運動量を decouple していない) を受けて起案。打ち手 = 既存 axis_flag in-Jacobian decouple 機構を壁の運動量3行へ拡張 (SU2 `DeleteValsRowi` 相当)。
- `2026-06-22` — **実装・検証完了 (done)**。`implicit_defect_correction_block_d` に `wall_flag` 引数 + 運動量3行 (index1-3) decouple を追加 ([timeIntegration_d.cu](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu))、起動側で `(node && nodeWallDirichlet) ? msh.wall_flag_d : nullptr` を配線。`build-native` ビルド成功。
  - **効果 (case/36 node 層流 implicit, 3000 step)**: 壁ノード `|roUx|` mean **5.6→34 (旧) → 厳密 0 (新)**、`enforceWallNoSlip` の KE 除去 mean **~4000/stage → 0** (エネルギーシンク消滅)。近壁 P 隣接偏差 max **5.81e4→3.73e4 (−36%)**、mean 1141→1068。残差は旧と同等 (rms_ro 0.0188 vs 0.0203、plateau は pseudoshock 非定常で旧同様=本 BC 起因でない)。
  - **回帰**: cell implicit は baseline (git stash で本変更のみ除いた再ビルド) と比較し、差は **atomicAdd 自己揺れ (baseline 自身 run 間 3.2%) と同オーダー (3.2%)** = 実質無影響 (`wall_flag=nullptr` ガード)。node+Dirichlet は差が出る (decouple 有効)。case/29 の node 主レシピは explicit で block-DPLUR を呼ばないため構造上無影響。
  - **残課題 (別 plan)**: plateau (pseudoshock 非定常 → dual-time)、SST omega 近壁ソース (打ち手 A)、近壁勾配退化 ([diffusion-node-wall-viscous-distance])、弱形式 half-face flux (ghostless Phase 2)。本 plan 後は `enforceWallNoSlip`/`zeroWallMomentumResidual` は実質 no-op だが explicit・IC 経路の保険として残置。

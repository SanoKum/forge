# node 入口 Dirichlet の保存量整合と残差除外 (k/omega)

## メタ

- **area**: `turbulence` (boundary 横断)
- **status**: `done`
- **related_docs**:
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) (§3.5 境界条件 / §3.7 wall treatment)
- **related_plans**:
  - [`turbulence-node-wall-function-coverage.md`](../accepted/turbulence-node-wall-function-coverage.md) (壁ピン+res除外の被覆。本 plan は入口側の同種欠落)
  - [`turbulence-sst-freestream-init.md`](../accepted/turbulence-sst-freestream-init.md)
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 1. 目的

node-centered で `case/18.backstep` の **rms_roOmega プラトー (~9.5)** の主因が、入口×壁コーナーノードの
omega 残差スパイク (|res_roOmega|=568 = 入口内部 111 の5倍) であると局在化済み
([notes/investigations/backstep-lowmach-checkerboard.md](../../notes/investigations/backstep-lowmach-checkerboard.md))。
真因は node 入口 Dirichlet (`rans_dirichlet_scalar_boundary_d`) が **primitive `k`/`omega` のみピンし、
保存量 `roK`/`roOmega` を整合させず・残差も除外しない**こと (壁 BC は両方やる)。結果、入口ノードで
`roOmega/ro` vs `omega` が最大 888727、`roK/ro` vs `k` が 7729 不整合し、その残差が rms を汚染する。
本 plan は入口 Dirichlet を壁ピンと同形にして不整合と残差汚染を断つ。

## 2. スコープ

- **やる**: node 入口 (および任意の Dirichlet スカラー境界) で保存量 `roK[ic]=ρ·kb`/`roOmega[ic]=ρ·omegab` を
  整合ピン + `res_roK`/`res_roOmega`/`src_jac_*` を 0 化 (ransSource)。新フラグ `scalarDirichletPin` で識別。
- **やらない**: cell モード (ghost 経由で正しく課されるため不変)。omega/k の Dirichlet 値自体の変更
  (freestream `kb`/`omegab` は据え置き)。コーナーで壁 vs 入口の優先順位の再設計 (現状の last-writer のまま)。
  k の物理残差 (再付着生産) は本 plan の対象外 (入口不整合のみ是正)。

## 3. 関連 docs と前提

壁 scalar BC は `rans_wall_scalar_boundary_d` が `omega[ic]`/`roOmega[ic]` をピンし、ransSource が
`wall_flag[ic]==1` で `res_roOmega`/`src_jac_omega` を 0 化する
([implementation.md §3.7](../../methods/turbulence/implementation.md))。入口は `rans_dirichlet_scalar_boundary_d`
が node で primitive のみピン (保存量・res 未処理)。本 plan で入口側に同機構を追加する。

## 4. 設計方針

新フラグ `scalarDirichletPin` (cell 配列, 既定 0, node のみ使用)。
- `rans_dirichlet_scalar_boundary_d` (node 分岐): primitive ピンに加え
  `roK[ic]=ro[ic]·kb`・`roOmega[ic]=ro[ic]·omegab` を整合し、`scalarDirichletPin[ic]=1`。
- `rans_sst_source_d`: `scalarDirichletPin[ic]==1` (node) で `res_roK`/`res_roOmega`/`src_jac_k`/`src_jac_omega` を 0 化。
  壁の `omega_pinned` 除外と OR で合流 (壁ノードと入口ノードはどちらも除外対象)。
- cell モードは node 分岐を通らずフラグ 0 のまま・wrapper で nullptr 渡し → **ビット不変**。

## 5. 実装ステップ

1. `variables.hpp`: `cellValNames` に `scalarDirichletPin` を追加 (既定 0 確保)。
2. `cuda_forge/ransBoundary_d.cu`: `rans_dirichlet_scalar_boundary_d` に `ro`/`roK`/`roOmega`/`scalarDirichletPin`
   引数追加、node 分岐で保存量ピン+フラグ立て。`ransBoundary_d_wrapper` の呼び出しを更新。
3. `cuda_forge/ransSource_d.cu`: `rans_sst_source_d` に `scalarDirichletPin` 引数追加、node で res/src_jac 0 化。
   `ransSource_d_wrapper` で node 時のみ配列を渡す (cell は nullptr)。

## 6. 検証

- **ビルド**: native `build/` フルリビルド。
- **検証ケース**: `case/18.backstep` (node `run_0067` 設定, planar mesh)。
  - 収束場から restart し、入口ノードの `roOmega/ro`-`omega` 不整合 (888727) と `roK/ro`-`k` (7729) が解消するか。
  - **node rms_roOmega プラトー (~9.5) が低下するか** (`check_convergence` VERDICT)。
  - **x_R・場 (vis_t/l・再循環泡) が不変** (派生量を壊さない)。
- **cell ビット不変**: `run_0057` 設定で 2000 step、vis_t/l mean/max が baseline 一致。
- 併せて `verify-node-and-cell-both` に従い両 discretization で確認。

## 7. 影響範囲

- ファイル: `variables.hpp`, `cuda_forge/ransBoundary_d.cu`, `cuda_forge/ransSource_d.cu`。
- 既存ケース: cell は不変。node SST は入口残差が除外され rms_roOmega が改善 (場は不変の想定)。
- docs: `methods/turbulence/implementation.md` §3.5/§3.7 に入口 Dirichlet の保存量整合+res除外を追記。

## 8. 完了条件

- [x] `methods/turbulence/implementation.md` 更新済み (§3.5 に node 入口 Dirichlet の保存量整合+残差除外を追記)
- [x] 実装・検証完了 (§6)
- [x] `status` を `done`、§9 に変更ログ
- [x] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) 同期

## 9. 変更ログ

- `2026-06-28` — 初稿 (計画)。backstep node の rms_roOmega プラトー局在 (入口×壁コーナー) を受けて起票。

- `2026-06-28` — **実装・検証完了 (done)**。`scalarDirichletPin` (cell 配列) を新設し、
  `rans_dirichlet_scalar_boundary_d` (node 分岐) で `roK[ic]=ρ·kb`/`roOmega[ic]=ρ·omegab` を整合ピン+フラグ立て、
  `rans_sst_source_d` で `scalarDirichletPin==1` のノードの `res_roK`/`res_roOmega`/`src_jac_*` を 0 化。
  ファイル: `variables.hpp`, `cuda_forge/ransBoundary_d.cu`, `cuda_forge/ransSource_d.cu`。
  - **検証 (node, `run_0091` 75k 収束)**: **rms_roOmega プラトー 9.5→5.1e-3 (−約3桁、falling)** = 入口×壁コーナー残差
    (568→0) を除去。入口 `roOmega/ro`-`omega` 不整合 888727→2.5e-5、`roK/ro`-`k` 7729→20。場は物理的・NaN無。
  - **cell ビット不変 (`run_0088` vs 旧 binary `run_0089/0090`)**: 差は atomicAdd 非決定性ノイズ床と同等
    (旧vs新 ≈ 旧vs旧)。cell 経路は無影響。
  - **解への影響 (重要・許容)**: 本修正は残差ノイズ除去だけでなく **node の x_R を 5.63→6.67 と変える** (vt/l 207→183)。
    旧実装が入口近傍乱流を保存量不整合で汚していたためで、修正後 x_R は実験 6.26・従来 node 値 6.71 に接近 =
    **正味改善**。「x_R 不変」を満たさないのは旧実装のバグ是正の帰結 (user 確認のうえ確定)。
  - 診断用に `res_roOmega`/`src_jac_omega`/`transport_diag_omega`/`res_roK` を output に追加 (一時診断、wf_pk 系と同列)。
    run: `case/18.backstep/run_0083–0091`。

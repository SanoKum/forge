# 化学種の陰解法緩和整合 (matched-relaxation scalar-DPLUR)

## メタ

- **area**: `time_integration` / `thermophysics`
- **status**: `done` (案B 実装・検証完了。**結論=案B は要因2 の打ち手として無効**。詳細 §9)
- **related_docs**:
  - `docs/thermophysics/theory.md` (§3.1 陰解法での緩和整合)
  - `docs/thermophysics/implementation.md` (§4b matched-relaxation scalar-DPLUR)
  - `docs/time_integration/`
- **related_plans**:
  - `thermophysics-multicomponent-tpgas.md` (§10 残課題「多成分 implicit 結合の恒久修正 (要因2)」の具体化)
  - `time_integration-implicit-stable-cfl.md` (k/ω の transport_diag 陰化と同型の手法)
- **created**: `2026-06-19`
- **owner**: `CFD Dev`

## 1. 目的

多成分 TP (`thermalMethod==2`) の陰解法 (block-DPLUR) で、流れ 5 変数 (5×5 block・`nStepInner` sweep) と
化学種 `ρY_s` (commit 後に 1 回だけ点陰的更新) の**擬似時間緩和ミスマッチ** (要因2) を、化学種を
**流れと同一の緩和** (同一凍結残差・同一 `dt_local`・同一 `implicitRelax`・同一 `nStepInner` sweep) で
前進させて構造的に解消する。狙いは安定 `cfl_pseudo` 上限を現状 2 (sensible-enthalpy datum 後) から
引き上げ、残差プラトーを破ること。

## 2. スコープ

- **やる**: 化学種 `ρY_s` のスカラ DPLUR (緩和整合)。流れ DPLUR と同じ M 行列構造 (対角=流出質量流束 +
  `transport_diag`、非対角=流入質量流束) を Jacobi sweep で `nStepInner` 回緩和。config トグル
  `speciesImplicitCoupling` (既定 0=従来 segregated 点陰的・ビット不変)。
- **やらない**: 完全結合 (5+N の単一 block・固有値ベース |A| 一般化) は本 plan では行わない (案A は緩和整合で
  上限が上がらない場合の後続)。化学種拡散の非対角結合 (点陰的のまま据え置き)。化学反応源項。

## 3. 関連 docs と前提

- 要因2 の診断: `case/16.nozzle_wys/README.md`「多成分 TP 発散の再検証」(run_0069〜0109) と
  `thermophysics-multicomponent-tpgas.md` §9 (2026-06-18/19 ログ)。
- 緩和整合の理論は docs を先に更新済 (theory §3.1 / impl §4b)。
- 流れ scalar-DPLUR (`implicit_defect_correction_d`)・点陰的化学種更新 (`scalarTimeIntegration_d` の
  `timeIntegration==11` 分岐)・k/ω の `transport_diag` 陰化 (`time_integration-implicit-stable-cfl.md`) を
  ミラー元とする。

## 4. 設計方針

詳細は [implementation.md §4b](../../docs/thermophysics/implementation.md)。要点:

- **緩和整合**: `ρY_s` は受動移流 (接触波速のみ)。流れ 5×5 block に同梱せず、スカラ DPLUR で同一緩和に
  そろえる。凍結 Jacobian の対角 `D=V/Δτ + transport_diag_Y{s}`、非対角 = 流入 `max(∓ṁ,0)/ρ_nbr`。
- **退化保証**: `nStepInner=1`・`implicitRelax=1`・非対角 0 で従来点陰的に一致。`speciesImplicitCoupling=0`
  または `nSpecies<2` で従来経路 (ビット不変)。
- **timing**: `res_roY{s}`/`transport_diag` は `assembleResidual` で凍結。`ρ_nbr` 正規化は `roN`(=ρⁿ) を使い
  `transport_diag` の ρ と整合 (flow commit 後の ρ ではなく)。

## 5. 実装ステップ

1. **config** (`input/solverConfig.{hpp,cpp}`): `int speciesImplicitCoupling=0` を追加し `time.deltaT` から読む。
2. **バッファ** (`variables.cpp registerSpecies`, `speciesCellVarNames`): 化学種ごとに `dq_roY{s}` /
   `dq_roY{s}_old` を登録。
3. **カーネル+ドライバ** (`cuda_forge/speciesTransport_d.{cu,cuh}`): `species_dplur_sweep_d` (流れ scalar-DPLUR を
   1 本にミラー)、`species_commit_correction_d` (ρY=ρY_N+dq, floor 0)、`speciesImplicitDPLURSolve_d_wrapper`
   (memset→nStepInner sweep+swap→commit)、判定 `speciesImplicitCoupled(cfg,var)`。
4. **配線** (`main.cpp implicitNonlinearUpdate`): `speciesImplicitCoupling==1` で
   `speciesUpdateOuter→speciesImplicitDPLURSolve→speciesRenormalize→speciesPrimitive`、0 で従来経路。
5. **full rebuild** (`solverConfig.hpp` 変更のため): `ninja -t clean && ninja` (native)。

## 6. 検証

- **ビルド**: native full rebuild 成功 (stale build trap 回避)。
- **検証ケース**: `case/16.nozzle_wys` hot N2+H2O 1% (全>200K)。IC `ic_hot_href.h5`、源場 `run_0081/res_400.h5`。
  A/B 先例 `run_0101〜0109`。`cfl_pseudo` sweep (1,2,4,5,8,10) × `speciesImplicitCoupling` on/off。
- **判定基準**:
  - 安定 `cfl_pseudo` 上限が緩和整合で 2 超へ上がるか (主目的)。
  - `tools/check_convergence.py` の VERDICT を必須 (残差プラトー突破=PASS への接近)。
  - NaN/Inf 無し、`ΣY=1`、`Y∈[0,1]`、`T>200K`。
- **回帰**: ① 単成分 TP・CPG ビット不変 (`speciesImplicitCoupling` 経路を通らない)、② `case/28.cutler` 安定、
  ③ `case/08.bump` / `case/26.flat_plate_sst` 不変。

## 7. 影響範囲

- `input/solverConfig.{hpp,cpp}`、`variables.cpp`、`cuda_forge/speciesTransport_d.{cu,cuh}`、`main.cpp`。
- 既存ケース: `speciesImplicitCoupling` 既定 0 で従来動作不変。多成分 TP 陰解法ケースのみ 1 で有効化。
- docs: `docs/thermophysics/{theory,implementation}.md` 更新済、`docs/index.md` は項目変更なし (既存ファイル内追記)。

## 8. 完了条件

- [x] `docs/thermophysics/theory.md` §3.1 更新
- [x] `docs/thermophysics/implementation.md` §4b 更新
- [x] 実装・検証完了 (§6) — **案B は無効と判明 (§9)**
- [x] `.github/plans/README.md` の状態を同期
- [x] 本 plan の `status` を更新し §9 に変更ログ

## 9. 変更ログ

- `2026-06-19` — 初稿。docs (theory §3.1 / impl §4b) 先行更新。案B (緩和整合 scalar-DPLUR) を実装着手。
- `2026-06-19` — **実装完了 + 検証 = 案B は要因2 の打ち手として無効 (切り分け完了)**。
  - **実装**: config `speciesImplicitCoupling` (既定 0=従来 segregated・ビット不変)、`variables.cpp` に化学種ごと
    `dq_roY{s}/_old`、`speciesTransport_d.cu` に `species_dplur_sweep_d`/`species_commit_correction_d`/
    `speciesImplicitDPLURSolve_d_wrapper`/`speciesImplicitCoupled`、`main.cpp implicitNonlinearUpdate` に分岐配線。
    native full rebuild 成功。
  - **検証** (`case/16.nozzle_wys` README「案B 検証 run_0110〜0120」): hot N2+H2O(1%, href ON)。
    1. **緩和整合は cfl 上限を一切動かさない**。matched(=1) と off(=0) の発散 step が全 cfl で完全一致
       (cfl≤2 完走 / cfl3 step20=step20 / cfl4 step7=step7 / cfl5,8,10 同一)。cfl1/2 のプラトーも matched≡off
       (`check_convergence.py`=NOT CONVERGED、rms_roOmega 上昇律速)。
    2. **cfl4 の壁は化学種起因でない**: 単成分 N2 (H2O 無) も cfl4 で step7 (run_0117_cflp4) = 2成分と同一。
       href ON 後の cfl 上限 (~2) は**流れ+SST block-DPLUR 律速**。
    3. **無効の機構**: 本ケースの組成は完全一様 (Y_N2=0.98905/Y_H2O=0.010950 が全セル同値, min=max) → `res_roY≈0`
       → 化学種補正が matched/segregated とも ~0。要因2 の H2O ペナルティは化学種**移流**の緩和ではなく
       **EOS エネルギー再構成** (Newton 反転の組成依存 h_mix) にあり、移流緩和をそろえても触れない。
    4. cfl3 で 2成分(step20)<単成分(step55)=残存する組成-エネルギー結合ペナルティはあるが緩和率統一では消えない
       (cross-Jacobian ∂roe/∂roY=案A が必要)。ただし cfl4 で単成分も壁=案A の上積みは流れ/SST 上限に頭打ち。
  - **回帰**: default 2sp 経路 (=0) が旧バイナリと ~atomicAdd 一致 (run_0116 vs run_0104, step100 同値)、単成分も同様
    (run_0117 vs run_0081)。新経路は `speciesImplicitCoupling==1 && nSpecies>=2` でのみ起動。
  - **処置**: `speciesImplicitCoupling` は既定 0 で温存 (本ケースで無益だが組成勾配の大きいケース向けに flag-gated 残置)。
    **後続の高優先 = 流れ/SST 陰解法の cfl 上限引き上げ** (cf. `time_integration-implicit-stable-cfl.md` が flat_plate で
    cfl120 達成。hot ノズルで効かない理由=rms_roOmega 律速 の調査)。要因2 の残ペナルティ恒久修正には EOS cross-Jacobian
    (案A) が必要だが流れ/SST 上限が先に律速するため優先度は低い。

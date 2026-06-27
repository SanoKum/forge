# scalar DPLUR の軸対称ソース Jacobian 整合 (TP 対応)

## メタ

- **area**: `time_integration`
- **status**: `done`
- **related_docs**:
  - `docs/time_integration/implementation.md`
  - `docs/axisymmetric/implementation.md`
- **related_plans**: `gpu-implicit-plan.md` (block 陰解法本体), `architecture-axisym-axis-singularity.md` (block 粘性対角の幾何是正), `thermophysics-multicomponent-tpgas.md` (陰解法 Jacobian の per-cell γ 化)
- **created**: `2026-06-14`
- **owner**: `CFD Dev`

## 1. 目的

scalar DPLUR (`timeIntegration==11`, `blockDPLUR==0`) の LHS が、軸対称フープ源
`res_roUy += (P − τ_θθ)·A_planar` ([axisymmetricSource_d.cu]) の Jacobian を持たず、block
(`implicit_defect_correction_block_d` の `diag_block[2][*]`) と非整合だった。block と同形の対角成分を
scalar の roUy 方程式対角に陰化して整合させる。**thermally perfect (TP) でも整合**するよう係数 γ は
per-cell `gamma_arr[ic]` (TP=γ_mix(T) / CPG=cfg.gamma) を使う (block と同方針)。

## 2. スコープ

- **やる**: scalar カーネル `implicit_defect_correction_d` に `gamma_arr`/`isAxisymmetric`/`A_planar` を渡し、
  roUy 方程式の対角に `A_pl·((γ−1)u_y + 2μ/(ρ r_eff))` の非負側を加える (式別対角)。
- **やらない**: ソース Jacobian の**非対角 (式間連成) 項**は scalar の対角近似では表現できないため対象外
  (block のみ)。scalar DPLUR を 2 次精度の実用陰解法に昇格させること (下記の通り別要因で不可)。

## 3. 関連 docs と前提

- residual 側の軸対称源: `docs/axisymmetric/implementation.md` / [axisymmetricSource_d.cu] L41。
- block の対応する陰化: [timeIntegration_d.cu] `implicit_defect_correction_block_d` L879-891。
- per-cell γ 方針: `thermophysics-multicomponent-tpgas.md` (TP 整合の陰解法 Jacobian)。

## 4. 設計方針

block の `diag_block[2][2] += A_pl·((γ−1)u_y + hoop)`, `hoop = 2μ_total/(ρ r_eff)`, `r_eff = V/A_pl` と
同形を scalar の roUy 対角に加える。scalar は全 5 式で単一対角 `diag` を共有していたので、roUy 用に
`diag_roUy = diag + max(src_diag, 0)` を別途持ち `inv_diag_roUy` で roUy のみ割る。**非負クランプ**は
scalar 対角の正値性 (対角優位) を保つため。defect-correction なので**定常不動点 (ΔQ→0 の解) は不変**。

## 5. 実装ステップ

1. `implicit_defect_correction_d` シグネチャに `gamma_arr` (implicit_relax の直後)・`isAxisymmetric`/`A_planar`
   (unsteady_diag の直前) を追加。
2. body で `isAxisymmetric==1` のとき `src_diag` を計算し `diag_roUy` を作って roUy のみ適用。
3. 呼び出し側 (`timeIntegration_d.cu` の dispatcher) に `var.c_d["gamma"]`・`cfg.isAxisymmetric`・
   `(isAxisymmetric? A_planar : volume)` を渡す。

## 6. 検証

`case/29.bell_vs_conical` conical メッシュ + 等エントロピー IC, SLAU 軸対称 laminar, build-viscfix。

- **before/after でほぼ不変**: cfl0.5 NaN@557→557, cfl1 @246→246, cfl2 @85→82
  (`run_disent_scalar_*` vs `run_disent_scalarfix_*`)。軌跡は摂動するが発散段は同一。
- 空間切り分け (`run_disent_scalarfix_cfl2_diag`, dense 出力): 発散は**軸ではなく出口面・外壁コーナー**
  (x≈0.085, y≈0.032≈R_e) で P→1Pa 過膨張→ρ<0→逆流発散。
- **1 次精度では scalar DPLUR が安定収束** (`run_disent_scalarfix_o1_cfl0p5` −2.4dec, `_o1_cfl2` −1.6dec)。

## 7. 結論 (重要)

本整合 (軸対称ソース Jacobian の対角追加) は **block との整合・TP 対応として正しい改善**だが、
**case/29 の scalar DPLUR 発散の主因ではない**。真因は **2 次 MUSCL の出口コーナーでのオーバーシュートを
scalar の弱い (式分離・空間 1 次凍結) LHS が減衰できない**こと (block の真の 5×5 Jacobian は減衰可、
陽 RK3 も安定)。よって **2 次精度の陰解法は block DPLUR を使う**。scalar DPLUR は **1 次** (起動・ロバスト用)
に限る、が実用指針。詳細な切り分けは `case/29.bell_vs_conical/README.md`。

## 変更ログ

- `2026-06-14`: scalar カーネルへ軸対称ソース Jacobian 対角 (per-cell γ, TP 対応) を追加。検証の結果、
  block 整合・他軸対称ケース向けの正しい改善だが case/29 の発散 (出口コーナー 2 次 MUSCL 起因) は救えないと確定。
  status `done` (整合改善として)。

# 陽解法 RK のスカラー (k/ω) 源項 point-implicit 化

## メタ

- **area**: `time_integration`
- **status**: `done`
- **related_docs**:
  - `docs/time_integration/theory.md`
  - `docs/time_integration/implementation.md`
  - `docs/turbulence/implementation.md`
- **related_plans**: `architecture-rans-sst.md`（親 SST 実装）, `gpu-implicit-plan.md`（block 陰解法側の point-implicit）
- **created**: `2026-06-07`
- **owner**: `CFD Dev`

## 1. 目的

陽解法 RK (`timeIntegration==1/3`) で RANS (SST) を回せるようにする。現状、SST 消散項の
point-implicit 処理は block 陰解法 (`timeIntegration==11`) パスにしか無く、陽解法 RK では
スカラー (k/ω) を純陽的に積分するため、壁近傍の大きな $\omega$ と $D_\omega=\beta\rho\omega^2$ の
stiff 性で CFL を下げても初回サブステップで発散する。これを解消し、平均流の時間積分法
（陽 RK / block 陰）を揃えた RANS 比較を可能にする。

## 2. スコープ

- **やる**:
  - `runge_kutta_exp_scalar_d`（`timeIntegration==1/3` のスカラー RK）に、消散ヤコビアン
    `src_jac_k`/`src_jac_omega` を使った残差増分の point-implicit 減衰を追加。
  - `ScalarTransportDesc` に `src_jac` を持たせ、k→`src_jac_k`、ω→`src_jac_omega` を割り当て。
- **やらない**:
  - `runge_kutta_exp_scalar_4th_d`（4 次 RK）の RANS 対応（平均流 4 次自体が RANS 未検証のため別途）。
  - 生産・移流・拡散項の陰化（lagged のまま。block 陰解法と同じ方針）。
  - 近傍結合（純 point-implicit、消散 stiff 性のみ陰化）。

## 3. 関連 docs と前提

- 消散ヤコビアンの定義と block 陰解法側の point-implicit は
  [`docs/time_integration/theory.md`](../../docs/time_integration/theory.md) §"スカラー (k/ω) の陰解法"、
  [`docs/turbulence/implementation.md`](../../docs/turbulence/implementation.md) §6 を参照。
- `src_jac_k = β*ω`、`src_jac_omega = 2βω` は [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu)
  `rans_sst_source_d` が毎 `assembleResidual` で出力（RANS-SST 有効時のみ）。

## 4. 設計方針

陽解法 RK のスカラー更新を

$$
\rho\phi \leftarrow c_N\rho\phi^N + c_M\rho\phi^M
  + \frac{1}{1 + c_{\text{res}}\Delta t_{\text{loc}}\,\partial D/\partial(\rho\phi)}
    \cdot c_{\text{res}}\frac{\Delta t_{\text{loc}}}{V}\text{res}_{\rho\phi}
$$

とする。減衰係数 $1 + c_{\text{res}}\Delta t_{\text{loc}}\,\text{src\_jac}\ge 1$ は block 陰解法対角
$D_\phi = V/\Delta\tau + V\,\text{src\_jac}$ に $\Delta\tau/V$ を掛けた形と整合（消散のみ陰化）。
`src_jac=0`（乱流なし／LES WALE）では `fac=1` で従来の純陽的更新に一致するため LES は無影響。

呼び出し順序: `advanceExplicitRK` → `assembleResidual`（`scalarTransport`+`ransSource` で res と src_jac 確定）
→ `ransTimeIntegration_d_wrapper`（ここで src_jac 利用可能）。

## 5. 実装ステップ

1. `ScalarTransportDesc` に `flow_float* src_jac;` を追加。
2. `buildScalarDescs` で k に `src_jac_k`、ω に `src_jac_omega` を割り当て。
3. `runge_kutta_exp_scalar_d` に `src_jac` 引数を追加し、残差増分を `fac` で割る。
4. `ransTimeIntegration_d_wrapper` の `timeIntegration==1/3` 呼び出しで `desc.src_jac` を渡す。
5. docs (theory/implementation) 更新、本 plan 追加。
6. `case/18.backstep` 2D メッシュで RK3 RANS を実行し、発散しないこと・残差が下がることを確認。

## 変更ログ

- 2026-06-07: 実装完了。`scalarTransport_d.cu` の `runge_kutta_exp_scalar_d` に point-implicit 減衰 + realizability
  floor を追加。スカラー (k/ω) 源項側の stiff 性は本変更で陰的に扱える。

- 2026-06-07 (追記・重要): `case/18.backstep` 2D で RK3 RANS を検証した結果、**本変更だけでは低マッハ RANS を
  陽解法 RK で回せない**ことが判明した。発散には 2 系統あり、切り分けると:
  1. **SST 源項の stiff 性** — 本 plan の point-implicit 減衰 + floor で対処済み（壁近傍 ω の即時発散は解消）。
  2. **低マッハ平均流の陽解法定常不安定（別問題・本 plan 外）** — `unsteady=0`（局所時間刻みで定常へ反復）では、
     M≈0.1 の音響 stiffness によりエネルギー方程式 (roe) が先行して指数発散する。convMethod 0/1/2・cfl 0.2/0.3/0.8・
     cold/restart のすべてで発散（roe が ~1.4 倍/step で増大、乱流は巻き添え）。block 陰解法は LHS の音響固有値
     による減衰で収束するが、陽解法には相当機構が無い。`lowMachPrecond` はフラックス散逸 (Phase 1) のみで
     陽解法定常を安定化できない（散逸を下げる方向）。過去の explicit SST 検証 (run_0087-0090) は超音速ノズルで
     低マッハ stiffness が無かったため成立していた。
  → 結論: 低マッハ RANS 定常は **block 陰解法が必須**（`run_0001` で収束確認済み）。陽解法 RK での低マッハ
     平均流安定化（時間微分の完全前処理）は本 plan の範囲外で、別 plan が必要。

- 2026-06-07: なお `unsteady=1`（時間精度マーチ）では事情が異なり、`case/18.backstep` の LES (WALE, RK3,
  dt=1e-4, 音速 CFL≈0.68<1) は `run_0003_slau_les_rk3` で 2000 step を発散なく完走（残差有界の非定常乱流状態）。
  本 plan の point-implicit 減衰は LES (src_jac=0) には無影響で、LES の安定性とは独立。

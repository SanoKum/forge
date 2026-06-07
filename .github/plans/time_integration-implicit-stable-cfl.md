# 陰解法 (block DPLUR) の安定 cfl_pseudo 引き上げ

## メタ

- **area**: `time_integration`
- **status**: `done`
- **related_docs**:
  - `docs/time_integration/theory.md`
  - `docs/time_integration/implementation.md`
  - `docs/turbulence/implementation.md`
- **related_plans**: `time_integration-explicit-pointimplicit-sst.md`（同じ src_jac を共有）, `gpu-implicit-plan.md`（block 陰解法本体）
- **created**: `2026-06-07`
- **owner**: `CFD Dev`

## 1. 目的

定常 RANS (SST) 陰解法 (`timeIntegration==11`, `blockDPLUR=1`) の安定 `cfl_pseudo` が、発達場 restart でも
MUSCL ≈5 に制限される。診断の結果、律速は壁近傍の k/ω **輸送項（移流+拡散）の陽的扱い**だった
（平均流は静止し k/ω のみ先に発散、`scalarDiffusion=0` で発散消失）。k/ω 輸送スペクトル半径を
point-implicit 対角に陰化して安定 `cfl_pseudo` を引き上げる（目標 10〜20、達成 **120**）。定常解
（壁法則・$C_f$）は defect-correction により不変に保つ。検証ケースは `case/26.flat_plate_sst`。

## 2. スコープ（最終・実装結果反映）

切り分けの結果、`case/26.flat_plate_sst` の安定 `cfl_pseudo` 律速は当初想定（SST 生産項 stiff /
平均流粘性ヤコビアン）ではなく **k/ω 輸送項（移流+拡散）の陽的（lagged）扱い**だった
（§9 変更ログの診断参照）。最終スコープは次の通り。

- **採用（本命）**: k/ω 輸送項を point-implicit 対角に陰化。`scalarTransport_d.cu` の advection/diffusion
  カーネルで輸送スペクトル半径 $\Lambda^{T}_\phi=\sum_f[\max(\pm\dot m,0)+\mu_{\text{face}}|\delta|/dcc]/\rho$
  [m³/s] を面ループ集計し（`transport_diag_k`/`transport_diag_omega`）、`applySSTPointImplicit_d`（block 経路）
  と `runge_kutta_exp_scalar_d`（陽 RK 経路）の対角 $D_\phi=V/\Delta\tau+V\cdot\text{src\_jac}+\text{transport\_diag}$
  に加える。defect-correction のため定常解は不変。
- **不採用（#1 生産項）**: 生産振幅 $P_k/(\rho k)$, $P_\omega/(\rho\omega)$ を正対角に足す案を実装・計測したが、
  安定 `cfl_pseudo` を一切上げず（120 で不変）、ω は $P_\omega=\alpha\rho S^2$ が $\omega$ に依らず勾配ゼロ。棄却。
- **不採用（#2 平均流粘性ヤコビアン）**: 発散時も平均流 5 式は残差フロアに静止しており律速ではないため、
  block 版 `implicit_defect_correction_block_d` の粘性近傍結合は本ケースに無効。未着手（別ケースで必要時）。
- **やらない**: 壁垂直方向ライン陰解法（非構造でのライン抽出・ブロック三重対角は新規実装コスト過大、別 plan）。
  k/ω の近傍 ΔQ 結合（純 point-implicit を維持。対角化のみで目標達成のため不要）。

## 3. 関連 docs と前提

- 輸送項 point-implicit 対角 $\Lambda^{T}_\phi$ の定義・定常解不変性・検証は
  [`docs/time_integration/theory.md`](../../docs/time_integration/theory.md) §"輸送項 (移流+拡散) の point-implicit 対角"。
  segregated point-implicit 本体は同 §"スカラー (k/ω) の陰解法"。実装は
  [`docs/time_integration/implementation.md`](../../docs/time_integration/implementation.md) §`applySSTPointImplicit_d`。
- `transport_diag_k`/`transport_diag_omega` は [`scalarTransport_d.cu`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu)
  の `scalar_advection_first_order_d`/`scalar_diffusion_first_order_d` が面ループで集計。消費は
  [`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) `applySSTPointImplicit_d`（block 経路）と
  同 `scalarTransport_d.cu` `runge_kutta_exp_scalar_d`（陽解法 RK 経路）。消散項 `src_jac_k`/`src_jac_omega` は
  従来通り [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) `rans_sst_source_d` が出力。

## 4. 設計方針（採用案: 輸送項 point-implicit 対角）

k/ω の 1 次風上移流・2 点拡散の流束を $\rho\phi$ で微分した対角寄与 $-\partial\text{res}/\partial(\rho\phi)\ge0$ を
面ループで集計し、point-implicit 対角に加える（$\dot m_f$=面質量流束、$\mu_{\text{face}}=\mu+\sigma_\phi\mu_t$、
$\delta_f$=拡散幾何係数）:

$$
\Lambda^{T}_\phi = \sum_f \frac{\max(\dot m_f,0)+\max(-\dot m_f,0)}{\rho}
  + \sum_f \frac{\mu_{\text{face}}}{\rho}\frac{|\delta_f|}{|\Delta\mathbf r_f|}\quad[\mathrm{m^3/s}],
\qquad
D_\phi = \frac{V}{\Delta\tau} + V\frac{\partial D}{\partial(\rho\phi)} + \Lambda^{T}_\phi
$$

単位は $[\mathrm{m^3/s}]$ で $V/\Delta\tau$ と同じため $D_\phi$ にそのまま加える（src_jac のように $V$ を掛けない）。
平均流 block DPLUR が $A^\pm$・粘性半径を対角に持つのと同じ思想を k/ω に適用したもの。近傍 $\Delta Q$ 結合は
持たない純 point-implicit（対角化のみで目標達成）。defect-correction のため収束解（res=0）は不変。

## 5. 実装ステップ（実施済み）

1. `variables.hpp` に `transport_diag_k`/`transport_diag_omega`（cell 値）を追加。
2. `scalarTransport_d.cu`: `scalar_advection_first_order_d`/`scalar_diffusion_first_order_d` に `ro`・`transport_diag`
   引数を追加し対角を atomicAdd 集計。`scalarTransport_d_wrapper` 冒頭で `transport_diag_*` をゼロ初期化。
   `runge_kutta_exp_scalar_d` の `fac` に `transport_diag/v` を加算（陽 RK 経路）。
3. `update_d.cu`/`.cuh`: `applySSTPointImplicit_d` に `transport_diag_k`/`omega` を渡し $D_k,D_\omega$ に加算。
4. docs（theory/implementation）更新、本 plan・README 追記。
3. Docker 再ビルド（baseline 確認 → #1 適用 → 再ビルド）。
4. **検証 #1**: `case/26.flat_plate_sst` の発達場（`run_0007/res_120000.h5`）から restart し、
   `cfl_pseudo` を 5→10→15→20 と上げて安定上限を測る（短 run で発散有無を判定）。
5. #1 で目標（10〜20）に届かなければ **#2** を実装し、`timeIntegration_d.cu` の block 粘性近傍結合を追加、
   theory/implementation に #2 節を追記、再ビルド・再スイープ。
6. 到達した安定 `cfl_pseudo` で十分 step 回し、`postprocess_wall_law.py` で u⁺-y⁺・Cf-Re_x が
   既存 baseline と一致することを確認。`residual_history.png` 生成。

## 6. 検証

- **ビルド**: Docker `forge-solver:cuda-dev` で `tools/build.sh` 成功。
- **検証ケース**: `case/26.flat_plate_sst`（本番）。補助 `case/20.naca_ml/001.test/run_slau`、`case/08.bump`。
- **判定基準**:
  - 安定 `cfl_pseudo` が baseline（≈5）から上昇（目標 10〜20、到達上限を報告）。
  - 同一反復数での残差降下が baseline より速い（または同等以上）。
  - 壁法則（u⁺=y⁺ 粘性低層・log 則 $(1/0.41)\ln y⁺+5.0$）と Cf-Re_x（Schlichting $0.0592Re_x^{-1/5}$ と数%）が
    改良後も保たれる（収束解不変の確認）。

## 7. 影響範囲

- 触るファイル: `solver_density_cuda/cuda_forge/ransSource_d.cu`（#1）、必要なら
  `solver_density_cuda/cuda_forge/timeIntegration_d.cu`（#2）。
- `src_jac` を共有するため block 陰解法・陽解法 RK の両 RANS 経路に影響（LES は `src_jac=0` で無影響）。
- docs: `docs/time_integration/theory.md`・`docs/time_integration/implementation.md`（#2 実装時に粘性節追記）。
  新規 docs ファイルは作らないため `docs/index.md` は不変。

## 8. 完了条件

- [x] 関連 `docs/time_integration/theory.md` 更新済み（輸送項 point-implicit $\Lambda^{T}_\phi$）
- [x] 関連 `docs/time_integration/implementation.md` 更新済み（`transport_diag`）
- [x] 実装・検証完了（§6 を満たす：安定 `cfl_pseudo` 5–6→120、壁法則不変）
- [x] `.github/plans/README.md` の状態を `done` に更新
- [x] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-06-07` — 初稿。当初は #1 生産項陰化を本命に docs（生産項減衰 $J^{P}_\phi$）を先行更新。

- `2026-06-08` — **診断と方針転換**。`case/26.flat_plate_sst` 発達場（`run_0007/res_120000`）restart で
  baseline 安定 `cfl_pseudo` を計測すると **5 安定 / 7–8 発散**（README の ≈5 と一致）。発散の onset を追うと、
  平均流 5 式（`rms_ro`〜1e-9, `rms_roe`〜4e-4）は残差フロアに**静止したまま k/ω だけが先に指数発散**。
  `scalarDiffusion=0`（k/ω 拡散を切る）と発散が消えることから、律速は **k/ω 輸送項の陽的扱い**と確定
  （壁第一セル高さ 4e-6 m の陽的拡散 $\Delta t<\delta^2/2\nu$ 限界・外層の陽的移流 CFL 限界）。
  → #1 生産項は ceiling を動かさず（実装・計測のうえ棄却）、#2 平均流粘性は律速でないため無効と判断。

- `2026-06-08` — **実装（採用案）**。k/ω 輸送項を point-implicit 対角に陰化（`transport_diag`）。
  `scalar_advection_first_order_d`（移流: $\max(\pm\dot m,0)/\rho$）と `scalar_diffusion_first_order_d`
  （拡散: $(\mu_{\text{face}}/\rho)|\delta|/dcc$）の面ループで集計し、`applySSTPointImplicit_d` と
  `runge_kutta_exp_scalar_d` の対角に加算。変数 `transport_diag_k`/`transport_diag_omega` を追加。

- `2026-06-08` — **検証結果（`case/26.flat_plate_sst`, MUSCL, 発達場 restart, 3000 step スイープ）**:
  | 構成 | 安定 `cfl_pseudo`（安定 / 発散） |
  |---|---|
  | baseline（消散のみ陰化） | 5 / 7 |
  | + 拡散のみ陰化 | 13 / 15 |
  | + 移流も陰化（最終） | **120 / 130** |

  → 安定 `cfl_pseudo` を **5–6 → 120（約 20–24 倍）** に向上（目標 10–20 を大幅超過）。
  発散の最終律速は ω 方程式（残る陽的項：交差拡散・近傍 ΔQ・生産）。
  - **定常解不変**: `cfl_pseudo=50` を発達場から 3000 step 回した壁法則は baseline（`cfl_pseudo=5`, 120k step）と
    $u_\tau$ 2.707→2.708、$C_f(x{=}0.3)$ 3.117e-3→3.116e-3（差 <0.1%）、Schlichting 比 同一。defect-correction の
    結果不変性を実測で確認。
  - **収束加速（cold start, 1次風上, 4000 step）**: baseline `cfl_pseudo=2`（README の spinup 設定）に対し
    implicit は `cfl_pseudo=50` でも安定。同 step で `rms_roUx` は impl(cfl20)=1.3e-5 < base(cfl5)=2.1e-5 < base(cfl2)=4.1e-5、
    $C_f(x{=}0.3)$ は impl(cfl50)=4.23e-3 が base(cfl5)=4.73e-3 より収束値 3.12e-3 に近い（同 step でより発達）。
  - LES（`LESorRANS=1`）は `scalarTransport`/`applySSTPointImplicit` を呼ばないため無影響。

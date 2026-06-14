# 時間積分 — 実装

forge の時間積分・更新カーネルの実装とソース対応をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/cuda_forge/timeIntegration_d.cuh`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cuh) | カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/timeIntegration_d.cu`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) | RK / 陰解法カーネル本体 (約 950 行) |
| [`solver_density_cuda/cuda_forge/update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) | 外側・内側ループ前後の Q コピー、陰解法補正反映 |
| [`solver_density_cuda/cuda_forge/implicitCorrection_d.cu`](../../solver_density_cuda/cuda_forge/implicitCorrection_d.cu) | dual-time 陽 (簡易) スキームの補正 |
| [`solver_density_cuda/cuda_forge/setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu) | 局所 $\Delta t_{\text{loc}}$ 計算 |
| [`solver_density_cuda/update.cpp`](../../solver_density_cuda/update.cpp) | 旧 CPU 経路 (現状未使用) |

## エントリポイント

[`timeIntegration_d_wrapper`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L780)
が共通入口。`cfg.timeIntegration` の値に応じて次を呼ぶ。

| 値 | 呼び出し |
| --- | --- |
| `4` | `runge_kutta_exp_4th_d` |
| `1`, `3` | `runge_kutta_exp_d` |
| `11`, `blockDPLUR == 1` | `implicit_defect_correction_block_d`（5×5 block DPLUR、LU-SGS $A^\pm$。**推奨既定**） |
| `11`, `blockDPLUR == 0` | `implicit_defect_correction_d`（scalar 対角版＝スペクトル半径。軽量だが擬似 CFL が低く収束が遅い） |

> `solverConfig::initTimeIntegrationScheme` の `case 11` は `blockDPLUR ∈ {0,1}` を受理（それ以外を throw）。両者とも古典 DPLUR 制御フロー（残差固定 + `nStepInner` sweep + 単一 commit）に対応。`unsteady == 1`（dual-time）は dispatcher 側で throw する（本体未実装）。

> **軸対称ソース項の Jacobian**: scalar 版 `implicit_defect_correction_d` も block と整合させ、軸対称フープ源 `res_roUy += (P − τ_θθ)·A_planar` の Jacobian 対角成分 `A_pl·((γ−1)u_y + 2μ/(ρ r_eff))`（非負側、per-cell γ で TP 整合）を roUy 方程式の対角に陰化する。ただし scalar は式間連成 (源 Jacobian の非対角) を表現できないため、軸近傍以外が律速の強膨張ケース（例 `case/29` 出口コーナーの 2 次 MUSCL オーバーシュート）は救えない。**2 次精度の陰解法は block DPLUR 推奨**、scalar DPLUR は 1 次（起動・ロバスト用）に限るのが実用指針（切り分けは [`time_integration-scalar-dplur-axisym-source.md`](../../.github/plans/time_integration-scalar-dplur-axisym-source.md) / `case/29.bell_vs_conical/README.md`）。

## ループ全体 ([`main.cpp`](../../solver_density_cuda/main.cpp))

`advanceOneStep` は巨大 lambda を廃し、`StepContext`（cfg, cuda_cfg, msh, mat_ns, var, fluct,
pprobes, profiler, residual_logger, implicit_diag_logger, iStep を束ねた参照集約構造体）を
受け取る自由関数群に分解する。スキームは 3 階層構造で、dual-time が後付けできる形にする。

```
advanceOneStep(ctx):                       // dispatcher
  if (isImplicit):
     if (unsteady) advanceImplicitDualTime(ctx)   // 後続フェーズ: 今は throw
     else          advanceImplicitSteady(ctx)
  else            advanceExplicitRK(ctx)

assembleResidual(ctx, stage):              // 残差組み立ての単一情報源（旧 assembleCurrentState）
  updateInner → dependentVariables → gasProperties → applyBconds → applyRansScalarBoundaries
  → calcGradient → axisymmetricGeomTerms → limiter → ducrosSensor → turbulent_viscosity
  → convectiveFlux → ransTransport → ransGradient + ransSource
  → axisymmetricSource → viscousFlux
  → [addUnsteadyTimeTerm(ctx)]            // dual-time の BDF 物理時間項フック（定常は no-op）

advanceExplicitRK(ctx):                    // tI 1/3/4（挙動不変）
  for iloop in perStepIterationCount():
     updateVariablesInner; assembleResidual(ctx,iloop+1); logResidualSnapshot
     timeIntegration_d_wrapper(iloop); ransTimeIntegration_d_wrapper(iloop)
  updateVariablesOuter; writeStepOutputs; setDT; logOuterEnd

implicitNonlinearUpdate(ctx):              // 定常・dual-time 共有の核
  assembleResidual(ctx, 1)                 // 残差・フラックスは 1 回（ransSource が src_jac_k/ω も出力）
  setDT_d_wrapper                          // 局所擬似時間 dτ（diag の V/Δτ）
  blockDPLURSolve(ctx)                     // 下記（平均流 5 式）
  applyBlockImplicitCorrection(ctx)        // Q = Q_baseline + dq を 1 回 commit
  if scalarResidualEnabled:                // RANS(SST) のとき
     applySSTPointImplicit(ctx)            // k/ω を segregated point-implicit で更新（凍結解除）

blockDPLURSolve(ctx):                      // 古典 DPLUR 線形ソルバ（res・Q 固定）
  for iSweep in nStepInner:
     implicit_defect_correction_block_d(...)   // 固定 res + lagged dq_old → dq_new
     swapBlockImplicitCorrectionBuffers(var)

advanceImplicitSteady(ctx):                // 定常: 擬似時間=メインループ、1 更新/step の縮退形
  updateVariablesInner; logOuterBegin
  implicitNonlinearUpdate(ctx); logResidualSnapshot
  updateVariablesOuter; writeStepOutputs; setDT; logOuterEnd
```

`updateVariablesOuter_d` は外側ループ開始時に `Q_N`, `Q_M` の両方を現在値に揃え、
`updateVariablesInner_d` は内側ステージ後の `Q_M` のみを更新する。

## RK カーネル詳細

### `runge_kutta_exp_d` (Jameson 多段)

ステージ係数 `coef_N, coef_M, coef_Res` を内部ループ index `loop` で参照し、

```cpp
Q[ic] = coef_N * Q_N[ic] + coef_M * Q_M[ic]
      + coef_Res * res[ic] * dt_local[ic] / vol[ic];
```

を全成分に適用。`dt_local` は `setDT_d_wrapper` で事前に書き込まれている。

### `runge_kutta_exp_4th_d`

低 storage 4 段 4 次 RK。`loop == 0` で残差累積バッファ `res_*_m` をゼロクリアし、
各段で `res_*_m += coef_Res * res * dt_local / vol`。
`loop < 3` の中間段では `Q = Q_N + coef_DT * res * dt_local / vol`、
最終段 (`loop == 3`) で `Q = Q_N + res_*_m` として確定。

### `runge_kutta_exp_scalar_d`（スカラー k/ω の陽解法 RK ＋ point-implicit 源項）

[`scalarTransport_d.cu`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu) の
`ransTimeIntegration_d_wrapper`（[`ransTransport_d.cu`](../../solver_density_cuda/cuda_forge/ransTransport_d.cu)）が平均流 RK と同じ段で k/ω を別カーネルで積分する
（`timeIntegration==1/3` は `runge_kutta_exp_scalar_d`、`==4` は `runge_kutta_exp_scalar_4th_d`）。

RANS (SST) の消散項・輸送項は stiff なため、`timeIntegration==1/3` の更新は**残差増分のみ源項+輸送ヤコビアンで減衰**する:

```cpp
const flow_float fac = 1.0 + coef_Res * dt_l * (src_jac[ic] + transport_diag[ic] / v); // ≥ 1
rho_phi[ic] = coef_N * rho_phi_N[ic] + coef_M * rho_phi_M[ic]
            + (coef_Res * res_rho_phi[ic] * dt_l / v) / fac;
```

- `src_jac`（消散 $\beta^\*\omega,\,2\beta\omega$）は `ScalarTransportDesc` 経由で k=`src_jac_k`、ω=`src_jac_omega`
  （[`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) が毎 `assembleResidual` で出力）。
- `transport_diag`（移流+拡散の対角 $\Lambda^{T}_\phi$ [m³/s]）は `scalar_advection`/`scalar_diffusion` カーネルが
  面ループで集計（k=`transport_diag_k`、ω=`transport_diag_omega`）。$V$ で割って $[1/s]$ 化して `fac` に入る。
- `LESorRANS!=2`（乱流なし／LES WALE）では両者 0 で `fac=1`、従来の純陽的更新と一致（LES は無影響）。
- 減衰係数は `applySSTPointImplicit_d` の対角 $D_\phi=V/\Delta\tau+V\cdot\text{src\_jac}+\text{transport\_diag}$ に $\Delta\tau/V$ を掛けた形と整合。
- `runge_kutta_exp_scalar_4th_d`（4 次）は本減衰未適用＝RANS 非対応のまま（平均流 4 次自体が RANS 未検証）。
- 理論は [theory.md](theory.md) §"陽解法 RK での point-implicit 源項"。

## 陰解法カーネル詳細

### `implicit_defect_correction_block_d`（block DPLUR、LU-SGS $A^\pm$ 分割）

5×5 ブロック版。`block_dplur` 名前空間に補助 device 関数群を持つ。カーネルは線形 solve の内部精度
`ST` (float/double) で `template<typename ST>` 化され、状態/残差 (float) を `ST` へキャストして取り込む
(混合精度。下記「閉形式 FVS と混合精度」参照)。

- `accumulate_split_jacobian_cf<T>`（**既定の閉形式 FVS**）— $R,\Lambda,L$ の 5×5×5 三重積を作らず、
  **音響右/左固有ベクトルのみ**で $M(g)=g_2 I+(g_1-g_2)\,r_1\!\otimes l_1+(g_5-g_2)\,r_5\!\otimes l_5$
  ($=R\,\mathrm{diag}(g)\,L$) を構成し、$a^{+}=M(\Lambda^{+})$ を対角へ・$k_{\rm off}=M((-\Lambda)^{+})$ を近傍へ
  **直接畳み込む**（$a^{+}/k_{\rm off}$ を materialize しない）。`build_jacobian_split` (下記) と軸対称 (nz=0) で
  数値等価かつ ~10% 高速 (レジスタ・スピル削減)。
- `build_jacobian_split` — （旧経路・precond 版が使用）固有分解 $R,\Lambda,L$ から **対角用 $A^{+}=R\,\Lambda^{+}L$** と
  **RHS 近傍用 $K=-A^{-}=R\,(-\Lambda^{-})L=\tfrac12(|\widetilde A|-\widetilde A)$** を同時に返す
  （$\Lambda^{+}=\max(\Lambda,0)$、$-\Lambda^{-}=\max(-\Lambda,0)$、共に非負）。
- `add_identity_scaled`, `add_scaled_5x5` — $V/\Delta\tau\,I$・粘性対角・面寄与の加算。
- `solve_5x5` — 部分ピボット付き Gauss 消去（`diag` を破壊して in-place）。`|pivot| < 1e-20` でゼロ解にフォールバック。
- `multiply_add_5x5_vec` — 行列ベクトル積。

各セルで近傍寄与を集約し、対角 $D_i = V/\Delta\tau\,I + \sum_f A^{+}_f S_f + \sum_f \Lambda^{\nu}_f\,I$、
RHS = $-\mathbf R + \sum_f K_f S_f \cdot \Delta\mathbf Q_{\text{nbr}}^{\text{old}}$（$K_f=-A^{-}_f$）を構築し、
$\Delta\mathbf Q_{\text{new}} = D_i^{-1}\,\text{RHS}$ を解く。`cfg.implicitRelax` で $\Delta\mathbf Q$ を緩和。

ここで粘性対角は $\Lambda^{\nu}_f = 2\nu_f\,\dfrac{|S_f|^2}{\Delta\mathbf{cc}_f\cdot S_f}$（$\nu_f=(\mu_{\rm lam}+\mu_t)/\rho$）。
これは粘性流束 residual ([`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu)) の
法線拡散項 $\mu_f(\Delta U/|\Delta\mathbf{cc}|)\,\delta$（$\delta=|\Delta\mathbf{cc}|\,|S_f|^2/(\Delta\mathbf{cc}_f\cdot S_f)$）の
**Jacobian の大きさと整合**する（$\Lambda^{\nu}_f = 2\nu_f\,\delta/|\Delta\mathbf{cc}|$）。

> **2026-06-14 修正（粘性対角の幾何是正）**: 旧コードは粘性対角を $|S_f|\cdot(2\nu_f/\delta)$ と書いていたが、
> $\delta$ は**面積次元** ($\approx|S_f|$) なので $|S_f|$ が約分されて $\approx 2\nu_f$ に潰れ、(1) 軸対称近軸で
> 本来 $\propto r$ で消えるべき内側面の寄与を過大評価し、(2) residual に無いゼロ面積(対称/軸)面にも
> スプリアス項を載せていた（residual 側は `ip<nNormalPlanes` で境界面を除外）。これが **float block-DPLUR が
> 軸対称近軸第一セルの $U_r$ を収束させきれず固着する真因**だった。上記の residual 整合形
> $2\nu_f\,|S_f|^2/(\Delta\mathbf{cc}_f\cdot S_f)$ に是正すると $\propto r$ でゼロ面積面では消え、**float のまま固着が解消**
> （case 29 laminar conical 第一セル $U_r$: $+1.4\to+17.9$ で double solve と一致、1 次では未収束→収束）。
> 修正は `timeIntegration_d.cu` の scalar (`implicit_defect_correction_d`) / block (`implicit_defect_correction_block_d`) /
> precond (`implicit_defect_correction_block_precond_d`) の 3 箇所。**LHS のみの変更**で defect-correction の
> 定常解は不変（planar 回帰 bump で base/fix 場が $L2\sim10^{-5}$ 一致・RANS で残差レベル同一を確認）。

> **2026-06 修正**: 旧コードは対角に $A^{+}$ ではなく $|\widetilde A|$ を、近傍に $-A^{-}$ ではなく $+|\widetilde A|$ を
> 使っていた（符号付き分割でなく絶対値の誤用）。対角が upwind 自己 Jacobian と不一致・近傍結合が逆符号となり、
> block DPLUR は収束せず発散していた。`build_jacobian_split` による $A^{+}/{-}A^{-}$ 分割でこれを修正。
> `res_*` の符号は $-\mathbf R$（陽解法 `runge_kutta_exp_d` の `Q=Q_N+\mathrm{res}\,\Delta t/V` と整合）なので、
> カーネルでは `rhs` を `res_*` で初期化し近傍寄与 $K_f S_f\,\Delta\mathbf Q_{\text{nbr}}$ を**加算**する。

**古典 DPLUR の構造（重要）**: カーネルは `dq_block_new` の生成のみを行い、
`Q`（ro..roe）を**インライン更新しない**。`blockDPLURSolve` が `nStepInner` 回 sweep を回し、
各 sweep 後にドライバ側で [`swapBlockImplicitCorrectionBuffers`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cuh)
（block 専用 swap、`.cuh` に公開）で `dq_block_old <-> dq_block_new` を入れ替える。
全 sweep 後に [`applyBlockImplicitCorrection`](../../solver_density_cuda/cuda_forge/update_d.cu)
（`applyScalarImplicitCorrection` と対称、`update_d.cu` に新設）で `Q = Q_baseline + dq_block` を
**1 度だけ** commit する。残差 `res_*` と $|\widetilde A_f|$ は sweep 中固定（matrix-free のため固定 Q から毎 sweep 再構築してよい）。

#### 閉形式 FVS と混合精度 (`implicitSolvePrecision`)

`accumulate_split_jacobian_cf<T>` は固有ベクトル行列 $R,L$ を陽に作らず、$\mathrm{diag}(g)-g_2 I$ が
shear/entropy の 3 モードを消すことを使って **音響右/左固有ベクトル $r_1,r_5,l_1,l_5$ のみ**で
$M(g)=R\,\mathrm{diag}(g_1,g_2,g_2,g_2,g_5)\,L = g_2 I+(g_1-g_2)\,r_1 l_1^\top+(g_5-g_2)\,r_5 l_5^\top$
を構成し、$a^{+}=M(\Lambda^{+})$ を対角へ・$k_{\rm off}=M((-\Lambda)^{+})$ を近傍へ直接畳み込む。
$R\,\Lambda\,L$ の 5×5×5 三重積と $a^{+}/k_{\rm off}/\text{solve\_mat}$ 保持を排除し、float 陰解法で ~10% 高速
(レジスタ・スピル削減)。**軸対称・平面 (nz=0) では legacy `build_jacobian_split` と数値厳密一致**
(一般 3D は legacy の $R,L$ が厳密逆行列でない=$RL\neq I$ ため僅差だが、固有値を厳密に $\max(\lambda,0)$ とする
valid な FVS で defect-correction の定常解は不変)。

カーネルは線形 solve の内部精度 `ST` で `template<typename ST>` 化。`solverConfig` の
**`implicitSolvePrecision`** (`time.deltaT`、既定 `0`) で切替える:

- `0` (float): 状態/残差 (float) のまま float で組立・solve（既定・高速）。
- `1` (double): 状態/残差を **double へキャスト**して Jacobian 構築・5×5 solve・近傍 sweep を **double** で行い、
  補正 $\Delta\mathbf Q$ を float `dq_new` へ書戻す（混合精度 iterative refinement）。

**動機と位置づけ（重要・2026-06-14 更新）**: 当初 float32 の block-DPLUR が軸対称 近軸第一セルで平均速度 `Uy` を
収束させきれず偽固着する (laminar conical で `Uy` が物理値 $-15$ でなく $-0.6$) 問題に対し、`implicitSolvePrecision=1`
(線形 solve を double) を root-fix とした。**その後、真因は精度ではなく上記「粘性対角の幾何不整合」であることが判明**
(粘性が無い Euler は float でも固着せず、固着は粘性 LHS 由来。詳細は本節冒頭「粘性対角の幾何是正」)。
粘性対角を residual 整合形に直せば **float のまま固着が解消**し double solve は不要。
したがって `implicitSolvePrecision=1` は**根治ではなく、悪条件 LHS を倍精度で押し切る検証/保険用の手段**として残す
(幾何是正後は通常 `0` で良い)。RTX 3060 では FP64=FP32 の 1/32 ゆえ ~×2.8 遅い。
切り分け・速度の詳細は [`.github/plans/precision-mixed-axisym.md`](../../.github/plans/precision-mixed-axisym.md)
と [`.github/plans/architecture-axisym-axis-singularity.md`](../../.github/plans/architecture-axisym-axis-singularity.md)。
現状 `blockDPLUR=1`・`lowMachPrecond` 0/1 経路のみ対応 (precond=2 / scalar 版は float のまま)。

**低マッハ前処理 (LHS 固有値) は不採用** (2026-06 検証・実装後 revert)。`build_jacobian_split` の
固有値 `lambda[5]={U+c,U,U,U,U-c}` を前処理固有値 `U±c'` に差し替える案を実装・検証したが、
block DPLUR では**対角優位性の源である大きい音響固有値を縮めてしまい有害**（フラックス散逸前処理
単独で安定だった `eps=0.15` すら発散させ、安定 `eps` 範囲を狭めた。収束加速も根治もなし）。よって
LHS は従来の $A^\pm$（物理音速 `sonic`）のまま。低マッハ前処理は `SLAU_d` の散逸スケールにのみ適用する。
根拠・データは [`theory.md`](theory.md#低マッハ前処理固有値-weisssmith--試行したが不採用) と計画
[`time_integration-lowmach-preconditioning.md`](../../.github/plans/time_integration-lowmach-preconditioning.md) §9。

### `implicit_defect_correction_block_precond_d`（Phase 4: 完全 $\Gamma^{-1}A$ 前処理・`lowMachPrecond=2`）

上の LHS 固有値前処理 (不採用) と異なり、**前処理を一貫させた別カーネル**。`blockDPLUR==1 && lowMachPrecond==2`
のとき wrapper がこちらを起動する (既存 `implicit_defect_correction_block_d` は 0/1 専用・ビット/レジスタ不変)。

- **Sherman-Morrison 解法 (FP64 回避・高速)**: $\Gamma_c=I+\alpha g r^\top$ がランク 1 ゆえ対角ブロックは
  $D=D_0+\gamma g r^\top$ ($D_0$=V/Δτ'·I+物理 FVS+粘性+軸対称=既存 block と同形・良条件、$\gamma=(V/\Delta\tau')\alpha$)。
  $D_0$ を **float** で 2 RHS 同時 ([`solve_5x5_2rhs`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)) に解き
  $y=D_0^{-1}b,\ z=D_0^{-1}g$、$x=y-[\gamma(r^\top y)/(1+\gamma(r^\top z))]z$。悪条件 $\sim1/\beta$ は分母スカラーのみ double。
  consumer GPU (FP64=FP32/64) でも物理 block とほぼ同速 (per-step 16.9 vs 17.8ms)。$g,r,\alpha$ はカーネル内インライン。
- **時間項**: $\Gamma_c\,V/\Delta\tau'$ (上の $\gamma g r^\top$ 寄与)。$\Delta\tau'$ は `setDT_d` の
  `setDTlocal_precond_scale_d` が `dt_local *= (|u|+c)/ρ'` で拡大済 ($\rho'$=前処理スペクトル半径)。
- **フラックス分割**: **物理の厳密 FVS** をそのまま使う (対角に $a^{+}=A_c^{+}$、近傍に $k_{\rm off}=-A_c^{-}$。既存 block と同一)。
  保存形 $(\Gamma_c V/\Delta\tau' + A_c)\Delta Q=-R$ より前処理は**時間項のみ**で、$A_c$ は残差の真のヤコビアンゆえ非前処理が正しい
  ([`theory.md`](theory.md) 参照)。収束は $(\Delta\tau'/V)\Gamma_c^{-1}A_c$ の固有値 $\lambda'$ で一様に前処理される。
  ※当初フラックスも $\hat A_\Gamma=\Gamma_c^{-1}A_c$ のスペクトル半径分割で前処理したが、別系で過散逸ゆえ撤回 (上記が正)。
- **その他**: 粘性スペクトル半径・軸対称ソースヤコビアン・dual-time 物理 BDF 項 (非前処理) も倍精度で踏襲。
- **`SLAU_d`** は `lowMachPrecond>=1` で `c'` 散逸を使う (==2 も散逸是正を併用)。$\beta=1$ で $\Gamma_c=I$・$\Delta\tau'=\Delta\tau$・
  フラックス同一ゆえ現行カーネルと同一組み立て (倍精度の丸め差 ~1e-7、解一致)。

> **検証結果 (2026-06-09・採用)**: `case/23.axi_nozzle` で**低マッハ自励振動を根治**。Phase1 (物理 LHS) が発散した
> $\epsilon=0.05$ を前処理 LHS が安定化し、chamber 圧振幅 (M<0.08, 4k–20k) を 0.882%→**0.087% (定常収束・振動消滅)**。
> 安定 `cfl_pseudo` も m1~1→m2~5-7 と拡大 (収束加速は per-step 2.54× で等 wall-clock 互角。価値は根治)。データは計画 §9。

### `implicit_defect_correction_d`（scalar 対角版＝スペクトル半径）

平均流 5 式を、5×5 ブロックの代わりに**スカラー対角**で解く軽量版（`blockDPLUR == 0`）。
各セルで対角 $D = V/\Delta\tau + \sum_f(|U_n|+c+\rho^\nu)S_f$、off-diagonal $\tfrac12(|U_n|+c+\rho^\nu)S_f$ で
`dq_*_new = relax·(res_* + Σ offdiag·dq_*_nbr^{old}) / D`（5 変数とも同じスカラー $D$）を作り、
[`blockDPLURSolve`](../../solver_density_cuda/main.cpp) が `nStepInner` 回 sweep（各 sweep 後に
`swapScalarImplicitCorrectionBuffers`）、最後に [`applyScalarImplicitCorrection`](../../solver_density_cuda/update.cpp)
で `Q = Q_N + dq_*_old` を 1 回 commit する（block 版と同じ古典 DPLUR 制御フロー）。
スペクトル半径は符号不変なので block の $A^\pm$ 法線符号問題は持たない。

**block 版との比較（2026-06, `case/20.naca_ml`）**: 収束先は同一（収束場から再開すると同じ解を保持、
壁面静圧 平均 0.02% 一致）。ただし近似ヤコビアンが粗いため**安定 `cfl_pseudo` が大幅に低い**
（scalar ≲ 1〜2、block は 20〜50）。supercritical 始動では block が cfl_pseudo=20 で 4000 step / 25s で
roe→0.5 に収束する一方、scalar は cfl_pseudo=1 で 12000 step でも収束せず大きな過渡オーバーシュート
（roe ピーク ~186 vs block/explicit ~85）を示す。よって **block DPLUR が既定、scalar 対角版は
5×5 を避けたい軽量用途・低レジスタ用途向けのフォールバック**と位置づける。

### `applySSTPointImplicit_d`（SST k-ω の segregated point-implicit）

平均流 commit 後に呼ぶスカラー陰解法。源項+輸送項のヤコビアン対角を陰化して k/ω の凍結を解き、
壁近傍の陽的輸送 stiff 性を緩和して安定 `cfl_pseudo` を一桁以上引き上げる。

- 消散ヤコビアン `src_jac_k=β*ω`・`src_jac_omega=2βω` は
  [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) `rans_sst_source_d` が出力。
- 輸送ヤコビアン `transport_diag_k`/`transport_diag_omega` [m³/s] は
  [`scalarTransport_d.cu`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu) の `scalar_advection_first_order_d`
  （1次風上 $\sum_f\max(\pm\dot m,0)/\rho$）と `scalar_diffusion_first_order_d`（$\sum_f(\mu_{\text{face}}/\rho)|\delta|/dcc$）が
  面ループで atomicAdd 集計する（`ransTransport_d_wrapper` 冒頭で毎 `assembleResidual` ゼロ初期化）。
- [`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) の `applySSTPointImplicit_d` が各セルで
  $D_\phi = V/\Delta\tau + V\cdot\text{src\_jac}_\phi + \text{transport\_diag}_\phi$、
  $\delta(\rho\phi)=\text{relax}\cdot\text{res}_{\rho\phi}/D_\phi$、$\rho\phi=\max(\rho\phi^{N}+\delta,\ \text{floor})$ を適用
  （$\rho k\ge0$, $\rho\omega>0$）。`dt_local` は平均流と共用。
- [`main.cpp`](../../solver_density_cuda/main.cpp) `implicitNonlinearUpdate` で `scalarResidualEnabled` のとき
  `applyBlockImplicitCorrection` 直後に `applySSTPointImplicit` を呼ぶ。生産・近傍 ΔQ は `res` に含む lagged。
- 近傍 ΔQ 結合なしの純 point-implicit（消散+輸送の対角のみ陰化）。defect-correction のため定常解は不変。
  理論・検証は [theory.md](theory.md) §"輸送項 (移流+拡散) の point-implicit 対角"。

## 局所時間刻み

`setDT_d_wrapper` ([`setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu)) が
セル中心スペクトル半径を集計して `dt_local[ic] = CFL * V / λ_max` を書き込む。
`cfg.dt`, `cfg.cfl` で挙動を制御。

**setDT の低マッハ前処理は不採用** (2026-06 検証)。`setCFL_pln_d` のスペクトル半径の音速 `sonic` を
前処理音速 `c'` に置換する案は、低マッハ域で `dt_local` を増大させ陰解法対角 `V/Δτ` を縮め、block DPLUR の
対角優位性を崩して発散させた。よって `setCFL_pln_d` は従来の `sonic` のまま。低マッハ前処理は対流フラックス
([`convectiveFlux_d.cu`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) `SLAU_d`) の散逸スケールにのみ
適用する。詳細は計画 [`time_integration-lowmach-preconditioning.md`](../../.github/plans/time_integration-lowmach-preconditioning.md) §9。

## 入出力

入力: 残差 `res_ro, res_roUx, …`、$\Delta t_{\text{loc}}$ (`dt_local`)、
過去ステージの `Q_N, Q_M`、対角ヤコビアン構築用の `Ux, Uy, Uz, sonic, Ht, vis_turb`、
ステージ係数 (`cfg.coef_*`)。

出力: 更新後の `Q = (ro, roUx, roUy, roUz, roe)`。
陰解法では補助バッファ `dq_*_old/new`, `dq_block_old/new_k`, `diag_block_*`, `rhs_block_k`。

## 非定常 dual-time 陰解法（実装済み 2026-06）

`advanceImplicitDualTime` ([`main.cpp`](../../solver_density_cuda/main.cpp)) が 1 物理ステップを担当する。
使用条件 `unsteady=1, dualTime=1, timeIntegration=11, blockDPLUR=1, time.deltaT.control=0`（それ以外は throw）。

1 物理ステップの流れ:
1. **時間レベルシフト** `shiftDualTimeLevels_d_wrapper`: `roNN←roN`, `roN←ro`（$\mathbf Q^{n-1}\!\leftarrow\!\mathbf Q^n$, $\mathbf Q^n\!\leftarrow$現在）。
2. **BDF 係数**: 初回ステップ or `bdfOrder==1` は BDF1 $(a,b,c)=(1,1,0)$、以降 BDF2 $(\tfrac32,2,\tfrac12)$。
   `cfg.unsteadyDiagCoef = a/\Delta t` を設定（陰解法カーネルが対角に $V\cdot$この係数を加える）。
3. **擬似時間サブ反復** `nSubIterDualTime` 回:
   - `assembleResidual` → `addUnsteadyTimeTerm_d_wrapper(a,b,c)` で `res_* -= (V/\Delta t)(a\mathbf Q - b\mathbf Q^n + c\mathbf Q^{n-1})`
     （`include_scalar` で k/ω も）。
   - `setDT`（擬似 $\Delta\tau$）→ `blockDPLURSolve`（対角に $V\,a/\Delta t$ 込み）。
   - **in-place commit** `applyBlockImplicitCorrectionInPlace_d_wrapper`（$\mathbf Q\mathrel{+}=\delta\mathbf Q$。`roN`=$\mathbf Q^n$ 固定のため）。
     RANS のとき `applySSTPointImplicit`（こちらも in-place）。
4. `cfg.totalTime += \Delta t`、`unsteadyDiagCoef=0` リセット、出力。

実装した CUDA（[`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu)）: `addUnsteadyTimeTerm_d`,
`applyBlockImplicitCorrectionInPlace_d`, `shiftDualTimeLevels_d`。block/scalar/SST 各カーネルに対角の物理時間係数
`unsteady_diag` (`cfg.unsteadyDiagCoef`) を追加。第 2 時間レベル `roNN`/`roKNN` 系は [`variables.hpp`](../../solver_density_cuda/variables.hpp) に登録済。

## 既知の TODO / 注意点

- 非定常 dual-time 陰解法（`tI==11 && unsteady==1 && dualTime==1`）は実装済（2026-06、`blockDPLUR==1` のみ、物理 $\Delta t$ 固定 `control=0`）。`implicitCorrection_d.cu` の `dualtime_explicit_d` は SLAU/Roe 用の別系統補助で本流とは独立（未使用）。
- scalar 対角陰解法（`tI==11 && blockDPLUR==0`）は有効（2026-06）。block より低 `cfl_pseudo` で収束も遅いため既定は block（上記比較参照）。
- block 陰解法でも SST(k/ω) は `applySSTPointImplicit` で segregated point-implicit 更新され、**凍結しない**（2026-06）。消散+輸送(移流+拡散)の対角を陰化、生産・近傍 ΔQ は lagged。
- `matrix mat_ns` は陰解法では未使用だが非陰解法のシグネチャに残るため `StepContext` に保持（除去は別途）。
- 旧 CPU の [`update.cpp`](../../solver_density_cuda/update.cpp) は使用されていない。
- 局所時間刻みの粘性スペクトル半径寄与は `setDT_d` 側で個別に実装されている。詳細は同ファイルを参照。

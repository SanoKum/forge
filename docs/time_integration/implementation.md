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
  → convectiveFlux → scalarTransport → calcScalarGradient + ransSource
  → axisymmetricSource → viscousFlux
  → [addUnsteadyTimeTerm(ctx)]            // dual-time の BDF 物理時間項フック（定常は no-op）

advanceExplicitRK(ctx):                    // tI 1/3/4（挙動不変）
  for iloop in perStepIterationCount():
     updateVariablesInner; assembleResidual(ctx,iloop+1); logResidualSnapshot
     timeIntegration_d_wrapper(iloop); scalarTimeIntegration_d_wrapper(iloop)
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

## 陰解法カーネル詳細

### `implicit_defect_correction_block_d`（block DPLUR、LU-SGS $A^\pm$ 分割）

5×5 ブロック版。`block_dplur` 名前空間に補助 device 関数群を持つ。

- `build_jacobian_split` — セル中心状態から固有分解 $R,\Lambda,L$ を作り、**対角用 $A^{+}=R\,\Lambda^{+}L$** と
  **RHS 近傍用 $K=-A^{-}=R\,(-\Lambda^{-})L=\tfrac12(|\widetilde A|-\widetilde A)$** を同時に返す
  （$\Lambda^{+}=\max(\Lambda,0)$、$-\Lambda^{-}=\max(-\Lambda,0)$、共に非負）。
- `add_identity_scaled`, `add_scaled_5x5` — $V/\Delta\tau\,I$・粘性対角・面寄与の加算。
- `solve_5x5` — 部分ピボット付き Gauss 消去。`|pivot| < 1e-20` でゼロ解にフォールバック。
- `multiply_add_5x5_vec` — 行列ベクトル積。

各セルで近傍寄与を集約し、対角 $D_i = V/\Delta\tau\,I + \sum_f A^{+}_f S_f + \sum_f \rho^{\nu}_f S_f\,I$、
RHS = $-\mathbf R + \sum_f K_f S_f \cdot \Delta\mathbf Q_{\text{nbr}}^{\text{old}}$（$K_f=-A^{-}_f$）を構築し、
$\Delta\mathbf Q_{\text{new}} = D_i^{-1}\,\text{RHS}$ を解く。`cfg.implicitRelax` で $\Delta\mathbf Q$ を緩和。

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

平均流 commit 後に呼ぶスカラー陰解法。消散項のヤコビアン対角を陰化して k/ω の凍結を解く。

- [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) の `rans_sst_source_d` が
  残差 `res_roK`/`res_roOmega` に加え、消散ヤコビアン `src_jac_k = β*ω`・`src_jac_omega = 2βω` を出力。
- [`update_d.cu`](../../solver_density_cuda/cuda_forge/update_d.cu) の `applySSTPointImplicit_d` が
  各セルで $D_\phi = V/\Delta\tau + V\cdot\text{src\_jac}_\phi$、$\delta(\rho\phi)=\text{relax}\cdot\text{res}_{\rho\phi}/D_\phi$ を作り
  $\rho\phi = \max(\rho\phi^{N}+\delta,\ \text{floor})$ を適用（$\rho k\ge0$, $\rho\omega>0$）。`dt_local` は平均流と共用。
- [`main.cpp`](../../solver_density_cuda/main.cpp) `implicitNonlinearUpdate` で `scalarResidualEnabled` のとき
  `applyBlockImplicitCorrection` 直後に `applySSTPointImplicit` を呼ぶ。移流・拡散・生産は `res` に含む lagged。
- 近傍結合なしの純 point-implicit（消散 stiff 性のみ陰化する最小構成）。理論は [theory.md](theory.md) §"スカラー (k/ω) の陰解法"。

## 局所時間刻み

`setDT_d_wrapper` ([`setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu)) が
セル中心スペクトル半径を集計して `dt_local[ic] = CFL * V / λ_max` を書き込む。
`cfg.dt`, `cfg.cfl` で挙動を制御。

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
- block 陰解法でも SST(k/ω) は `applySSTPointImplicit` で segregated point-implicit 更新され、**凍結しない**（2026-06）。消散項のみ陰化、移流・拡散・生産は lagged。
- `matrix mat_ns` は陰解法では未使用だが非陰解法のシグネチャに残るため `StepContext` に保持（除去は別途）。
- 旧 CPU の [`update.cpp`](../../solver_density_cuda/update.cpp) は使用されていない。
- 局所時間刻みの粘性スペクトル半径寄与は `setDT_d` 側で個別に実装されている。詳細は同ファイルを参照。

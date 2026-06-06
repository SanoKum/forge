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
| `4` | `runge_kutta_exp_4th_d` ([L226](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L226)) |
| `1`, `3` | `runge_kutta_exp_d` ([L317](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L317)) |
| `11`, `blockDPLUR == 0` | `implicit_defect_correction_d` ([L372](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L372)) |
| `11`, `blockDPLUR == 1` | `implicit_defect_correction_block_d` ([L516](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu#L516)) |
| `10` (TODO) | `runge_kutta_dual_explicit_d` (コメントアウト) |

## ループ全体 ([`main.cpp`](../../solver_density_cuda/main.cpp))

```
updateVariablesOuter_d_wrapper()        // Q_N ← Q, Q_M ← Q
for iInner in 0 .. nInner-1:
    applyBconds()
    gradient + limiter + ducros
    convectiveFlux + viscousFlux       // res_* に累積
    [implicitCorrection_d_wrapper()]   // 任意
    timeIntegration_d_wrapper(iInner)  // Q 更新
    updateVariablesInner_d_wrapper()   // Q_M ← Q
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

### `implicit_defect_correction_d`

スカラ近似版。各セルで法線スペクトル半径から擬時間係数を作り、
`dq_*_new = (-res - implicit_off_diag(dq_old)) / (V/Δτ + spectral_radius)` 形式の
更新を 1 回行う。`loop == 0` で `dq_*_old` をゼロクリア。
反復終了後、[`applyScalarImplicitCorrection_d`](../../solver_density_cuda/cuda_forge/update_d.cu#L113)
で `Q = Q_N + dq_*_old` を適用する。

### `implicit_defect_correction_block_d`

5×5 ブロック版。`block_dplur` 名前空間に補助 device 関数群を持つ。

- `build_abs_jacobian` — Roe 平均状態から $|\widetilde A| = R |\Lambda| L$ を構築。
- `add_identity_scaled`, `add_scaled_5x5` — $V/\Delta\tau\,I$ や近傍寄与の加算。
- `solve_5x5` — 部分ピボット付き Gauss 消去。`|pivot| < 1e-20` でゼロ解にフォールバック。
- `multiply_add_5x5_vec` — 行列ベクトル積。

各セルで近傍 5 セルの寄与を集約し、LHS = $V/\Delta\tau\,I + \sum_f |\widetilde A_f| S_f / 2$、
RHS = $-\mathbf R + \sum_f |\widetilde A_f| S_f / 2 \cdot \Delta\mathbf Q_{\text{nbr}}^{\text{old}}$ を構築して解く。

`cfg.implicitRelax` で更新量を緩和。`swapImplicitCorrectionBuffers(var)` で
`dq_*_old <-> dq_*_new` をスワップして反復を進める。

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

## 既知の TODO / 注意点

- `timeIntegration == 10` (dual time + 陽 RK) は未実装 (`runge_kutta_dual_explicit_d` がコメントアウト)。
- `implicitCorrection_d.cu` の `dualtime_explicit_d` は SLAU/Roe 用補助の dual-time 補正で、
  本流の time integration とは独立。
- 旧 CPU の [`update.cpp`](../../solver_density_cuda/update.cpp) は使用されていない。
- 局所時間刻みの粘性スペクトル半径寄与は `setDT_d` 側で個別に実装されている。
  詳細は同ファイルを参照。

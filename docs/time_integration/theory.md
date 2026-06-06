# 時間積分 — 理論

forge の時間積分は次の 3 系統を提供する。

| `cfg.timeIntegration` | 概要 |
| --- | --- |
| `1` / `3` | Jameson 型多段陽解法 (擬時間ステップ) |
| `4` | 4 段 4 次 Runge–Kutta (低 storage、係数は `coef_DT_4thRunge`, `coef_Res_4thRunge`) |
| `11` | 陰解法 defect-correction (対角ヤコビアン近似)、`blockDPLUR == 1` で 5×5 ブロック版 |

## 共通枠組み

時間 1 ステップは内部ループ (`loop = 0 .. nInner-1`) に分かれ、
1 ループあたり次の手順を踏む。

```
applyBconds → gradient → limiter → ducrosSensor
   → convectiveFlux + viscousFlux  (res_* に累積)
   → implicitCorrection (任意)
   → timeIntegration_d_wrapper (Q を更新)
   → dependentVariables (派生量)
```

`res_*` は対流・粘性両方からの符号付き面寄与の合計。これを各セルで体積で
正規化したものが $-\partial \mathbf{Q}/\partial t \cdot V$ に相当する。

## 陽解法 (Jameson 多段)

低 storage 形式: 1 ステージあたり

$$
\mathbf{Q}^{(k+1)} = c_N \mathbf{Q}^N + c_M \mathbf{Q}^M
+ c_R \frac{\Delta t_{\text{loc}}}{V}\,\mathbf{R}^{(k)},
$$

ここで $\mathbf{Q}^N$ は前回外部ステップ、$\mathbf{Q}^M$ は前回内部ステージの値。
係数 $c_N, c_M, c_R$ は `solverConfig` の `coef_N`, `coef_M`, `coef_Res` に格納する
(段数 = `cfg.nInner`)。

定常計算では局所時間刻み $\Delta t_{\text{loc}}$ を使うことで収束を加速する。

## 4 段 4 次 Runge–Kutta (low-storage)

参考: Allampalli et al. (2009)。`coef_DT_4thRunge`, `coef_Res_4thRunge` を用い、
内部 4 段で

$$
\mathbf{Q}^{(k+1)} = \mathbf{Q}^N + c_{DT}^{(k)} \frac{\Delta t_{\text{loc}}}{V}\,\mathbf{R}^{(k)},
$$

累積残差 $\mathbf{R}_m = \sum_k c_R^{(k)} \mathbf{R}^{(k)} \,\Delta t_{\text{loc}}/V$ を保持し、
最終段で $\mathbf{Q}^{N+1} = \mathbf{Q}^N + \mathbf{R}_m$ とする。LES など時間精度の必要な解析向け。

## 陰解法 defect-correction

擬時間 $\tau$ を導入した擬時間スキームで、$\Delta \mathbf{Q} = \mathbf{Q}^{n+1} - \mathbf{Q}^n$ を未知に取る。
1 次風上に基づく対角ブロック近似ヤコビアン $\widetilde{A}$ で

$$
\left[\frac{V}{\Delta \tau} I + \widetilde{A}\right] \Delta \mathbf{Q}
= -\mathbf{R}(\mathbf{Q}^n) + (\text{defect})
$$

を、サブイタレーションを 1 回ずつ反復しながら更新する (Gauss–Seidel 系)。
`implicitRelax` で緩和量を制御。`blockDPLUR == 1` のとき 5×5 ブロック実装
(`implicit_defect_correction_block_d`) を用い、フル 5×5 行列を各セルで構築・LU 分解する。

### 対角ヤコビアン構築

`build_abs_jacobian` でセル中心の Roe 平均状態から $|\widetilde{A}| = R \,|\Lambda|\,R^{-1}$ を
構築 (右/左固有ベクトル `R, L`)。固有値は $\{U \pm c,\, U,\, U,\, U\}$ の絶対値。
これに $V/\Delta\tau$ をスケールした単位行列を足したものを LHS とする。

### 5×5 LU 分解

`solve_5x5` で部分ピボット付き Gauss 消去を実装。ピボットが
`1e-20` 未満の場合はゼロ解を返してフェイルセーフ動作する。

## 局所時間刻み

`setDT_d_wrapper` (`cuda_forge/setDT_d.cu`) が各セルの $\Delta t_{\text{loc}}$ を
CFL 条件から決定する。スペクトル半径

$$
\lambda_{\max} = \sum_f \big( |U_n| + c \big)_f A_f + \text{粘性項}
$$

から $\Delta t_{\text{loc}} = \mathrm{CFL}\,V / \lambda_{\max}$ を求める。
`cfg.dt` が指定された場合はそれを上限とする。

## 参考

- 各カーネル分岐の詳細: [implementation.md](implementation.md)。
- 境界条件適用 (`applyBconds`): [`docs/boundary/`](../boundary/)。
- 残差を作る対流・粘性: [`docs/convection/`](../convection/), [`docs/diffusion/`](../diffusion/)。

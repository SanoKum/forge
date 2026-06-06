# 粘性・熱伝導フラックス — 理論

forge は Navier–Stokes の粘性応力と熱伝導束を、セル中心勾配の面平均と
セル中心値の法線差分を組み合わせた "over-relaxed" スキームで離散化する。
実装の対応は [implementation.md](implementation.md) を参照。

## 連続方程式系

粘性応力テンソル (Newton 流体、Stokes 仮定) と熱伝導束:

$$
\tau_{ij} = \mu \!\left( \frac{\partial u_i}{\partial x_j}
+ \frac{\partial u_j}{\partial x_i} \right)
- \frac{2}{3}\mu (\nabla \cdot \mathbf{u})\,\delta_{ij},
\qquad
q_i = -\kappa\, \frac{\partial T}{\partial x_i}.
$$

セル境界フラックスは

$$
\mathbf{F}^{\text{visc}} =
\begin{pmatrix}
0\\
\tau_{xj}\,n_j\\
\tau_{yj}\,n_j\\
\tau_{zj}\,n_j\\
\tau_{ij} u_i n_j + \kappa\, (\partial_j T)\, n_j
\end{pmatrix}.
$$

## Over-relaxed 離散化

セル $c_0, c_1$ をまたぐ面 $f$ で速度・温度の勾配を求めるとき、
セル中心間ベクトル $\mathbf{d} = \mathbf{r}_{c_1} - \mathbf{r}_{c_0}$ と面ベクトル
$\mathbf{S}$ の非直交を補正する。$|\mathbf{S}| = S$ として

$$
\boldsymbol{\Delta} = \frac{|\mathbf{S}|^2}{\mathbf{d}\cdot\mathbf{S}}\,\mathbf{d},
\qquad \mathbf{k} = \mathbf{S} - \boldsymbol{\Delta}.
$$

法線方向成分はセル中心差分から、接線成分は両側セルの勾配を面平均
($f_x$ 重み) して構成する:

$$
(\partial_j u_i)\, S_j \;\approx\;
\underbrace{\frac{u_i^{c_1} - u_i^{c_0}}{|\mathbf{d}|} \, \Delta_j}_{\text{法線差分}}
\;+\;
\underbrace{\overline{(\partial_j u_i)}_f \, k_j}_{\text{接線補正}},
$$

$\overline{(\cdot)}_f = f_x (\cdot)_{c_0} + (1-f_x)(\cdot)_{c_1}$。

この組み合わせを直接 $\tau$ に代入することで、粘性束を面ごとに評価し、
両側セルに $\pm$ で加算する。エネルギ束も同様に $\tau_{ij} u_i$ と
$\kappa (\partial_j T)$ を組み合わせる。

## 粘性係数

- 層流粘性 $\mu_{\text{lam}}$: 設定値 (`cfg.visc`、定数) を使用。
- 乱流粘性 $\mu_{\text{turb}}$: `vis_turb[ic]` (SGS / RANS モデル別ファイルで生成)。
- 有効粘性 $\mu = \mu_{\text{lam}} + \mu_{\text{turb}}$。
- 熱伝導率 $\kappa$: 設定値 `cfg.thermCond` を使用 (Pr 数からの導出は呼び出し側)。

## 壁面寄与

非滑り条件下の壁面では速度がゼロとなる。壁面用カーネル `viscousFlux_wall_d` で
壁面距離 (ghost cell までの法線距離) を使って粘性応力を片側差分で評価する。
これにより内部面ループでは未処理の壁面寄与を補う。

## 参考

- 対流フラックスとの加算: [`docs/convection/`](../convection/) も参照。
  両者の残差は同じ `res_*` バッファに累積される。
- 勾配の生成: [`docs/gradient/`](../gradient/) を参照。

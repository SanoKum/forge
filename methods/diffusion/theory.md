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

> **重要 — 法線項はスカラー係数で評価する。** 上式の法線差分項
> $(u_i^{c_1}-u_i^{c_0})/|\mathbf{d}|\cdot\Delta_j$ で各成分 $i$ に掛かるのは、
> over-relaxed の**スカラー** $\beta = |\mathbf{S}|^2/(\mathbf{d}\cdot\mathbf{S})$ である。
> すなわち $\mu\,\nabla u_i\cdot\mathbf{S}$ の法線寄与は
> $\mu\,\beta\,(u_i^{c_1}-u_i^{c_0})$(全成分共通のスカラー $\beta$)であり、
> 熱伝導束 $\kappa\,\beta\,(T^{c_1}-T^{c_0})$ と同型である。
> 速度成分ごとに $\boldsymbol\Delta$ の**自身の座標成分** $\Delta_i$ を掛けると
> (例: $\tau_x$ に $\Delta_x$)、軸平行な $y$ 法線面で $\Delta_x=0$ となり
> 流れ方向運動量の横方向拡散 $\mu\,\partial u_x/\partial y$ が落ちる。
> Laplacian 接線補正 $\overline{\nabla u_i}_f\cdot\mathbf{k}$ は**同一成分** $u_i$ の勾配
> $\overline{\partial u_i/\partial x_j}\,k_j$ を用いる。転置寄与 $\mu\,\partial u_j/\partial x_i\,S_j$ は
> 別項として面平均勾配にフル $\mathbf{S}$ を内積して加える。実装は
> [implementation.md](implementation.md) を参照。

## 粘性係数

- 層流粘性 $\mu_{\text{lam}}$: 設定値 (`cfg.visc`、定数) を使用。
- 乱流粘性 $\mu_{\text{turb}}$: `vis_turb[ic]` (SGS / RANS モデル別ファイルで生成)。
- 有効粘性 $\mu = \mu_{\text{lam}} + \mu_{\text{turb}}$。
- 熱伝導率 $\kappa$: **有効値 $\kappa_{\text{eff}} = \kappa_{\text{lam}} + c_p\,\mu_{\text{turb}}/Pr_t$**。
  $Pr_t$ は `turbulence.turbulentPrandtl` で設定 (既定 0.85)。層流 $\kappa_{\text{lam}}$ は設定値
  `cfg.thermCond`、乱流寄与は $\mu_{\text{turb}}$ から渦熱伝導として加える。$c_p$ は
  thermally-perfect の温度依存を反映するためセル配列 (`var.c_d["cp"]`、`thermalMethod==2` で $c_p(T)$)
  を面平均して使う。応力が $\mu=\mu_{\text{lam}}+\mu_{\text{turb}}$ を使う以上、熱伝導も対応する乱流寄与を
  含めないと乱流境界層でエネルギーが保存せず、断熱壁の静温が回復温度 ($\le$ 全温) を超えて発散的に上昇する。

## 壁面寄与

非滑り条件下の壁面では速度がゼロとなる。壁面用カーネル `viscousFlux_wall_d` で
壁面距離 (ghost cell までの法線距離) を使って粘性応力を片側差分で評価する。
これにより内部面ループでは未処理の壁面寄与を補う。

壁面でも粘性束は内部面と**同じ** $\tau_{ij} n_j$ の評価でなければならない。
すなわち運動量 $i$ 成分に対し

$$
F^{\text{visc}}_i = \tau_{ij}\, S_j
= \mu\!\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) S_j
- \frac{2}{3}\mu (\nabla\cdot\mathbf{u})\, S_i .
$$

平板チャネル等の軸平行壁 (法線 $\mathbf{n}=\pm\mathbf{e}_y$) では、壁摩擦の主項は
$\tau_{xy} S_y = \mu\,(\partial u_x/\partial y)\, S_y$ という**せん断 (off-diagonal) 成分**であり、
これがストリーム方向運動量 $\rho u_x$ への no-slip 抗力を与える。
したがって法線差分 $(u_i^{g}-u_i^{c})/|\mathbf{d}|$ は内部面と同じく
over-relaxed 係数 $\Delta_i$ (直交格子では $\Delta_i = S_i$、すなわち面積 $S$ に法線
単位ベクトルを掛けたもの) と組み合わせ、$\tau_{ij}S_j$ として評価する必要がある。
速度成分ごとに自身の座標成分 $S_i$ のみを掛ける ($\tau_x \propto S_x$ 等) と、
軸平行壁では主せん断項が落ちて壁摩擦がゼロになる (この不具合は 2026-06-06 に修正済み。
実装は [implementation.md](implementation.md) を参照)。

## 参考

- 対流フラックスとの加算: [`methods/convection/`](../convection/) も参照。
  両者の残差は同じ `res_*` バッファに累積される。
- 勾配の生成: [`methods/gradient/`](../gradient/) を参照。

# リミッタ — 理論

MUSCL 再構成 (2 / 3 次) は不連続近傍で振動を生む。
forge ではセル中心勾配にスカラ係数 $\psi \in [0, 1]$ を掛けて再構成を抑制する。

$$
\phi_L = \phi_C + \psi_C \,\nabla \phi_C \cdot (\mathbf{r}_f - \mathbf{r}_C).
$$

実装は [implementation.md](implementation.md) を参照。

## 評価対象

セルごとに $\rho, U_x, U_y, U_z, P$ の 5 成分について独立にリミッタ係数を持つ。
温度 $T$ と全エンタルピ $H_t$ は派生量として扱い、独立リミッタは持たない。

## $\delta^+, \delta^-, \delta_m$ の定義

セル $C$ について隣接セル値の振れ幅を

$$
\delta^+_C = \phi_C^{\max} - \phi_C ,\qquad
\delta^-_C = \phi_C^{\min} - \phi_C
$$

($\phi_C^{\max/\min}$ は近傍セル値の最大・最小)。
面 $f$ への線形再構成による予測差分を

$$
\delta_m = \tilde\phi_f - \phi_C
= \nabla \phi_C \cdot (\mathbf r_f - \mathbf r_C)
$$

とする。$\delta_m > 0$ なら $\delta^+$ を、$\delta_m < 0$ なら $\delta^-$ を比較対象に取る。

## Barth–Jespersen リミッタ

$$
\psi_f = \min\!\left(1,\, \frac{\delta^\pm}{\delta_m}\right),
\qquad \psi_C = \min_{f \in \partial C} \psi_f .
$$

単純で TVD 性が強く、滑らかな領域でも 2 次精度の頭打ちを招きやすい。

## Venkatakrishnan リミッタ

Barth–Jespersen を滑らか化したもの。$\epsilon^2 = (K |V_C|^{1/3})^3$ (`K=1.0`) として

$$
\psi(\delta^+, \delta_m) =
\frac{1}{\delta_m}\,
\frac{(\delta^{+\,2} + \epsilon^2)\,\delta_m + 2\delta_m^2 \delta^+}
     {\delta^{+\,2} + 2\delta_m^2 + \delta^+ \delta_m + \epsilon^2}.
$$

$\delta_m \to 0$ で $\psi \to 1$ (連続)。$|\delta_m|$ が体積スケール以下では実質
無リミットになり、滑らかな領域での精度低下を避ける。

## Nishikawa R1 リミッタ (未有効)

`nishikawa_r1_limiter` の実装は残されているがコメントアウト済み。
$\delta'_C$ と幾何比 $r_{ik}$ から不連続適応する形式。

## Ducros センサとの併用

[`docs/convection/`](../convection/) で記述したように、
Roe カーネルでは Ducros センサ値 $\mathrm{duc}$ で

$$
\psi \leftarrow \begin{cases}
\psi & (\mathrm{duc} \le 0.8) \\
\max(0,\, (1 - \mathrm{duc})\,\psi) & (\mathrm{duc} > 0.8)
\end{cases}
$$

と二段がけする。

## 全リミッタ 1 の場合 (`cfg.limiter == 0`)

リミッタを評価せず $\psi = 1$ を一律で配るモード。低次再構成 (`convMethod = 0`)
や smooth な計算ですべての面で線形再構成を使いたい場合に選ぶ。

## 参考

- リミッタの参照側 (`interp_dispatch`): [`docs/convection/implementation.md`](../convection/implementation.md)。
- 勾配の生成: [`docs/gradient/`](../gradient/)。

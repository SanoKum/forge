# SST 乱流モデルの理論メモ

## 1. 目的

本章では forge に導入する RANS ベースの **Menter SST**
($k$-$\omega$ Shear Stress Transport) モデルの支配方程式と、
本コードで採用する前提を整理する。

最終的な導入対象は SST 全系だが、初回マイルストーンは
scalar transport 配線確認のため `k`, `\omega` の advection-only とする。

初期スコープは次のとおりとする。

- 低 Re の壁解像型 SST
- 圧縮性 NS ソルバ上の保存形輸送
- まず generic scalar transport 基盤上で `k`, `\omega` を扱う
- 初回は advection-only を対象とし、その後 diffusion / source を追加する
- 軸対称は既存幾何の上で advection を確認し、SST 固有幾何項は別 plan で拡張
- 陽解法 (`timeIntegration = 1, 3, 4`) を対象とし、陰解法は別 plan

## 2. 保存変数

既存の 5 保存変数

$$
Q_{NS} = [\rho, \rho u_x, \rho u_y, \rho u_z, \rho E]^T
$$

に加えて、SST では次の 2 変数を導入する。

$$
Q_{turb} = [\rho k, \rho \omega]^T
$$

ここで $k$ は乱流運動エネルギー、$\omega$ は比散逸率である。

## 3. 輸送方程式

圧縮性保存形では、2 つの追加方程式を generic scalar transport の特殊例として扱う。
最終形は次である。

$$
\frac{\partial (\rho k)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} k)
= \nabla \cdot \left[
\left(\mu + \sigma_k \mu_t\right) \nabla k
\right]
+ P_k - \beta^* \rho \omega k
$$

$$
\frac{\partial (\rho \omega)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} \omega)
= \nabla \cdot \left[
\left(\mu + \sigma_\omega \mu_t\right) \nabla \omega
\right]
+ \alpha \frac{\omega}{k} P_k
- \beta \rho \omega^2
+ D_{\mathrm{cross}}
$$

ここで $P_k$ は乱流生成項、$D_{\mathrm{cross}}$ は SST 特有の
cross-diffusion 項である。

初回マイルストーンでは右辺を落とし、

$$
\frac{\partial (\rho \phi)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} \phi) = 0
$$

で $\phi \in \{k, \omega\}$ を輸送し、保存量配線と境界条件・時間積分を先に確立する。

### 3.1 source term の全体像

本来の SST では、対流・拡散だけでなく生成項と散逸項が入る。実装上は
`k` と `\omega` を独立な scalar と見なしつつ、source は同じ 2 本の式へ加える。

標準形は次である。

$$
S_k = P_k - \beta^* \rho \omega k
$$

$$
S_\omega = \alpha \frac{\omega}{\max(k, k_{\min})} P_k
- \beta \rho \omega^2
+ D_{\mathrm{cross}}
$$

ここで

$$
P_k = \min\left(\tau_{ij} \frac{\partial u_i}{\partial x_j},\; 10\,\beta^*\rho\omega k\right)
$$

であり、$\tau_{ij}$ は Reynolds stress を近似する eddy-viscosity closure から作る。
SST では stagnation point などで $P_k$ が過大になりやすいため、上式のように
production limiter を入れるのが基本である。

cross-diffusion 項は、SST が $k$-$\varepsilon$ 側と $k$-$\omega$ 側を混ぜるために入る補正で、
一般には

$$
D_{\mathrm{cross}} = 2 (1 - F_1) \rho \sigma_{\omega 2}
\frac{1}{\max(\omega, \omega_{\min})}
\nabla k \cdot \nabla \omega
$$

と書ける。`F_1` は壁近傍で 1 に寄る blending function で、壁から離れるほど 0 に近づく。
したがって cross-diffusion は壁近傍では消え、自由流側で効いてくる。

### 3.2 blending による係数切替

SST では係数を 2 組持ち、$F_1$ で混合する。

$$
\sigma_k = F_1 \sigma_{k1} + (1 - F_1) \sigma_{k2},
\qquad
\sigma_\omega = F_1 \sigma_{\omega1} + (1 - F_1) \sigma_{\omega2}
$$

$$
\alpha = F_1 \alpha_1 + (1 - F_1) \alpha_2,
\qquad
\beta = F_1 \beta_1 + (1 - F_1) \beta_2
$$

`F_1` と `F_2` は壁距離 $y$、$k$、$\omega$、勾配を使って定義される。
forge では既存の `wall_dist` をその入力として再利用する。

source と blending の実装上の注意は次のとおりである。

- $k$ と $\omega$ は正値制約を壊しやすいので、source でも下限値を入れる
- $\omega$ の発散は壁近傍で特に起こりやすいため、$\omega_{\min}$ を明示する
- production limiter と blending を分離して実装し、検証しやすくする

## 4. 渦粘性

SST の渦粘性は、$k$-$\omega$ 系をベースに blending function を通して
壁近傍と自由流側を切り替える。

代表形は次である。

$$
\mu_t = \rho \nu_t,
\qquad
\nu_t = \frac{a_1 k}{\max(a_1 \omega, S F_2)}
$$

ここで $S$ はひずみ速度の代表量、$F_2$ は壁距離を含む blending function
である。forge では既存の `wall_dist` を再利用し、追加の壁距離ソルバは
導入しない。

## 5. 係数と blending

SST では $k$-$\omega$ 壁近傍モデルと、自由流側の $k$-$\varepsilon$ 的な
振る舞いを blending function $F_1$, $F_2$ で接続する。

実装では次を守る。

- 係数はコード内に固定値として持つ初期実装とする
- 将来、YAML から係数上書きを許す場合も、既定値は Menter SST を維持する
- $k$, $\omega$, $\mu_t$ の positivity を壊さないために下限値を設ける

## 6. 壁条件

初期スコープは壁関数なしの壁解像型とする。SST の壁条件は、
`k` と `\omega` を壁で別々に閉じる必要がある。

### 6.1 no-slip wall / wall_isothermal

最小構成では、壁面で乱流エネルギーを消し、比散逸率は壁距離から与える。

$$
k_w = 0
$$

$$
\omega_w = \max\left(\frac{60\,\nu_w}{\beta_1 y_w^2},\; \omega_{\min}\right)
$$

ここで $y_w$ は最初のセル中心から壁までの距離、$\nu_w$ は壁面近傍の動粘性係数である。
`wall_isothermal` でも turbulence の考え方は同じで、温度境界は熱方程式側の責務とし、
SST 側は同じ $k_w$, $\omega_w$ を使う。

SST の理論上の ghost-cell closure は、壁値 $\phi_w$ を first-order の線形補間で与える形で書ける。
壁面がセル中心と ghost セルのちょうど中点にあるとき、

$$
\phi_w = \frac{\phi_c + \phi_g}{2}
$$

なので、Dirichlet 値を ghost へ戻す式は

$$
\phi_g = 2\phi_w - \phi_c
$$

である。したがって SST では

$$
k_g = -k_c
$$

$$
\omega_g = 2\omega_w - \omega_c
$$

となる。つまり、`k` は wall value が 0 なら ghost 側は負値になるが、これは「負の乱流エネルギーを物理的に意味づける」ものではなく、壁値 $k_w=0$ を中点補間で満たすための数値的な反射である。

壁からの距離が中点対称でない場合は、一般には距離比で係数が変わる。
現行コードは `kb` / `omegab` を壁値として持ち、ghost にはこの反射式を適用する形に更新した。

### 6.2 slip / axis / periodic

これらの境界は、SST scalar に対しては基本的に zero normal flux とする。

- `slip`: 法線方向勾配を持たせず、ghost へ copy する
- `axis`: 軸対称の対称条件として mirror / copy を使う
- `periodic`: 対向面の値をそのまま使う

このグループでは、壁のような強い散逸条件は入れない。

SST の立場では、`slip` / `axis` / `periodic` は壁面 Dirichlet ではなく、
\(\partial n\) に対する零流束条件として扱う。したがって ghost 値は通常

$$
\phi_g = \phi_c
$$

の copy で十分である。

### 6.3 inlet / outlet

流入は自由流乱流量を直接与える。実装上は boundary config から `k` と `omega` を受け取り、
必要に応じて turbulence intensity と length scale から換算する。

$$
k_\infty = \frac{3}{2} (U_\infty I)^2,
\qquad
\omega_\infty = \frac{\sqrt{k_\infty}}{C_\mu^{1/4} L_t}
$$

流出は 1 次外挿または zero-gradient を基本とし、逆流時のみ安定化のための clamp を入れる。

### 6.4 実装メモ

forge では `wall_dist` を既に持っているので、壁条件と source の両方で再利用する。
SST の壁条件は、単なる ghost copy ではなく、壁距離に基づく `\omega_w` の決定が必要になる。
実装上はまず boundary-plane 側で `kb`, `omegab` を持ち、そこから ghost セルへ

$$
k_g = 2k_w - k_c,
\qquad
\omega_g = 2\omega_w - \omega_c
$$

を与えるか、あるいは wall カーネルの中で ghost 値を直接埋める。
現行のコードは `wall` / `wall_isothermal` の boundary value を `kb` / `omegab` として保持し、
wall kernel 側で ghost へ反射して入れる。
等距離でない場合は、壁からセル中心までの距離を $d_c$、壁から ghost 中心までの距離を $d_g$
として

$$
\phi_g = \frac{(d_c + d_g)\phi_w - d_g\phi_c}{d_c}
$$

を使う。
したがって、壁条件は source よりも先に正しく閉じる必要がある。

## 7. 2D / 3D / 軸対称

2D と 3D は、既存ソルバと同様に同一の離散化カーネルを共有する。
2D は単に $u_z = 0$ と z 方向幾何が退化した特別な 3D として扱う。

軸対称では事情が異なる。既存の軸対称実装は B 流儀

$$
V = \bar{r} A_{planar},
\qquad
\mathbf{S}_f = \bar{r}_f \, dl_f \, \hat{\mathbf{n}}
$$

を用いている。advection-only 段階では既存の体積・面積重みで scalar の対流を
通せるが、SST 特有の diffusion / source には幾何項の検討が必要である。
そのため軸対称 SST の完成は別 plan とする。

## 8. 初期実装の範囲外

次は本章の対象外とする。

- 壁関数
- 遷移モデル
- 圧縮性補正の詳細オプション化
- 陰解法 Jacobian への SST 連成
- 軸対称 SST の幾何 source 詳細

それぞれ別 plan または後続フェーズで扱う。
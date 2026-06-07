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
詳細は子 plan
[`architecture-axisym-sst.md`](../../.github/plans/architecture-axisym-sst.md)
で扱い、以下に理論的な要点を整理する。

### 7.1 対流・拡散には幾何 source が残らない

$k$, $\omega$ の対流・拡散は円筒座標でも発散形を保つ。
$\phi \in \{k, \omega\}$ について

$$
\frac{1}{r}\partial_r\!\big(r\,\rho u_r \phi\big) + \partial_z(\rho u_z \phi)
=
\frac{1}{r}\partial_r\!\big(r\,\Gamma_\phi\,\partial_r \phi\big)
+ \partial_z(\Gamma_\phi\,\partial_z \phi)
+ S_\phi
$$

であり、両辺に $r$ を掛けて planar セル $A$ で積分すると、対流項・拡散項は
ともに $\partial_r(r\,\cdot) + \partial_z(r\,\cdot)$ の純粋な発散になる。
運動量方程式の $-\partial_r p$ のように「発散形に直せない残差」が無いため、
$r$ 重み付き面積 $\mathbf{S}_f = \bar r_f\,dl_f\,\hat{\mathbf n}$ を使う限り、
**対流・拡散からは軸対称特有の幾何 source は発生しない**。
拡散の $\theta\theta$ 寄与
$\frac{1}{r}\partial_r(r\,\Gamma_\phi\,\partial_r\phi)$ も $r$ 重み面積に
畳み込まれており、別途のソース項は不要である（子 plan で単体検証する）。

### 7.2 生産項のフープひずみ（幾何 source の本体）

幾何 source の本質は **ひずみ速度テンソルの周方向成分** にある。
軸対称流れ ($u_\theta = 0$, $\partial_\theta = 0$) でも、円筒座標のひずみ速度は
平面に現れない対角成分

$$
S_{\theta\theta} = \frac{u_r}{r} = \frac{U_y}{r}
$$

を持つ（フープひずみ, hoop strain）。SST の生産項は
$P_k = \mu_t S^2$, $P_\omega = \alpha\rho S^2$ で
$S = \sqrt{2\,S_{ij}S_{ij}}$ を用いるため、正しい大きさは

$$
S^2 = 2\!\left(S_{rr}^2 + S_{zz}^2 + S_{\theta\theta}^2\right)
+ \left(\partial_z u_r + \partial_r u_z\right)^2,
\qquad
S_{\theta\theta} = \frac{u_r}{r}
$$

となる。planar 速度勾配 ($\partial_x U_x \dots \partial_z U_z$) だけから組んだ
$S^2$ には $2\,S_{\theta\theta}^2 = 2\,(u_r/r)^2$ が欠落しており、軸近傍や急拡大部で
生産が過小評価される。同じ $S$ は渦粘性の制限子
$\mu_t = \rho a_1 k / \max(a_1\omega,\,S F_2)$ にも入るため、両方を一貫して
補正する必要がある。

### 7.3 圧縮性（dilatation）補正

超音速ノズルでは膨張・圧縮による発散 $\nabla\!\cdot\!\mathbf u$ が大きい。
軸対称の完全な発散は

$$
\nabla\!\cdot\!\mathbf u
= \partial_z u_z + \frac{1}{r}\partial_r(r u_r)
= \partial_x U_x + \partial_y U_y + \frac{u_r}{r}
\quad(=\texttt{axisym\_divU})
$$

であり、planar の発散に $u_r/r$ を加えた量となる。

#### 圧縮性 Boussinesq と生産項の正確形

圧縮性では乱流応力の Boussinesq 近似は等方成分を分離した形が正確である。

$$
\tau_{ij} = 2\mu_t\left(S_{ij} - \tfrac13(\nabla\!\cdot\!\mathbf u)\,\delta_{ij}\right)
- \tfrac23\rho k\,\delta_{ij}
$$

生産項 $P_k = \tau_{ij}\,\partial_j u_i = \tau_{ij} S_{ij}$ を展開すると

$$
P_k = \underbrace{2\mu_t\!\left(S_{ij}S_{ij} - \tfrac13(\nabla\!\cdot\!\mathbf u)^2\right)}_{\text{(A) deviatoric}}
\;\underbrace{-\;\tfrac23\rho k\,(\nabla\!\cdot\!\mathbf u)}_{\text{(B) isotropic}}
$$

となる。非圧縮形 $P_k = \mu_t S^2 = 2\mu_t S_{ij}S_{ij}$ (config `dilatationCorrection: 0`)
は (A) の $-\tfrac13(\nabla\!\cdot\!\mathbf u)^2$ と (B) を欠き、正確形との差は

$$
P_k^{\rm 現行} - P_k^{\rm 正確}
= \tfrac23\mu_t(\nabla\!\cdot\!\mathbf u)^2 + \tfrac23\rho k\,(\nabla\!\cdot\!\mathbf u).
$$

- **(A) deviatoric トレース除去**: 等方膨張・圧縮そのものはせん断生産に数えない。
  $S^2 \to S^2 - \tfrac23(\nabla\!\cdot\!\mathbf u)^2$（$S^2=2S_{ij}S_{ij}$ に対して
  $2\cdot\tfrac13(\nabla\!\cdot\!\mathbf u)^2 = \tfrac23(\nabla\!\cdot\!\mathbf u)^2$ を引く）。
  膨張コアでの spurious な $k$ 生成を抑える。低リスク。
- **(B) 等方項**: 膨張 ($\nabla\!\cdot\!\mathbf u>0$) で $k$ のシンク、
  圧縮 ($\nabla\!\cdot\!\mathbf u<0$, 衝撃波) でソース。$k$ 方程式固有で、
  $\omega$ 生産には入れない。膨張で $P_k<0$ になり得るため $P_k \ge 0$ クリップを併用する。

これら (A)(B) はいずれも $\nabla\!\cdot\!\mathbf u$ (= `axisym_divU`) を用いる。
軸対称では $u_r/r$ を含む完全発散を使わなければ強圧縮場で過大/過小評価が残る。

#### 別物: dilatation-dissipation 系（採用しない）

Sarkar / Zeman / Wilcox の乱流マッハ数補正
$\varepsilon = \varepsilon_s(1 + \alpha_1 M_t^2)$, $M_t=\sqrt{2k}/a$ は
**散逸**を増やす別系統の補正で、自由せん断層向けにキャリブレートされている。
壁境界層・ノズル内部では逆効果になりやすく、標準 SST には組み込まない。
本実装では (A)(B) の生産項補正のみを対象とし、$M_t$ 散逸補正は採用しない。

### 7.4 既存の軸対称診断量の再利用

7.2・7.3 で必要な量は、NS 用に
[`axisymmetricSource_d.cu`](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)
が既に計算・保持している。

| 量 | 変数 | 内容 |
| --- | --- | --- |
| フープひずみ $u_r/r$ | `axisym_uy_over_r` | $S_{\theta\theta}$ そのもの |
| 完全発散 $\nabla\!\cdot\!\mathbf u$ | `axisym_divU` | $\partial_x U_x+\partial_y U_y+u_r/r$ |

したがって SST 側は新規物理量を計算せず、これらを kernel 引数で受け取り、
`isAxisymmetric` のとき $S^2$ と発散補正に畳み込むだけでよい。

## 8. 初期実装の範囲外

次は本章の対象外とする。

- 壁関数
- 遷移モデル
- 圧縮性補正の詳細オプション化
- 陰解法 Jacobian への SST 連成

軸対称 SST の幾何 source (§7.1–7.4) は子 plan
[`architecture-axisym-sst.md`](../../.github/plans/architecture-axisym-sst.md)
で実装・検証する。その他は別 plan または後続フェーズで扱う。
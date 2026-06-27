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

### 6.5 automatic (enhanced / y⁺ 非依存) wall treatment

§6.1 の `\omega_w = 60\nu/(\beta_1 y^2)`, `k_w = 0` は **第一セルが粘性低層 (y⁺≲1) に
ある**ことを前提とした wall-resolved 型である。第一セルがバッファ層 (y⁺≈5〜30) や
対数層 (y⁺≳30) に落ちると、この `\omega_w` は過大評価になり、壁せん断応力 $\tau_w$ を
分子勾配 $\mu\,\partial U/\partial n$ から評価する運動量境界も $\tau_w$ を過小評価する。

automatic wall treatment は、粘性低層・バッファ層・対数層を **1 つの定式で滑らかに繋ぎ**、
第一セルの y⁺ 位置に依存せず妥当な $\tau_w$, $\omega_w$ を与える (Menter, *automatic
near-wall treatment*)。forge では opt-in (`wallTreatmentSST = 1`) とし、既定 `0` は §6.1 の
wall-resolved 型を保つ。

#### (a) 摩擦速度 $u_\tau$ — Reichardt 普遍速度則の逆解き

粘性低層・バッファ・対数を 1 式で繋ぐ Reichardt 則を用いる。

$$
u^+ = \frac{1}{\kappa}\ln(1+\kappa y^+)
      + 7.8\left[1 - e^{-y^+/11} - \tfrac{y^+}{11}e^{-y^+/3}\right]
$$

ここで $u^+ = U_t/u_\tau$, $y^+ = u_\tau y/\nu$, $U_t$ は壁接線速度の大きさ
(no-slip 壁では $U_t = |\mathbf{U}_c - (\mathbf{U}_c\cdot\hat{\mathbf n})\hat{\mathbf n}|$)、
$y$ は第一セル中心-壁距離 (`wall_dist`)、$\kappa = 0.41$。
$U_t$, $y$, $\nu$ は既知なので、

$$
f(u_\tau) = \frac{U_t}{u_\tau} - u^+\!\big(y^+(u_\tau)\big) = 0
$$

を Newton 法数回 (初期値 $u_\tau^{(0)} = \sqrt{\nu U_t / y}$、粘性則) で解く。
$y^+\to 0$ では $u^+ = y^+$ (粘性解) に、大 $y^+$ では対数則に自動収束する。

#### (b) $\omega$ 壁面境界 — Menter automatic ブレンド

粘性低層の解析漸近と対数層の値を二乗和の平方根で滑らかに繋ぐ。

$$
\omega_{vis} = \frac{6\nu}{\beta_1 y^2},
\qquad
\omega_{log} = \frac{u_\tau}{\sqrt{\beta^*}\,\kappa\,y},
\qquad
\omega_w = \sqrt{\omega_{vis}^2 + \omega_{log}^2}
$$

粘性低層では $\omega_{vis}$、対数層では $\omega_{log}$ が支配する。§6.1 の係数 60 ではなく
解析漸近の係数 6 を用いる ($y^+\to0$ の厳密漸近に一致させるため)。$\omega_w$ は全 $y^+$ で
妥当 ($y^+\to0$ で $\omega_{vis}$ 支配=壁解像値) なので、ghost 反射に加えて
**wall-adjacent セル値そのものを $\omega_c = \omega_w$ にピン留め**する (conserved $\rho\omega$ も整合)。
これにより渦粘性 $\mu_t = \rho a_1 k/\max(a_1\omega, S F_2)$ が近壁で適切に上限され、$k$ の暴走
(後述 (b')) を断つ。

#### (b') $k$ 壁面境界 — zero-gradient ($P_k$ の wall-function 化とセット)

wall-function では $k$ を **zero-gradient (Neumann, $\partial k/\partial n = 0$)** とするのが標準
(OpenFOAM `kqRWallFunction`、ANSYS automatic の $k$ zero-flux)。第一セルがバッファ・対数層に
あるとき $k$ は有限の対数層平衡値 $k \approx u_\tau^2/\sqrt{\beta^*}$ を取るべきだからである。

**ただし $k$ zero-gradient は近壁の生産項 $P_k$ の wall-function 化 (下記 (d)) とセットで初めて
成立する。** $P_k$ を標準 SST のまま**解像速度勾配**から計算したまま $k$ を zero-gradient にすると、
粗メッシュでは誤った解像生産で近壁 $k$ が暴走し $\mu_t$ 過大 → $u_\tau$/$C_f$ 過大になる
(実測 `case/26`: $k$ zero-gradient 単独で y⁺≈30 が $C_f/C_{f,\text{corr}}$ 0.91→1.80, y⁺≈10 が
1.10→1.68 と悪化した)。(d) の $P_k$ 置換と (b) の $\omega$ ピン留めを同時に入れると、$k$ は
自律的に平衡値へ収束し暴走しない。したがって automatic モードでは
$$
k_w = k_c,\qquad k_g = k_c \quad(\partial k/\partial n = 0)
$$
とする (wall-resolved モード=既定 0 は $k_w = 0$, $k_g = -k_c$ を保つ)。

#### (c) 運動量壁せん断 — 有効壁粘性

第一セルが粗いと分子勾配は $\tau_w$ を過小評価するため、モデル化した壁せん断
$\boldsymbol{\tau}_w = \rho u_\tau^2\,\hat{\mathbf e}_t$ ($\hat{\mathbf e}_t$ は接線速度方向)
を運動量残差に直接課す。等価な有効壁粘性は

$$
\mu_{eff} = \frac{\rho u_\tau^2\, y}{U_t}
$$

であり、$y^+\to0$ では $u_\tau^2 = \nu U_t/y$ より $\mu_{eff}\to\mu$ (wall-resolved に一致)、
対数層では $\mu_{eff} > \mu$ (wall function 化) となり、全層で連続に繋がる。
no-slip 壁では壁せん断は仕事をしない ($\boldsymbol{\tau}_w\cdot\mathbf{U}_w = 0$) ため、
エネルギー方程式の壁せん断仕事項と熱流束は不変である。

#### (d) wall-adjacent セルの乱流生産 $P_k$ — wall-function 値への置換

(b') の $k$ zero-gradient を成立させる要。wall-adjacent セルでは標準 SST の解像勾配生産
$P_k = \mu_t S^2$ を、普遍速度則から評価した wall-function 生産に置換する。乱流応力 (= 全応力
− 粘性応力) と壁面速度勾配の積として
$$
P_k = (\tau_w - \tau_{vis})\,\frac{\partial U}{\partial y}
    = \rho\frac{u_\tau^4}{\nu}\,g\,(1-g),
\qquad
g \equiv \left.\frac{du^+}{dy^+}\right|_{y^+_1}
$$
ここで $\partial U/\partial y = (u_\tau^2/\nu)\,g$, $\tau_w = \rho u_\tau^2$, $\tau_{vis} = \mu\,\partial U/\partial y
= \tau_w\,g$ を用いた。$g$ は (a) の Reichardt 則の微分 $du^+/dy^+$ を第一セル $y^+_1$ で評価。

- **粘性低層** ($y^+\to0$): $g\to1$ より $P_k\to0$ — せん断は全て分子粘性が担い乱流生産ゼロ。
  標準 SST の wall-resolved 極限に一致する (細メッシュで既定 0 と同結果)。
- **対数層** ($y^+\gg1$): $g\to 1/(\kappa y^+)$ より $P_k\to \rho u_\tau^3/(\kappa y)$ — 標準的な log 則生産。

この $P_k$ と (b) の $\omega$ ピン留めにより、$k$ 方程式の生産・消滅平衡
$P_k = \beta^*\rho k\,\omega_w$ から $k\to u_\tau^2/\sqrt{\beta^*}$ (対数層平衡値) へ自律収束し、
$k$ を直接固定せずとも暴走しない。検証は plan / `case/26` を参照 (3 帯 y⁺≈0.35/10/30 で
$C_f/C_{f,\text{corr}}$ がいずれも 0.9–1.0 に収まる)。

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
[`architecture-axisym-sst.md`](../../plans/accepted/architecture-axisym-sst.md)
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

### 7.5 応力生産アノマリーと Kato–Launder 補正

標準 SST の生産項 $P_k = \mu_t S^2$ ($S=\sqrt{2S_{ij}S_{ij}}$) は **ひずみの大きさ**のみで
生産を見積もる。ところがひずみテンソルは

$$
S_{ij} = \underbrace{\tfrac13(\nabla\!\cdot\!\mathbf u)\delta_{ij}}_{\text{等方(膨張)}}
+ \underbrace{S'_{ij}}_{\text{偏差}},\qquad
S'_{ij} = \underbrace{\text{せん断成分}}_{\Omega\ \text{を伴う}}
+ \underbrace{\text{法線(伸長)異方性}}_{\text{非回転}}
$$

と分解され、偏差ひずみ $S'$ は**せん断**だけでなく**非回転の伸長(加速)ひずみ**も含む。
速度勾配の反対称部=渦度 $\Omega_{ij}=\tfrac12(\partial_j u_i-\partial_i u_j)$,
$\Omega=\sqrt{2\Omega_{ij}\Omega_{ij}}=|\boldsymbol\omega|$ は回転の不変尺度であり、
**単純せん断では $S\approx\Omega$、純伸長(よどみ点・加速)では $\Omega\to0$ かつ $S$ 大**となる。

§7.3 の `dilatationCorrection` は等方(トレース=$\nabla\!\cdot\!\mathbf u$)分のみを除くため、
**非回転の偏差伸長ひずみは残り**、強加速場（ノズル喉中心線・よどみ点）で
$P_k=\mu_t S'^2$ が偽の乱流を生む。これは strain-based SST の古典的欠陥
(**stagnation / round-jet anomaly**) であり、`dilatationCorrection` では消えない別系統の問題。

**Kato–Launder 修正** (Kato & Launder 1993) は生産の片方の $S$ を渦度 $\Omega$ に置換する:

$$
P_k = \mu_t\,S\,\Omega \qquad (\text{標準は } \mu_t S^2)
$$

- せん断層 ($S\approx\Omega$): $P_k\approx\mu_t S^2$ で**従来と一致**(壁 BL の本物の乱流は不変)。
- 非回転加速 ($\Omega\to0$): $P_k\to0$ で**偽生産を除去**。

不変量 $S,\Omega$ を使うため**座標系に依存しない**(「対角項除外」のような非不変操作は不可:
45° 回転でせん断と伸長が入れ替わるため)。

#### case 29 の「軸中心 k 過大」— 機構特定済み: **二層構造** (経緯の誤診は訂正済み)

詳細は [`architecture-axisym-axis-singularity.md`](../../plans/accepted/architecture-axisym-axis-singularity.md) /
[case 29 検証](../../case/29.bell_vs_conical/README.md)。結論:

**層 (a) — 実用格子 (軸セル ≥ 0.25mm) の軽度スパイク = 生産駆動。Kato–Launder で解決。**
平均流は滑らか。軸第1セル $k$=17 vs 核 4.5 は、大きな軸方向ひずみ + フープ + 生産リミタ正帰還 +
軸無フラックス滞留による。`katoLaunder: 1` で $k$ 17→1.9、壁 BL・推力不変 — 3D の平坦
プロファイル (第1セル/核 ~0.85) と整合する妥当な対処。

**層 (b) — 軸を ≤10µm に細分すると平均流の数値不安定 (未根治)。**
2 面クラスタ格子 (壁解像維持 + 軸細分) で **$U_y$ が近軸チェッカーボード (偶奇) 振動**
(物理 ~0.5 m/s のところ ±7 m/s 交番)。勾配計算は場の FD と一致 (正しい) — 場そのものが振動し、
$(\partial_r U_y)^2 \sim 10^{11}$ が $S^2$ 経由で $k$ を爆発させる。軸細分で悪化 (k 収束せず /
rms_roK 1e11 爆発)。機構: 半径方向は近軸で深い低 Mach ($u_r\to0$) → 風上流束の圧力-速度結合が
消え偶奇モードが減衰しない (古典的チェッカーボード)。3D は軸が内部点でこの構造を持たない。
`lowMachPrecond` 1/2 とも失敗 (超音速コア/細軸で発散)。

**実用指針**: 軸対称 RANS ノズルは `katoLaunder: 1` を使い、**軸方向の過細分を避ける**
(軸セル ≥ ~0.1mm; 壁解像は壁側クラスタで)。

> 経緯メモ (誤診の訂正): ①「planar 面積 vs $r$ 重み体積の不整合」は誤り (flux 面積・体積は両方
> $r$ 重み済, [`variables.cpp`](../../solver_density_cuda/variables.cpp))。②「勾配の planar 面積が
> bug」も誤り (勾配は planar 作用素で正しい)。③ uniform 格子のグリッド試験は $k$ 非収束で無効
> だった。最終特定は 2 面クラスタ格子での $U_y$ チェッカーボード観測による。

#### 適用範囲・既定

- config `katoLaunder` (turbulence): `0`=標準 $\mu_t S^2$ (既定, 既存挙動とビット一致),
  `1`=Kato–Launder $\mu_t S\Omega$。
- $\omega$ 生産 $P_\omega=\alpha\rho S^2$ にも同じ置換 ($\alpha\rho S\Omega$) を適用し、$\mu_t$ の
  時間スケール整合を保つ。
- `dilatationCorrection` と**直交・併用可** (等方分は除いた上で $S\to\sqrt{S^2_{\rm corr}}$ と
  $\Omega$ の積を取る)。

## 8. DES / DDES (Detached-Eddy Simulation)

剥離支配の非定常流 (SBLI・ピントルバルブ・後退段) を RANS より高精度・wall-resolved LES より
低コストで解くため、SST に **DDES (Delayed DES, Spalart et al. 2006)** の length scale 修正を加える。
付着境界層は RANS のまま保護し、剥離せん断層・自由せん断層でのみ LES へ切り替える。
実装は [`implementation.md` §3.8](implementation.md) と plan
[`turbulence-iddes-sst.md`](../../plans/active/turbulence-iddes-sst.md) を参照。

### 8.1 RANS length scale と k 消滅項の整合

SST の k 方程式の消滅項は
$$
D_k = \beta^{*}\rho k\omega = \frac{\rho k^{3/2}}{l_\mathrm{RANS}}
$$
である。この恒等式から RANS length scale は
$$
l_\mathrm{RANS} = \frac{\sqrt{k}}{\beta^{*}\omega}\qquad(\beta^{*}=0.09)
$$
**でなければならない**。$\rho k^{3/2}/l_\mathrm{RANS} = \rho k^{3/2}\cdot\beta^{*}\omega/\sqrt{k} = \beta^{*}\rho k\omega$
と一致する。これは Strelets (2001) / Gritskevich (2012) の SST-DES で用いられる定義であり、
DES は **この $l_\mathrm{RANS}$ を $l_\mathrm{DDES}$ に置換するだけ**で、$f_d=0$ の完全シールド域では
$l_\mathrm{DDES}=l_\mathrm{RANS}$ となり $D_k$ が標準 SST に厳密に戻る。

> **注意 (実装で判明した落とし穴)**: $l_\mathrm{RANS}=\sqrt{k}/(\beta^{*1/4}\omega)$ と書くと上の恒等式を
> 満たさず、$D_k=\beta^{*1/4}\rho k\omega$ となって標準 SST より係数が $\beta^{*-3/4}\approx 6.1$ 倍
> 大きい消滅項になる。この場合 $f_d=0$ の付着 BL でも mode1 が標準 SST に縮退せず、
> **モデル応力枯渇 (Modeled Stress Depletion)** で BL の $k$ が崩壊する (検証で `roK` が ~100% 変化)。
> 必ず $\beta^{*}$ (=0.09) を使うこと。

### 8.2 LES グリッドスケール $\Delta_\mathrm{max}$

$$
\Delta = \Delta_\mathrm{max} = \max_{\text{隣接面}} |\mathbf{cc}_{ic} - \mathbf{cc}_{jc}|
$$

隣接 CV 重心間距離の最大値。等方メッシュなら $V^{1/3}$ で代用できるが、**境界層の高アスペクト比
セル**では $V^{1/3}$ が接線スケールを大幅に過小評価し、シールドがあっても DES リミッタが
BL 内で誤発火しうる。$\Delta_\mathrm{max}$ は Spalart 系で標準。幾何セットアップ時に実行時の
`volume`/`ccx` (= CV 重心) から 1 回だけ計算するので、cell / node (median-dual) 双方で自動的に
双対 CV の $\Delta$ になる。

### 8.3 DDES シールド関数 (Spalart et al. 2006)

$$
r_d = \frac{\nu_t + \nu}{\kappa^2 d^2\,\sqrt{\partial U_i/\partial x_j\,\partial U_i/\partial x_j}}\qquad(\kappa=0.41)
$$
$$
f_d = 1 - \tanh\!\bigl([8\,r_d]^3\bigr)
$$

| $r_d$ | 物理的意味 | $f_d$ | 動作 |
|-------|-----------|-------|------|
| $\gtrsim 1$ | BL 内部 (渦粘性大 / 近壁 $d^2$ 小) | $\approx 0$ | 純 RANS (保護) |
| $\to 0$ | 剥離域・自由せん断層 | $\to 1$ | DES limiter 有効 |

分子の $(\nu_t+\nu)$ により、対数層では $\nu_t\!\sim\!\kappa u_\tau y$, $S\!\sim\!u_\tau/(\kappa y)$ から
$r_d\!\to\!1$ となり付着 BL が保護される。粘性低層では $\nu$ 項により $r_d\!\sim\!1/(\kappa^2 y^{+2})$ で
やはり保護される。

### 8.4 ハイブリッド length scale

$$
l_\mathrm{DDES} = l_\mathrm{RANS} - f_d\cdot\max\!\bigl(0,\; l_\mathrm{RANS} - C_\mathrm{DES}\Delta\bigr)
$$
$$
D_k^\mathrm{DDES} = \frac{\rho k^{3/2}}{l_\mathrm{DDES}},\qquad C_\mathrm{DES}=0.78\ (\text{初版 }k\text{-}\omega\text{ 固定})
$$

リミッタは **2 重に保護**される: ① シールド $f_d\approx 0$、または ② $l_\mathrm{RANS}<C_\mathrm{DES}\Delta$
(グリッドが BL 長スケールに対して粗い near-wall) で $\max(0,\cdot)=0$ となり、どちらでも
$l_\mathrm{DDES}=l_\mathrm{RANS}$ に戻る。これにより付着 BL は $f_d$ が多少立っても撹乱されない。
$\omega$ 方程式・渦粘性式は **変更しない** (標準 SST-DES は $k$ 消滅項のみ修正)。

$C_\mathrm{DES}$ は本来 $F_1$ ブレンド $C_\mathrm{DES}=F_1 C_{\mathrm{DES},k\omega}+(1-F_1)C_{\mathrm{DES},k\varepsilon}$
($0.78/0.61$) だが、渦粘性カーネルが $F_1$ を持たないため初版は $C_{\mathrm{DES},k\omega}=0.78$ 固定とする。

### 8.5 float32 ガード

forge は float32。$r_d$ は $d^2$ で効くため近壁・よどみ・一様流で飽和しやすい。
勾配大きさ・分母に floor を入れて 0 割を防ぎ ($\nabla u\to 0$ の一様流で $r_d\to$ 大 $\to f_d\to 0$ =
RANS が正しい極限)、$r_d$ を $[0,10]$ に clamp して $f_d$ の 0/1 張り付きを抑える。
診断のため `r_d`・`f_d`・`l_des`・`delta_les` を HDF5 出力する (一様な $f_d$ 張り付きは NaN より
発見しにくい DES 不全の兆候)。

### 8.6 IDDES (Phase 2・未実装)

$f_B$, $f_e$, $f_{dt}$ を加えた $l_\mathrm{IDDES}$ で WMLES モードへ自動切替する Shur et al. (2008) の
拡張は Phase 2。WMLES モードの定量評価には流入乱流生成が前提のため別 plan とする。

---

## 9. 初期実装の範囲外

次は本章の対象外とする。

- 壁関数
- 遷移モデル
- 圧縮性補正の詳細オプション化
- 陰解法 Jacobian への SST 連成

軸対称 SST の幾何 source (§7.1–7.4) は子 plan
[`architecture-axisym-sst.md`](../../plans/accepted/architecture-axisym-sst.md)
で実装・検証する。その他は別 plan または後続フェーズで扱う。
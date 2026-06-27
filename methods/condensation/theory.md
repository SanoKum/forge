# 非平衡凝縮 (4 モーメント方程式) — 理論

極超音速ノズルの強膨張で試験ガス (窒素 N$_2$ 等) が**非平衡凝縮**し、潜熱放出で静圧・静温が
上昇する現象を、気相 Euler/NS に**液相のモーメント輸送方程式**を結合して記述する。

主参照: L. Lin, W. Cheng, X. Luo, F. Qin, *On nitrogen condensation in hypersonic nozzle flows*,
Shock Waves 24:179–189 (2014)。物理係数・導出の詳細は
[`papers/.../_summary.md`](../../papers/on%20nitrogen%20condensation%20in%20hypersonic%20nozzle%20flows_summary.md)
にまとめてある。実装対応は [implementation.md](implementation.md) を参照。

---

## 1. 支配方程式 — Euler/NS + method of moments

液滴雲はサイズ分布関数 $f(r;x,t)$ の有限個の**モーメント** $Q_n=\int_0^\infty r^n f\,dr$ で表す
(Hill 1966)。凝縮種ごとに次の保存ベクトルを気相に追加する:

$$
U = (\rho,\ \rho u,\ \rho v,\ \rho w,\ \rho E,\ \underbrace{\rho g,\ \rho Q_2,\ \rho Q_1,\ \rho Q_0}_{\text{凝縮種ごと 4 本}})^T
$$

- $\rho Q_0$ : 単位質量あたり液滴数密度
- $\rho Q_1,\ \rho Q_2$ : サイズ分布のモーメント
- $\rho g$ : 液相質量分率 $g = \tfrac{4}{3}\pi\rho_l Q_3$ (実質 $Q_3$ モーメント、$\rho_l$ は液密度)

各モーメント方程式は**移流項 + 凝縮ソース $S$ (核生成 + 成長)** を持つ:

$$
\frac{\partial(\rho Q_n)}{\partial t} + \nabla\!\cdot(\rho Q_n \mathbf u) = S_{Q_n},\qquad
\frac{\partial(\rho g)}{\partial t} + \nabla\!\cdot(\rho g\, \mathbf u) = S_g
$$

気相の質量・運動量・エネルギー式のソースは**ゼロ** (混合気体全体は保存)。気相と液相の結合は
**熱力学 (状態方程式・温度関係) 経由**で、潜熱が静温を上げる。モーメントは気相速度 $\mathbf u$ で
運ばれる受動スカラー (拡散なし) として移流される。

ソース項は核生成項 ($T$ 依存) と成長項 ($T_d$ 依存) の加算分離:

$$
S_{Q_n} = \underbrace{J\, r_*^{\,n}}_{\text{核生成}} + \underbrace{n\,\rho Q_{n-1}\frac{dr}{dt}}_{\text{成長}},\qquad
S_g = \frac{4}{3}\pi\rho_l\Big( J\, r_*^{\,3} + 3\rho Q_2\frac{dr}{dt}\Big)
$$

---

## 2. 核生成モデル (nucleation) — **気相温度 $T$ で評価**

古典核生成理論 (CNT, Becker–Döring, capillarity 近似):

$$
J_{CNT} = K_{CNT}\exp\!\left(-\frac{\Delta G^*}{k_B T}\right),\quad
K_{CNT}=\Big(\frac{2\sigma}{\pi M\upsilon_l}\Big)\Big(\frac{p_v}{k_B T}\Big)^2,\quad
\Delta G^*=\frac{4}{3}\pi r_*^2\sigma
$$

臨界半径 (Kelvin–Thomson):

$$
r_* = \frac{2\sigma}{\rho_l R T\,\ln(p_v/p_{sat})}
$$

過飽和 $S=p_v/p_{sat}(T)$。**核生成は気相温度 $T$ で評価**する (臨界核は sub-nm で自己加熱なしとみなす)。

### 種ごとの補正 (切替可能)

| モデル enum | 補正 | 用途 |
| --- | --- | --- |
| `CNT` | 補正なし | 基準 |
| `CNT_Iland` | $J=J_{CNT}\exp(A+B/T)$, $A=-55,\ B=4270$ K | **N2** (Iland et al. 2009) |
| `CNT_Kantrowitz` | $J_{noniso}=J_{iso}/(1+b\,q^2)$ 系の非等温補正 | **H2O** (生成中クラスタ自己加熱) |

- **Iland 補正** (N2): 膨張チェンバー実験から決めた経験補正。本論文の圧力域 (2.8–6.0 kPa) で検証済み。
- **Kantrowitz 非等温補正** (H2O): 生成中の臨界核が放出潜熱を捨てきれず自己加熱し $J$ が下がる効果を
  **核生成項そのもの**に織り込む。N2 では使わず (自己加熱は成長側 $T_d$ で扱う、下記)。

---

## 3. 成長モデル (droplet growth) — **液滴温度 $T_d$ で評価**

### 修正 Gyarmathy / Goodheart (N2)

純窒素蒸気で広い Knudsen 数域に適用可能 (Goodheart 2004):

$$
\frac{dr}{dt} = \frac{1}{\rho_l}\,\frac{f_{F\text{-}S}(Kn)\,k\,R T^2}{L^2}\,\ln\!\left(\frac{p}{p_{sat}(T_d)}\right),\qquad
f_{F\text{-}S}(Kn)=\frac{1+2Kn}{r(1+3.42Kn+5.32Kn^2)}\Big(1-\frac{r_*}{r}\Big)
$$

$k$: 熱伝導率、$L$: 潜熱、$Kn=\lambda/2r$ (核生成の nm から成長後の µm までカバー)。
**成長は液滴温度 $T_d$ で評価**する (有限サイズ液滴は蓄積潜熱で $T_d>T_v$ になる)。

### Hertz–Knudsen (H2O、切替候補)

質量フラックスベースの成長則。H2O 検証 (Wyslouzil) 用に同じ枠組みで切替可能とする。

---

## 4. 液滴温度 $T_d$ — Hill のエネルギーバランス (純蒸気簡約)

$T_d$ は液滴表面のエネルギーバランス (Hill モデル, Johnson 2001 の比較で採用) から反復で解く。
純窒素では不活性キャリアガス寄与 $\beta_c$ 項を落とし ($\zeta=1.0$):

$$
K_v\!\left(1-\frac{\beta_d T_d}{\beta T_v}\right) + \left(1-\frac{\beta_d}{\beta}\right)(L_v-1)\frac{\gamma}{\gamma-1} = 0
$$

$$
K_v=\frac{1}{2}\frac{\gamma+1}{\gamma-1},\quad
L_v=\frac{(\gamma-1)L}{\gamma R T_v},\quad
\beta=\frac{p_v}{\sqrt{2\pi R T_v}},\quad
\beta_d=\frac{p_d}{\sqrt{2\pi R T_d}}
$$

液滴表面の蒸気圧は **Kelvin 効果込み**:

$$
p_d = p_\infty(T_d)\,\exp\!\left\{\frac{2\sigma}{\rho_l R T_d\, r}\right\}
$$

$\beta_d$ が $p_d/\sqrt{T_d}$ を通じて $T_d$ に依存し、$p_d$ も $T_d,r$ に依存するため、$T_v,p_v,r$ を所与に
$T_d$ を反復で求める。$T_d$ は成長則 (式 3) と Kelvin 式にのみ入り、**核生成には直接入らない** (核生成は
バルク $T$ 上昇を介した間接ループでのみ抑制される)。

### 実装した準定常 $T_d$ 閉包 (carrier + Hertz–Knudsen, `condTwoTemp=1`)

上の Hill 式は純蒸気 (N2) 向けの整理である。希薄水/N2 (carrier) では、**準定常の液滴エネルギーバランス**
「凝縮による潜熱解放 = 液滴から気相 (キャリア) への熱伝導」を直接解く形で実装した。Hertz–Knudsen
(質量律速) 成長に対し、質量フラックス $j$ と熱フラックスを連立する:

$$
\underbrace{L\,j(T_d)}_{\text{潜熱解放}} = \underbrace{h\,(T_d-T_g)}_{\text{気相への熱伝導}},\qquad
j(T_d)=\frac{\alpha\,[\,p_v-p_d(T_d)\,]}{\sqrt{2\pi R T_g}},\qquad
h=\frac{\lambda_g}{r\,(1+3.18\,Kn)}
$$

液滴表面の平衡蒸気圧は **液滴温度 $T_d$ で** Kelvin 補正: $p_d(T_d)=p_{sat}(T_d)\exp\!\{2\sigma/(\rho_l R T_d r)\}$。
$\lambda_g$ はキャリア (N2) の熱伝導率、$Kn=\lambda/2r$ (平均自由行程は全圧基準)。この 1 変数非線形方程式
$F(T_d)=L\,j(T_d)-h(T_d-T_g)=0$ を **Newton 法** ($\partial p_{sat}/\partial T_d$ と Kelvin 項の微分を含む解析勾配,
$T_d\ge T_g$ にクランプ) で解き、得た $T_d$ で成長率 $\dot r=j(T_d)/\rho_l$ を評価する。

- **効果**: 自己加熱で $T_d>T_g$ となり $p_d(T_d)>p_d(T_g)$、実効過飽和 $p_v-p_d$ が下がるため**成長が遅くなる**
  (一温度 $T_d=T_g$ は過冷却を過大評価し成長を過大評価していた)。$Kn$ が大きい (小液滴) ほど $h$ が小さく
  $T_d$ の上振れが大きい。
- **適用範囲**: 本フラグは **Hertz–Knudsen 経路 (`condGrowthModel=0`) のみ**に効かせる。Gyarmathy
  (`condGrowthModel=1`, 熱伝導律速) は元来「液滴が $T_s$ 付近まで自己加熱した極限」を仮定し過冷却
  $(T_s-T_g)$ を駆動力に使う式なので、$T_d$ を二重計上しないため非適用 (一温度のまま)。
- **簡約**: $T_d$ は**成長キネティクスのみ**に使う。二相 EOS のバルクエネルギーは液滴を $T_g$ のままとする
  (液相顕熱差 $\sim g\,c_l\,\Delta T_d \sim 1\%$ オーダーで潜熱 $\sim$ 数% に対し小、既知の近似)。
- **既定 off** (`condTwoTemp=0`, 一温度 $T_d=T_g$) でビット不変。

---

## 5. 二相 EOS と圧力 — 気相分圧近似

窒素を熱量的完全気体として扱う (実在気体効果 $<4\%$)。**液滴体積を無視**し、蒸気を理想気体と
みなすと、圧力に効くのは気相分圧だけ:

$$
\boxed{\ p = (\rho-\rho g)\,R\,T\ } \qquad (R: \text{蒸気の気体定数})
$$

混合の内部エネルギー (潜熱込み) から $T$ を逆算:

$$
\rho e = (\rho-\rho g)\,e_v(T) + \rho g\,e_l(T),\qquad \rho e=\rho E-\tfrac12\rho|\mathbf u|^2
$$

$e_v-e_l = L - RT$ (潜熱)。混合体積比熱 $C_v=(\rho-\rho g)c_{v,v}+\rho g\,c_l$。

### 圧力微分 (一般 EOS Roe の $(\chi,\kappa)$ + 凝縮新項 $\xi_g$)

$$
\kappa\equiv\left.\frac{\partial p}{\partial(\rho e)}\right|_{\rho,\rho g}=\frac{(\rho-\rho g)R}{C_v}
\quad(\gamma-1\ \text{に相当}),\qquad
\chi\equiv\left.\frac{\partial p}{\partial\rho}\right|_{\rho e,\rho g}=RT-\kappa\,e_v
$$

$$
\boxed{\ \xi_g\equiv\left.\frac{\partial p}{\partial(\rho g)}\right|_{\rho,\rho e}=-RT+\kappa(e_v-e_l)=-RT+\kappa(L-RT)\ }
$$

$\xi_g$ は二相化で**新規に出る圧力微分**で、「$\rho g$↑ = 蒸気→液化 = 潜熱放出 → $T$↑ → $p$↑」を表す
(N2 70 K 概算で $\xi_g\approx+5\times10^4>0$、凝縮で昇圧)。保存量での $\partial p/\partial U$ と flux
Jacobian への入り方、Roe/SLAU への影響は [implementation.md](implementation.md) の対流フラックス節を参照。

### 保存量での $\partial p/\partial U$ と semi-perfect (TP) 一般化

内部エネルギー密度 $\tilde\epsilon=\rho e=\rho E-\tfrac12|\mathbf m|^2/\rho$ ($\mathbf m=\rho\mathbf u$) の
全微分 $dp=\chi\,d\rho+\kappa\,d\tilde\epsilon$ と連鎖則
($\partial\tilde\epsilon/\partial\rho|_{\mathbf m,\rho E}=e_k\equiv\tfrac12|\mathbf u|^2$,
$\partial\tilde\epsilon/\partial m_k=-u_k$, $\partial\tilde\epsilon/\partial(\rho E)=1$) より、保存量
$q=(\rho,\rho u,\rho v,\rho w,\rho E)$ に対して

$$
\frac{\partial p}{\partial q_1}=\chi+\kappa\,e_k,\quad
\frac{\partial p}{\partial q_{2,3,4}}=-\kappa(u,v,w),\quad
\frac{\partial p}{\partial q_5}=\kappa
$$

凝縮種は $\partial p/\partial(\rho g)=\xi_g$ の列が加わる。

**semi-perfect (thermally-perfect) gas**: 比熱が $T$ 依存で $\epsilon(T)=\int_0^T C_v\,dT'\ne C_v(T)\,T$
のとき、$T,\rho,\tilde\epsilon$ の全微分
($(\partial T/\partial\rho)|_{\tilde\epsilon}=-\epsilon/(\rho C_v)$,
$(\partial T/\partial\tilde\epsilon)|_\rho=1/(\rho C_v)$) を用いて

$$
\chi=R\Big(T-\frac{\epsilon(T)}{C_v(T)}\Big)=\frac{p}{\rho}-(\gamma(T)-1)\epsilon(T),\qquad
\kappa=\gamma(T)-1
$$

となり、**$\chi\ne0$** となる (calorically-perfect は $\epsilon=C_v T$ の特殊形で $\chi=0$ に縮約)。
本書の $\chi=RT-\kappa e_v$ は $e_v$ を**積分内部エネルギー** $\epsilon_v(T)$ と取れば
($\kappa e_v=(R/C_v)\epsilon_v=(\gamma-1)\epsilon_v$) この一般形と一致する。forge は既に多成分 TP の陰解法
Jacobian で per-cell $\kappa=\gamma[ic]-1$ を実装済み ([thermophysics plan](../../plans/accepted/thermophysics-multicomponent-tpgas.md))
であり、凝縮 Phase 2 は同枠で $\kappa\to(\rho-\rho g)R/C_v$ (混合) と $(\rho g)$ 列の $\xi_g$ を足す最小拡張となる。

---

## 6. 数値解法 — fractional-step

支配方程式を「ソース無しの均質部 (対流)」と「凝縮ソース部」に**演算子分割** (fractional-step) する。
均質部は forge 既存の有限体積対流 (SLAU/Roe) + 時間積分で解き、凝縮ソース ($J,\ dr/dt$) は別途
積分する。核生成率 $J\propto\exp(-\Delta G^*/k_BT)$、$dr/dt\propto\ln(p/p_{sat})$ は過飽和に指数的・
対数的で極めて stiff なため、ソース積分は **point-implicit** (源項ヤコビアン $\partial S/\partial U$ を
対角に組む) とする。

### 精度・無次元化

数密度 $\rho Q_0\sim10^{18}$–$10^{22}$ /m$^3$、$r\sim$ nm のため $Q_0$ と $Q_3$ は $\sim$27 桁跨ぐ。
モーメントを基準量で無次元化して O(1) に抑える:

$$
\mu_n = \frac{Q_n}{N_{ref}\,r_{ref}^{\,n}}\quad(\text{例 } N_{ref}=10^{18}\,/\text{kg},\ r_{ref}=1\,\text{nm})
$$

平均半径 $r=(\mu_3/\mu_0)^{1/3} r_{ref}$ 等のモーメント間比も無次元量どうしで安定する。**まず float
で進め**、無次元化で桁を抑える。float で桁落ちが顕在化したらモーメントの double 化 (混合精度) を
フォールバックとして検討する (詳細は implementation.md)。

---

## 7. 検証

1. **dry 一致**: 凝縮ソースを切った状態 ($g\equiv0$) で case/34 Arthur ノズル run_0003 と一致 (回帰防止)。
2. **N2 凝縮 ON**: case/34 で貯気を振り、中心線静圧が dry 等エントロピー線より**上振れ**することを示す
   (論文 Fig.11)。
3. **H2O 凝縮**: case/16.nozzle_wys + `papers/condensation/wyslouzil2000.pdf` で水蒸気凝縮を検証し、
   N2 と H2O を同じ枠組みで切替できることを実証する。

---

## 8. N2 物性フィット (Lin 2014 Appendix 1)

実装で使う窒素の物性相関。臨界 $T_c=126.192$ K (表面張力は $T_c=126$ K)、臨界密度
$\rho_c=313.3$ kg/m³、三重点 $T_t=63.15$ K、$R_{N_2}=296.8$ J/(kg·K)。

### 相の方針 — 過冷却液 (supercooled liquid)、固体ではない

**Lin 2014 のモデルは液相 (過冷却液)**。非平衡凝縮は三重点以下でも準安定な**過冷却液滴**を作る
(均質核生成は液滴を生成、結晶化は膨張時間に対して遅い) ため、本実装も**全域で液相フィットを使う**
(固相昇華へは切替えない)。論文本文も "conservation of the **liquid** phase", "**liquid** mass
fraction $g$", "condensation of **liquid droplets**" とし、$r_*=2\sigma/(\rho_l RT\ln(p_v/p_{sat}))$
は**液密度 $\rho_l$** を使う。Appendix 1 が固相フィット ($p_{sat}^s,\rho_s,L_s,\sigma_s$) も併載するのは
図 (Fig.9,10) の比較用であり、**シミュレーションでは使わない**。

> **低温外挿の注意 (実装ガード)**: 液相フィット (Jacobsen $p_{sat}$、潜熱多項式) は **~40K 以下に外挿
> すると破綻**する ($p_{sat}$ が ~33K で非単調、潜熱が崩落)。有効域は概ね **45–126K**。凝縮 ON では
> 潜熱放出で活性域が ~40K 以上に留まる (論文 Ma6 ノズルで出口 40→56K) が、過渡/dry セルが低温に落ちても
> 破綻しないよう、**物性評価温度を $[45\text{K}, T_c)$ にクランプ**する (より低温=より低 $p_{sat}$ の
> 安全側=過剰核生成を起こさない)。既知の近似。実装 `COND_T_PROP_FLOOR`
> ([condensationProperties_d.cuh](../../solver_density_cuda/cuda_forge/condensationProperties_d.cuh))。

> **気相 thermo の制約**: 凝縮ケースの気相は **`thermalMethod 0` (熱量的完全気体, $c_p/\gamma$ 一定)**
> を使うこと。NASA-9/CEA は ~200K 未満で無効 (出口 ~27K) で、低温気相物性は CEA に無い。論文も
> calorically perfect gas。

以下のフィット式は参照用に液・固双方を併記するが、**既定実装は液 (過冷却) のみ**を使う。

### 飽和蒸気圧 $p_{sat}(T)$

- **液 (Jacobsen, 式22)** [atm], $T_c=126.192$:
$$
\ln p_{sat}^l = \frac{n_1}{T}+n_2+n_3T+n_4(T_c-T)^{1.95}+n_5T^3+n_6T^4+n_7T^5+n_8T^6+n_9\ln T
$$
$n_1=8394.409444,\ n_2=-1890.045259,\ n_3=-7.282229165,\ n_4=0.01022850966,$
$n_5=5.556063825\times10^{-4},\ n_6=-5.944544662\times10^{-6},\ n_7=2.715433932\times10^{-8},$
$n_8=-4.879535904\times10^{-11},\ n_9=509.5360824$ (→ ×101325 で Pa)。
- **固 (Frels, 式23)** [mmHg]: $\log_{10} p_{sat}^s = 7.614676 - 356.281/T$ (→ ×133.322 で Pa)。

### 凝縮相密度 $\rho_l/\rho_s(T)$

- **液 (Nowak, 式24)**, $\tau=1-T/T_c$, $\rho_c=313.3$:
$$
\ln(\rho_l/\rho_c)=n_1\tau^{0.3294}+n_2\tau^{4/6}+n_3\tau^{16/6}+n_4\tau^{35/6}
$$
$n_1=1.48654237,\ n_2=-0.280476066,\ n_3=0.0894143085,\ n_4=-0.119879866$。
- **固 (Scott, 式25)** [kg/m³]: $\rho_s = 1068.49 - 1.97830\,T$。

### 潜熱 $L(T)$

- **液 (式26)** [MJ/kg]: $L_l = p_1T^4+p_2T^3+p_3T^2+p_4T+p_5$,
$p_1=-2.137\times10^{-8},\ p_2=7.18\times10^{-6},\ p_3=-9.142\times10^{-4},\ p_4=0.05069,\ p_5=-0.809$
(→ ×10⁶ で J/kg)。
- **固 (昇華, 式27)** [J/kg, 定数]: $L_s = 2.43\times10^{5}$。

### 表面張力 $\sigma(T)$

- **液 (Stansfield, 式28)** [dyn/cm], $T_c=126$, $\sigma_0=29.06$: $\sigma_l=\sigma_0(1-T/T_c)^{1.247}$ (→ ×$10^{-3}$ で N/m)。
- **固 (Dotson, 式29)** [N/m]: $\sigma_s = 2.0\times10^{-4}(\rho_s)^{2/3}$ ($\rho_s$ は式25)。

### Kelvin 効果と $T_d$ (Appendix 2)

液滴表面蒸気圧 (式32): $p_d=p_\infty(T_d)\exp\{2\sigma/(\rho_l RT_d r)\}$。$T_d$ の Hill 陰的式は
[4 節](#4-液滴温度-t_d--hill-のエネルギーバランス純蒸気簡約) 参照。$R_{N_2}=296.8$ J/(kg·K)。

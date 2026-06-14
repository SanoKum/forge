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
(N2 70 K 概算で $\xi_g\approx+5\times10^4>0$、凝縮で昇圧)。単相帰着検算: $g=0$ で $\kappa=\gamma-1$,
$\chi=0$ となり既存形に一致。保存量での $\partial p/\partial U$ と flux Jacobian への入り方、Roe/SLAU
への影響は [implementation.md](implementation.md) の対流フラックス節を参照。

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

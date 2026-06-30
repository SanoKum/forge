# On nitrogen condensation in hypersonic nozzle flows: numerical method and parametric study

**著者**: L. Lin, W. Cheng, X. Luo, F. Qin
**所属**: University of Science and Technology of China (USTC), Hefei, China (W. Cheng は現 KAUST)
**掲載**: Shock Waves (2014) 24:179–189, DOI 10.1007/s00193-013-0490-3
**受理**: 2013年12月

---

## 1. 背景 (Introduction)

- 極超音速風洞では、試験ガスとして広く使われる窒素 (N₂) がノズル内の強い膨張で**凝縮**しうる。凝縮が起きると**潜熱の放出**が流れ場に強く影響し、試験データの精度を損なう。よって信頼できる空力データを得るには相変化を**避ける**必要がある。
- 1951–1952 年頃から空気・窒素の凝縮について準1次元計算と実験が行われてきた (R.A.E. Farnborough, R.A.R.D.E. 等の極超音速風洞での Ma 6.8〜12.9 の研究、Dolgushev らによる Ma≥20 の数値・実験など)。
- 既往の数値研究の多くは数十年前の**簡単な1次元法**にとどまる。一方、Crane & Marshall は「空気を Ma≈5 程度まで加速する場合、凝縮回避のため予熱が必要になる」と指摘していた。
- **本研究の動機**: 凝縮の数値解法と物理モデルの進展を取り込み、**2次元平面・軸対称**の極超音速ノズル流れ + 窒素凝縮を計算できる近代的な手法を構築すること。これを設計マッハ数 6.0 のノズルに適用し、貯気室の全圧・全温を変えて凝縮の影響を評価する。

---

## 2. 方法 (Physical model & Numerical method)

### 2.1 支配方程式 — Euler + モーメント法
- 気相 (窒素ガス) と液滴の2相系を、**Euler 方程式 (非粘性)** に相変化の記述を加えて表す。
- 液相 (液滴雲) は、サイズ分布関数 $f(r;x,t)$ の**有限個のモーメント** $Q_0, Q_1, Q_2, Q_3$ で記述する (**method of moments**, Hill 1966 に倣う)。
  - $Q_0$: 単位質量あたりの液滴数密度、$(4\pi/3)Q_3$: 単位質量あたりの液滴体積。液相質量分率 $g = 4\pi\rho_l Q_3/3$。
- 保存変数ベクトル: $U = (\rho, \rho u, \rho v, \rho E, \rho g, \rho Q_2, \rho Q_1, \rho Q_0)^T$。
- 軸対称ソース項 $\delta H/y$ ($\delta=0$ 平面 / $\delta=1$ 軸対称)、および相変化ソース項 $S$ を含む。質量・運動量・エネルギーの保存式 (先頭4式) のソース項はゼロ (混合気体全体として保存)、液相モーメント式に核生成・成長のソースが入る。

### 2.2 **核生成モデル (nucleation)** ← 質問への回答
- **古典核生成理論 (Classical Nucleation Theory, CNT)** を採用 (Becker–Döring 1935、capillarity approximation で核生成障壁を計算)。
  $$ J_{CNT} = K_{CNT}\exp\left(-\frac{\Delta G^*}{k_B T}\right) $$
  - 速度係数 $K_{CNT} = \left(\dfrac{2\sigma}{\pi M \upsilon_l}\right)\left(\dfrac{p_v}{k_B T}\right)^2$、障壁 $\Delta G^* = \frac{4}{3}\pi r^{*2}\sigma$、臨界半径 $r^* = \dfrac{2\sigma}{\rho_l RT\ln(p_v/p_{sat})}$。
- **CNT への経験的補正関数**を掛ける。これは **Iland et al. (2009, J. Chem. Phys. 130:114508「Homogeneous nucleation of nitrogen」)** が膨張チェンバー実験 (1.3–14.2 kPa) から決めた窒素用の補正:
  $$ J = J_{CNT}\,\exp\!\left(A + \frac{B}{T}\right),\quad A = -55,\ B = 4{,}270\ \mathrm{K} $$
  - 本論文の圧力域 (2.8–6.0 kPa) はこの補正の検証域に含まれる (ただし高圧域には必ずしも適用できない、と注記)。

### 2.3 **液滴成長モデル (droplet growth)** ← 質問への回答
- **修正 Gyarmathy モデル**を採用。具体的には **Goodheart (2004, TU München 博士論文)** が遷音速の窒素凝縮向けに修正したもの。純窒素蒸気で広い Knudsen 数域に適用可能:
  $$ \frac{dr}{dt} = \frac{1}{\rho_l}\,\frac{f_{F\text{-}S}(Kn)\,k\,R T^2}{L^2}\,\ln\!\left(\frac{p}{p_{sat}(T_d)}\right) $$
  - $k$: 熱伝導率、$L$: 潜熱、補正係数 $f_{F\text{-}S}(Kn)=\dfrac{1+2Kn}{r(1+3.42Kn+5.32Kn^2)}\left(1-\dfrac{r^*}{r}\right)$。
- Young のより精緻なモデルと比較したうえで、全 Knudsen 域で妥当かつ式が単純な **Gyarmathy/Goodheart** を選択。Knudsen 数 $Kn=\lambda/2r$ で連続領域〜自由分子領域をカバー (核生成時の nm スケールから成長後の µm スケールまで)。
- **液滴温度 $T_d$**: Hill のモデル (Appendix 2) を採用し、液滴表面のエネルギーバランスから求める。凝縮係数 $\zeta=1.0$。液滴表面の蒸気圧 $p_d = p_\infty(T_d)\exp\{2\sigma/(\rho_l RT_d r)\}$ (Kelvin 効果)。

### 2.4 熱力学物性・実在気体
- 窒素を**熱量的完全気体**として扱う。van der Waals 状態方程式で圧縮係数を見積もると、本研究の条件 (11 MPa/500 K〜2.7 kPa/30 K) で実在気体効果は **<4%** と小さいため近似は妥当。
- 混合気体の内部エネルギー $e = (1-g)e_v + g e_l$ から温度 $T$・圧力 $p$ を導出 (潜熱 $L$、各比熱を含む式)。
- 窒素の物性 (飽和蒸気圧・液体/固体密度・潜熱・表面張力) は文献データに基づき**新たに精度の高いフィット式**を Appendix 1 に整備 (Jacobsen, Nowak, Frels, Scott, Stansfield 等のデータ・式)。

### 2.5 数値解法
- **fractional-step method** で支配方程式を「ソース項なしの均質部」と「凝縮ソース項を持つ非均質部」に分離。
- 均質部は Sun & Takayama の **VAS2D** (2次元・軸対称 適応ソルバ) と同じ手法で解く。非均質部 (ソース項) は Mundinger / Prast の扱いに従う。
- **有限体積法 (FVM)**、フラックスは **2次精度の MUSCL-Hancock スキーム**。
- **非構造四辺形メッシュ + 時間依存の適応細分化** (密度の誤差センサ = テイラー展開の2階/1階微分比を判定に使用)。

---

## 3. 検証 (Validation)

- **Arthur (1952, Caltech 博士論文)** の極超音速ノズル流れ + 窒素凝縮の実験と比較。
- ノズル形状 (検証用): **2次元 source-flow ノズル**、スロート高さ 0.254 mm → 出口高さ約 25.4 mm、幅一定 25.4 mm、全開き角 11°。壁面静圧をスロートから 6.35 mm 間隔のタップで計測。
- 条件は Arthur の Run 34-1 (全圧 $p_0$=844,037 Pa, 全温 $T_0$=290 K)。
- **結果**: 本手法の壁面静圧分布が実験とよく一致。計算・実験ともに、**凝縮により静圧が等エントロピー (凝縮なし) 計算より上昇**することを再現。

#### 検証ノズル形状の一次出典 (Arthur 1952 本文を確認)
- 出典: **Arthur, P.D., *On Condensation in Hypersonic Wind Tunnels*, PhD thesis, California Institute of Technology (GALCIT), 1952** (`papers/Arthur_pd_1952.pdf`)。**印刷 p.4「II. TEST FACILITY / A. Apparatus」**(PDF 15ページ目) に記述あり。スキャン PDF (テキスト層なし) のため画像描画で確認した。
- 原文 (要約): *"The hypersonic tunnel … was a **two-dimensional source-flow tunnel** expanding from a **throat height of 0.010 inch** to a **final height of approximately one inch** with a **constant width of one inch**. The **total expansion angle was 11°**."* / *"Nozzle wall static pressures were measured by 0.040-inch diameter taps drilled perpendicular to the diverging walls **at 1/4-inch intervals**."*
- 換算: スロート 0.010 in = 0.254 mm、出口 ~1 in = 25.4 mm、幅 1 in = 25.4 mm、全開き角 11° (半角 5.5°)、タップ間隔 1/4 in = 6.35 mm。→ 凝縮論文 (2014) の検証ケース寸法と完全一致。
- **形状の性質**: 「**source-flow ノズル**」= 仮想ソース点から放射状に広がる流れを前提とするため、**側壁は半角 5.5° の直線 (直線ウェッジ)**。曲線輪郭ではないので、**座標テーブル形式の輪郭は論文に無い**(スロート位置 + 半角 5.5° で形状が一意に決まる)。製作の詳細は Arthur 1952 の Ref.5 に委ねられている。
- **図**: Fig.1 = プラント概略 (ガス供給系)、Fig.2 = 装置写真、Fig.3/4 = スロートからの距離 [inch] や面積比 $A/A^*$ (〜100) に対する圧力・マッハ数分布。**ノズル輪郭そのものの寸法図は無く**、本文の言葉による記述が一次情報。
- **施設の補足**: Caltech / GALCIT の小型2次元凝縮試験用トンネル。本文 p.2 に出てくる「GALCIT Hypersonic Wind Tunnel」は別件 (飽和流研究) の装置で、この検証ノズルとは別物。

#### forge での dry 再現 (case/34)
- この Arthur ノズル形状を forge で再現し、**乾いた (dry) 非粘性の超音速膨張**を計算した: [`case/34.arthur_n2_nozzle/`](../case/34.arthur_n2_nozzle/README.md)。
- forge には**凝縮 (相変化) モデルが無い**ため、再現できるのは「凝縮なし」の等エントロピー的膨張 (Arthur Fig.4/11 の dry 1D 理論線) まで。凝縮による静圧上昇そのものは対象外。
- 結果 (`run_0002_slau_dry_cfl1`): 出口 **M=6.75** (等エントロピー 6.93, A/A*≈100)、中心線・壁の静圧が等エントロピー線に一致し中心線と壁も互いに一致 (Arthur Fig.11 の知見を再現)。
- 形状再現メモ: source-flow=直線壁なので座標表は不要。発散部のみ・上半分をモデル化し、入口に超音速インレット (M=1.05) を置いて choking を回避。スロートは鋭頂点だと凸コーナーの膨張特異点で発散するため、双曲線 `y=√(yt²+(x tanα)²)` で滑らかにし下流は 5.5° 直線に漸近させた (スロート曲率は原典に記載が無く、ここでは R_t≈13.7mm の代表値)。

---

## 4. 結果と考察 (Results)

### 4.1 凝縮の影響 (設計 Ma=6.0 軸対称ノズル)
- 基準条件: $T_0$=400 K, $p_0$=9 MPa。凝縮あり/なしを比較。
- 核生成ゾーン下流で液体窒素は **<1.5%**、ノズル出口で全窒素の **7.5%** まで増加。
- 凝縮ありでは潜熱放出により静温・静圧がともに上昇:
  - 出口平均温度: 40 K → **56 K**
  - 出口平均圧力: 2,800 Pa → **3,400 Pa**

### 4.2 パラメトリックスタディ (全温 $T_0$・全圧 $p_0$ を変化)
- **凝縮なし**: 出口静圧はほぼ一定、出口静温は全温に比例 (33–50 K)、Ma はほぼ一定 (≈6.67)。
  - (非粘性 + 純窒素のため出口 Ma=6.67 が設計値 6 より大きく出る、と注記。)
- **凝縮あり**:
  - 全温↑ → 出口の静圧・静温・凝縮量はいずれも減少、Ma は増加。
  - 出口静圧・静温は常に凝縮なしより**高い** (潜熱放出)。Ma は小さくなる (5.40–6.67)。
  - 全圧 9→11 MPa で出口静圧は 600–1,000 Pa 増、凝縮量は平均約 5 g/kg 増。全圧は凝縮量にほとんど影響せず、**全温の影響が支配的**。
  - 最低全温 330 K では、凝縮により出口の静圧・静温・Ma が設計値の最大 **71%・85%・19%** も変化。

### 4.3 Ma=5 ノズル (室温運転の検討)
- $p_0$=9 MPa、$T_0$=270–310 K。凝縮量は $T_0$=285/290/295/300/305/310 K で **23.44 / 15.7 / 8.0 / 2.7 / 0.27 / 0 g/kg**。
- **室温 (~305 K) の窒素を Ma=5 まで膨張させると弱い凝縮が発生**する → 予熱が必要。Crane & Marshall の判断を裏付ける。

---

## 5. 結論 (Conclusion)
- 凝縮性ノズル流れの数値手法 (CNT+経験補正 / 修正 Gyarmathy / モーメント法 + Euler の FVM) を構築・検証した。
- 凝縮量は全温の低下とともに増加し、出口の静圧・静温上昇と Ma 低下を招く。最低全温 330 K で設計値の最大 71/85/19% 変化。
- 室温で Ma=5 まで膨張させると弱い凝縮が起きるため、**極超音速ノズルを室温運転する場合は窒素の予熱が必要**。

---

## 質問への回答

### Q1. ノズル形状はどこの風洞か？
- **パラメトリックスタディ用 (設計 Ma=6.0 / Ma=5.0 の軸対称ノズル) の具体的な形状 (輪郭・寸法) は論文中に記載がありません。** 設計マッハ数と軸対称であること、出口での平均値しか与えられていません。USTC グループの想定する風洞ノズルと思われますが、施設名・輪郭は明示されていません。
- **唯一、具体的な形状・施設が分かるのは「検証ケース」**で、これは **Arthur (1952) が California Institute of Technology (Caltech) の極超音速風洞**で実験した **2次元 source-flow ノズル** (スロート 0.254 mm → 出口 25.4 mm、幅 25.4 mm、全開き角 11°) です。
- 背景で言及される他の風洞 (R.A.E. Farnborough 7"×7"、R.A.R.D.E. 極超音速ガス風洞など) は**既往研究の引用**であって、本論文の計算対象ではありません。

### Q2. 核生成・核成長モデルは何を使っているか？
- **核生成 (nucleation)**: **古典核生成理論 CNT** (Becker–Döring、capillarity 近似) に、**Iland et al. (2009) の窒素用経験補正** $J = J_{CNT}\exp(A+B/T)$ ($A=-55$, $B=4270$ K) を掛けたもの。
- **核成長 (droplet growth)**: **修正 Gyarmathy モデル** (Goodheart 2004 が窒素向けに修正、Knudsen 補正係数 $f_{F\text{-}S}(Kn)$ 付き)。全 Knudsen 域に適用可能で式が単純なため、Young の精緻モデルより採用。
- **液相の輸送**: モーメント法 ($Q_0$–$Q_3$、Hill 1966) を Euler 方程式に結合し有限体積法で求解。
- **液滴温度**: Hill のモデル (Appendix 2)。

---

## 補足: 液滴温度 $T_d$ の計算方法 (Appendix 2 の詳細)

### 位置づけ
$T_d$ は液滴表面の**エネルギーバランス**から求める。凝縮性流れ向けの液滴温度モデルは Johnson (2001) が10種を比較しており、本論文はそのうち **Hill のモデル** (refs 11, 28) を採用する。$T_d$ は成長則(式9)の飽和圧 $p_{sat}(T_d)$ と Kelvin 式(32)に入り、凝縮の自己抑制(潜熱で液滴が温まり成長駆動力が落ちる)を表す。

### 支配式
一般形 (式30):

$$
K_v\zeta\!\left(1-\frac{\beta_d T_d}{\beta T_v}\right) + \zeta\!\left(1-\frac{\beta_d}{\beta}\right)(L_v-1)\frac{\gamma}{\gamma-1}
= K_v(1-\zeta)\alpha\!\left[\frac{T_d}{T_v}-1\right] + \frac{\beta_c}{\beta}K_v\alpha_c\!\left[\frac{T_d}{T_v}-1\right]
$$

- 右辺第2項 ($\beta_c$ 項) は**不活性キャリアガス**の寄与で、**純蒸気では無視できる**。本論文は純窒素なのでこれを落とし、式(31) に簡約:

$$
K_v\zeta\!\left(1-\frac{\beta_d T_d}{\beta T_v}\right) + \zeta\!\left(1-\frac{\beta_d}{\beta}\right)(L_v-1)\frac{\gamma}{\gamma-1}
= K_v(1-\zeta)\alpha\!\left[\frac{T_d}{T_v}-1\right]
$$

### 各係数の定義
- $\zeta$: 凝縮係数。本研究では $\zeta = 1.0$。
- $\alpha$: 熱適応係数 (thermal accommodation coefficient)。
- $K_v = \dfrac{1}{2}\dfrac{\gamma+1}{\gamma-1}$
- $L_v = \dfrac{(\gamma-1)L}{\gamma R T_v}$ (無次元化した潜熱)
- $\beta = \dfrac{p_v}{\sqrt{2\pi R T_v}}$ (蒸気側の分子フラックス因子)
- $\beta_d = \dfrac{p_d}{\sqrt{2\pi R T_d}}$ (液滴表面側の分子フラックス因子)
- 液滴表面の蒸気圧 $p_d$ は **Kelvin 効果込み** (式32):
  $$ p_d = p_\infty(T_d)\,\exp\!\left\{\frac{2\sigma}{\rho_l R T_d\, r}\right\} $$
  ここで $p_\infty(T_d)$ は温度 $T_d$ に対応する**平膜の飽和圧**。

### 解き方(実装上の要点)
- 式(31) は $T_d$ について**陰的** (非線形)。$\beta_d$ が $p_d/\sqrt{T_d}$ を通じて $T_d$ に依存し、さらに $p_d$ も $\exp\{2\sigma/(\rho_l R T_d r)\}$ と $p_\infty(T_d)$ を通じて $T_d$・液滴半径 $r$ に依存する。したがって $T_v, p_v, r$ を所与として **$T_d$ を反復的に求める**。
- $\zeta = 1.0$ とすると右辺 $(1-\zeta)$ 項が消え、式(31) は
  $$ K_v\!\left(1-\frac{\beta_d T_d}{\beta T_v}\right) + \left(1-\frac{\beta_d}{\beta}\right)(L_v-1)\frac{\gamma}{\gamma-1} = 0 $$
  という $T_d$ の陰的方程式になる。第1項は液滴⇄蒸気の正味の質量(分子)フラックスに伴う温度差、第2項は潜熱放出の寄与を表し、両者の釣り合いで $T_d$ が決まる(潜熱で $T_d > T_v$ となり、過飽和を下げて成長を自己抑制する)。

---

## 補足: 液滴温度 $T_d$ は「核生成」と「核成長」のどちらに効くか

### 区分けの原則 — 使う温度そのものを分けている
本論文は、$T_d$ を生成・成長に按分しているのではなく、**駆動力 $\ln(p/p_{sat})$ の飽和圧をどの温度で評価するかを過程ごとに変える**ことで区別している。

| 過程 | 駆動力 | 飽和圧の評価温度 | 理由 |
| --- | --- | --- | --- |
| 核生成 $J$ | $r^*=\dfrac{2\sigma}{\rho_l RT\ln(p_v/p_{sat})}$, $\Delta G^*\propto r^{*2}$ (式2–5) | **気相温度 $T$** | 臨界核は sub-nm。生成時の潜熱は即散逸し自己加熱が無いとみなすため、CNT は気相(浴)温度で定義 |
| 核成長 $dr/dt$ | $\dfrac{dr}{dt}\propto\ln\dfrac{p}{p_{sat}(T_d)}$ (式9) | **液滴温度 $T_d$** | 有限サイズの液滴は蓄積潜熱で $T_d>T_v$ になる |

- すなわち **$T_d$ は成長則(式9)・Kelvin 式(32)にしか直接は入らない**。核生成側(式5)は気相温度 $T$ のまま。
- $T_d$ を決めるエネルギーバランス式(31)・Kelvin 式(32)は有限半径 $r$ を前提($2\sigma/(\rho_l RT_d r)$、$\beta_d$ が $r$ 依存)とするため、**まだ存在しない臨界核には構造的に適用できない** → $T_d$ は成長専用の量になる。
- モーメント方程式のソース項 $S$ でも、**生成項 $\propto J\,r^{*n}$($T$ 依存)** と **成長項 $\propto n\,\rho Q_{n-1}\,dr/dt$($T_d$ 依存)** が加算的に分離している。

### $T_d$ が生成に効くのは間接的(バルク $T$ 経由)だけ
液滴温度が生成に及ぼす影響は、$J$ に $T_d$ を直接入れる形ではなく、

$$ \text{成長で潜熱放出}\;\to\;\text{バルク気相温度 }T\text{ 上昇(エネルギー式18)}\;\to\;p_{sat}(T)\uparrow\;\to\;\text{過飽和}\downarrow\;\to\;J\downarrow $$

という **バルク $T$ を介した間接ループ**で入る。局所の液滴自己加熱(①)は成長抑制、バルク昇温(③)は生成抑制、と棲み分けられている。

### Kantrowitz 補正との関係(重要)
- **Kantrowitz の非等温補正 (non-isothermal correction)** は、まさに「**生成中の臨界核が放出潜熱を捨てきれず自己加熱し、核生成率 $J$ が下がる**」効果を**核生成項そのものに**織り込む補正である。形式的には $J_{noniso}=J_{iso}/(1+b\,q^2)$ のような係数で $J$ を割り引く。ご認識のとおり、これは**核生成側に入れる**考え方。
- **本論文はこの Kantrowitz 補正を採用していない**(本文・付録に "Kantrowitz" / "non-isothermal" / "isothermal" の語は一度も登場しない)。代わりに核生成は **等温 CNT × Iland (2009) の経験補正 $\exp(A+B/T)$** で扱い、自己加熱の物理は**成長側の $T_d$ にのみ**持たせている。
- したがって整理すると:
  - **Kantrowitz** = クラスタ自己加熱を**生成率 $J$ に直接**織り込む(生成側の補正)。← 本論文は不採用。
  - **本論文の $T_d$** = 既存液滴の自己加熱を**成長則 $dr/dt$ に**織り込む(成長側のモデル)。生成への影響はバルク $T$ 経由の間接のみ。
- つまり「自己加熱を核生成へ入れる」発想(Kantrowitz)と、本論文が実装している「自己加熱は成長側 $T_d$ で扱う」発想は**別物**であり、本論文では前者を採らず後者だけを採っている、というのが答え。

---

*まとめ作成元 PDF: `papers/on nitrogen condensation in hypersonic nozzle flows.pdf`*

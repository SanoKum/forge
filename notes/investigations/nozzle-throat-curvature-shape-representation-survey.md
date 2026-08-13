# スロート近傍の二階微分 (曲率) を制御できる形状表現 — 文献調査

- 依頼元: [notes/sessions/2026-08-17-throat-curvature-shape-representation-survey.md](../sessions/2026-08-17-throat-curvature-shape-representation-survey.md)
- 対象: `case/41.wind_tunnel_design` ①軸対称風洞 (モード F)、$M_d=4$ / $r_t=10$ mm / 円弧 $R_u=R_d=2r_t$
- 関連計画: [plans/active/tooling-nozzle-phase3-windtunnel.md](../../plans/active/tooling-nozzle-phase3-windtunnel.md) §9.2 (B5–B10)
- 調査日: 2026-08-17 (実施 2026-08-13 セッション)。一次資料 PDF 7+3 本の精読 + Web 書誌調査
- 凡例: **[事実]** = 一次資料の記載 (ページ付き)、**[導出]** = 文献式からの本調査の計算、**[推測]** = 解釈・適用判断

---

## 0. 要旨

1. **古典 (Sivells/CONTUR) は「壁側で曲率を接合する」ことをそもそもしない**。曲率連続の壁は、軸中心線の速度/マッハ数分布を「遷音速解・radial flow・出口一様条件に **1–2 階微分まで整合させた区分多項式**」で規定し、壁を MOC の従属出力として一体生成することで**帰結として**得る。壁幾何側の $r,r',r''$ 突き合わせという操作自体が存在しない (§2.1)。本ケースの「円弧 + J で C2 接合した設計壁」は、CONTUR が序論で批判した Foelsch 型の系譜に構造が近い。
2. **円弧–下流壁の接合部が軸近傍に波を作る現象は 1963 年に文献化されている** (Darwell & Badham; Back & Cuffel 1966 が実験検出; Cuffel 1969 が接線接合直下流の壁圧上昇→Mach 線収束→弱斜め衝撃波を実測)。軸対称では壁擾乱が軸に集束する (Sivells) (§2.4)。
3. **$R=2r_t$ は遷音速理論の適用下限ちょうど** (Cuffel: "not less than about 2"、判定基準は壁静圧のみ)。Sauer 一次解の軸上速度は**厳密に $x$ の 1 次**で、軸 $M''$ に幾何情報が構造的に存在しない — 本ケースの「$M''$ 2.3 倍ずれ」は誤差ではなく**次数の欠落** (§2.3, §9)。CONTUR の実設計例は $R=5.5$–$6$、$R<12$ で 3 階微分補正が必要 ($R=1$ で 16%)。
4. **こぶ (x/rt≈1.43, 到達不能域内) は J 下流の壁表現をどう差し替えても消えない** — 波源の C⁻ 壁足は凍結円弧上 ($x_w=0.154 < x_j$) にあり、設計壁のどの自由度も原理的に届かない (B4 実測)。G3 化・CST 化・κ(s) 化を **J 下流だけに**適用するのは的外れ (§6)。効きうるのは (a) $R$ を上げる、(b) スロート下流の円弧自体を設計出力に置き換える (CONTUR 型接合廃止)、(c) 遷音速アンカーの高次化、の 3 経路のみ。
5. **推奨** (§8): 最小判別実験として **R=5 の単発 A/B** (表現・チェーン不変、既存 dv の変更のみ) を先行。文献 (Cuffel 下限・CONTUR 実例・Sauer 摂動 $\propto 1/\sqrt{R}$・自前 R スイープ) がすべて $R=2$ を「低 R 帯」と判定しており、こぶ振幅の $R$ 依存を測るのが最安の分岐点。表現差し替えが必要になった場合の本命は **κ(s) 直接パラメトライズ (Korakianitis CIRCLE 型)** で、円弧は「κ=const 区間」として自然に埋め込める。
6. **参考: 到達水準の相場** — Korte の実在気体 Sivells 設計でも CFD 検証で **~0.02 Mach の残留変動**があり、原因は「設計 MOC の遷音速モデルと実流れの差」に帰属されている **[事実]** (§2.5)。本ケースのこぶ振幅 0.019–0.023 は古典設計チェーンの実測残差と同水準。ただしこれは「不可避」を意味しない (Korte の帰属はまさに本ケースが CFD アンカーで潰しにいっている誤差源であり、B10 の撤回と整合)。

---

## 1. 前提の訂正 — ローカル `papers/` の棚卸し

依頼プロンプトの「Sivells / CONTUR の一次資料はローカルに無い」は**誤り**。以下が既にある:

| ファイル | 実体 |
| --- | --- |
| `papers/aerodynamic design program for axisymmetric and planar nozzles (CONTUR, Sivells AEDC-TR-78-63).pdf` | **CONTUR 本体レポート** (151p、精読済) |
| `papers/general characteristics of flow through nozzles at near critical speeds (Sauer NACA TM-1147).pdf` | Sauer 一次解の原典 |
| `papers/inviscid design of hypersonic wind tunnel nozzles for a real gas (AIAA 2000-0677).pdf` | Korte による Sivells 法の現代的解説 + 実在気体拡張 |
| `papers/preliminary nozzle design for a small-scale high-enthalpy facility (OSTI 1568032).pdf` | 卓上 HWT 予備設計 (Shope 文献へのポインタ源) |

また `papers/Korte explicit upwind algorithm ... (NASA TP-3050).pdf` は**ノズル設計論文ではない** (PNS ソルバのアルゴリズム開発のみ。全文走査で design/contour/throat の記載なし)。Korte の CFD-in-the-loop 設計の正本は AIAA 92-4009 (CAN-DO) / Korte & Hodge JSR 1995 で、ローカルに無い (§10)。`papers/Arthur_pd_1952.pdf` は N₂ 過飽和の凝縮実験 (Caltech PhD) で本件と無関係 (`case/16.nozzle_wys` の文脈では有用)。

---

## 2. Q1: 古典的ノズル設計法は throat 近傍の曲率をどう扱っているか

### 2.1 Sivells / CONTUR (AEDC-TR-78-63, 1978) — 曲率連続の実装法

ページは刷りページ (PDF ページ = +3)。

**[事実] 連続化の作用点は壁ではなく軸分布**。Abstract の宣言:

> "The continuous curvature is achieved through specification of a centerline distribution of velocity (or Mach number) which has first and second derivatives that 1) are compatible with a transonic solution near the throat and with radial flow near the inflection point and 2) approach zero at the design Mach number."

- 一致箇所は 3 つ (Summary p.41): **throat characteristic の軸交点 I** / radial flow 域の両端 E, B / test section 始端 C。
- 上流遷移 I–E は速度比 $W$ の 5 次多項式 (Eq. 35) で **両端 $W, W', W''$ の 6 条件** (p.16–17)。推奨は区間長を反復で決めた 3 次。**3 階微分の一致はオプション** (入力 `IX`: +1=遷音速解 / −1=source flow / 0=不一致)。下流遷移 B–C は $M$ の多項式で B 側 radial flow の $M,M',M''$、C 側 $M'=M''=0$ (p.17)。
- **曲率連続を要求する理由** (p.6): ゼロ渦度条件 $\partial q/\partial n=\kappa q$ から "the flow will be disturbed where a contour has a discontinuity in curvature." — plan §9.2 の訂正済み理解 ($\partial p/\partial n=\rho V^2\kappa$) と同型。
- **軸対称の擾乱集束** (p.7): "disturbances created by imperfections in the contour tend to be focused on the centerline." — 「壁の欠陥→軸のこぶ」という本ケースの観測経路そのもの。
- **軸勾配不連続 ⇔ 壁曲率不連続の対応** (p.5, 7): Foelsch 法は radial flow と一様流の接合で軸 M 勾配が不連続になり "This discontinuity produces a discontinuity in curvature of the contour at the inflection point."

**[事実] throat 壁に円弧などの解析形状を使わない**。遷音速域は Hall (1962) の $1/R$ 逆冪級数を Kliegel–Levine (1969) の $S=R+1$ 置換で低 $R$ 拡張した解 ($u$ は $x^3$ 打切り、p.10–14)。そこから**スロート壁点発の右行特性線 (throat characteristic)** を積分し、軸交点 I 以降を MOC で解く。壁は**各特性線を横切る質量流束の積分**で決まる点列 (p.23, 36) + 製作出力用 cubic spline (p.34)。1978 年版の主改良は「sonic line 接続で生じていた壁の未計算 gap ($R$ 小ほど拡大) を throat characteristic 接続で解消」(p.8–9)。

**[事実] 低 R の限界の定量**:

- $u$ 級数の $x^3$ 打切りにより "the accuracy of the transonic solution is limited, particularly for low values of R"。**$R<12$ では $d^3u/dx^3$ を補正** ("The correction is about 16 percent for R = 1 and decreases rapidly as R increases", p.14)。
- 音速点の $d^2W/dx^2$ の符号が **軸対称 $R=10.525$** (planar 11.767) を境に変わる (p.16)。
- 実設計例: "A recent (1975) design of a Mach 6 nozzle ... throat radius of curvature of about 5.5 times the throat radius" (p.9)、サンプル Mach 4 設計は **RC=6** (p.39–40)。極端な小 R は "would not normally be recommended"、逆に $R=\infty$ (uniform throat) はノズルを不必要に長くする (加速度 $\propto 1/\sqrt{R}$、p.10)。

**[推測] 本ケースへの読み替え**: CONTUR の観点では本ケースの接合問題は「壁幾何の C2/C3 の問題」ではなく「**J を通る特性線上の流れデータ ($M,\theta$ とその微分) の不整合**」であり、それが軸集束でこぶになる。B5 で「CFD 実測アンカー化が最も効いた (58% 減)」という実測はこの読みと正確に整合する。

### 2.2 Foelsch / Rao / Beckwith系

- **Foelsch (1949)** J. Aero. Sci. 16(3) 161–166: 点源 radial flow 域 + 一様流への遷移解析式。遷音速解を持たず radial flow がスロート付近から始まると仮定 — Beckwith らが「かなり不正確」と実証 (NACA TN 2711, 1952) **[事実、二次資料]**。接合は C1 級で、Sivells が曲率不連続を名指しで批判 (§2.1)。なお forge の `moc_inverse._CPlusMarch` docstring の「Foelsch/CONTUR 型 C⁺ マーチ」という記述は、CONTUR の実装 (throat characteristic + 質量流束積分) とは別物なので用語に注意 **[事実、コード照合]**。
- **Rao TOC (1958) / TOP (1960)**: ロケットノズルの推力最適化。定番の $R_u=1.5r_t$, $R_d=0.382r_t$ は**原典の明記箇所を Web では確定できず** (Sutton 教科書経由で定着した可能性) **[推測]**。ロケットは出口一様性が要求されないため小 $R_d$ で短縮化する — 風洞と定石が逆になる理由。
- **Beckwith/Moore/Chen (NASA Langley 静粛ノズル)**: 主敵は壁凹面曲率の Görtler 渦と壁 BL の radiated noise。緩膨張 + radial flow 域で quiet core を延長 (Chen-Malik-Beckwith 1990)。曲率分布は波だけでなく**境界層遷移経路も**支配する、という別軸の曲率制約 **[事実、抄録ベース]**。

### 2.3 遷音速スロート解と $R_d$ の定量 (Cuffel–Back–Massier / Sauer / Szaniszlo)

**Cuffel, Back & Massier (AIAA J 7(7) 1969, pp.1364–1366)** — 全文精読:

- **[事実]** 支配パラメータは $r_c/r_{th}$。"**Provided that rc/rth is not less than about 2**, existing two-dimensional flow theories [Sauer; Oswatitsch–Rothstein; Hall] adequately predict the transonic flowfield **as indicated by the wall static pressure measurements**" (p.1364)。→ $R=2$ は妥当域の**下限ちょうど**であり、しかも判定基準は壁静圧。starting line の Cauchy データに要る $M', M''$ の精度は保証外 **[推測]**。
- **[事実]** $r_c/r_{th}=0.625$ ではスロート面で軸 $M=0.8$ / 壁 $M=1.4$、軸静圧は壁の最大 3 倍。
- **[事実]** 円弧→円錐の**接線接合直下流で壁静圧が上昇**し ("compressive turning")、Mach 線が収束して弱斜め衝撃波を形成、軸を $z=3.15$ in で横切る (pp.1364–65)。
- **[事実]** 古典実務の starting line: "Kliegel used a variable-Mach-number start line that varied from **M = 1.3 along the wall to M = 1.1 along the centerline**; Shelton used a start line of M = 1.3" (p.1366) — 本ケースの $M=1.05$ 縦線より下流寄り。

**Sauer (NACA TM-1147, 1947)**:

- **[事実]** 軸上速度は**仮定として** $u_0(x)=ax$ (式 11) — べき級数の第 1 項のみ。$a=\sqrt{2/((k{+}1)\rho_S y_S)}$ (式 20b)。
- **[事実 (式の構造)]** よって軸上 $M''$ にノズル幾何の $x^2$ 情報は**ゼロ** ($u\to M$ の非線形変換由来の分のみ)。本ケースの「Sauer $M''$ が CFD と 2.3 倍ずれる」は**構造的必然であり、係数調整では直らない**。
- **[導出]** 幾何スロート面の壁摂動 $u_{wall}=1/(4R)$。$R=2$ で **0.125** — 「微小量」の線形化前提が既に苦しい。$M'$ +25% は一次打切りの想定内 (1D 近似は $\sqrt{(k{+}1)/2}\approx$ +10% 過大、と本文にも序列あり)。

**Szaniszlo (NASA TN D-7848, 1975)**:

- **[事実]** 連続壁曲率ノズル ($R_N=2.29$) の実測音速 $C_d$ = 0.989–0.991 (BL 込み)。"the sonic CD for any nozzle with a continuous and finite radius of curvature (**R_N > 2**) ... can be expected to be represented by the previous analytical method" — ここでも $R>2$ が扱える標準クラスの境界。
- **[導出]** Sauer 一次解を $R=2$, $\gamma=1.4$ で積分すると非粘性 $C_d=0.9934$。**本ケース CFD の mdot 比 0.994 は 0.06% 差で整合 — 妥当と確認**。
- **[事実]** "Any instability is most likely promoted by the local adverse pressure gradient generated by the **discontinuous radius of curvature at the juncture**" (p.8) — 曲率不連続接合が局所擾乱源という傍証。

**高次遷音速解の書誌** (参照文献リストから抽出): Hall, QJMAM 15(4) 1962, 487–508 (摂動パラメータは $1/R$ そのもの)。Kliegel & Levine, AIAA J 7(7) 1969, 1375–1378 ($1/(R{+}1)$ 展開へ変換、$R=2$ で摂動 1/3)。Kliegel & Quan, AIAA J 6(9) 1968。Hopkins & Hill, AIAA J 4(8) 1966。Back–Massier–Gier, AIAA J 3(9) 1965 ($R\ge2$ の壁圧検証の根拠文献)。

### 2.4 接合部起因の波 — 直接の古典文献

| 文献 | 内容 **[事実]** |
| --- | --- |
| **Darwell & Badham**, AIAA J 1(8) 1963, 1932–1934 | MOC 計算で「負族 Mach 線が**スロート円弧と円錐の接合点またはその直下流**から発し、**軸近傍で交差して衝撃波化**」。接合部近傍の壁形状修正で除去可能と結論。本ケースの現象の 1963 年版 |
| **Back & Cuffel**, AIAA J 4(12) 1966, 2219–2221 | 「円弧スロート+円錐」接合が作る斜め衝撃波の実験検出 |
| Cuffel–Back–Massier 1969 | §2.3 の通り (接線接合直下流の壁圧上昇を実測) |
| **Vallabh & Skews**, R&D Journal (S. Afr.) 33, 2017 (open access) | CSIR 超音速風洞の既存ノズルが「試験部に弱い波」を生む問題を「MOC 設計壁の変曲点曲率不連続」と特定し、**Sivells 法採用で解決** — 「曲率不連続→波→流れ品質劣化」の現代の実務実証 |

これらは G1 (接線) 接合 = **曲率が跳ぶ**ケースであり、本ケース (G2 済、$\kappa'$ が跳ぶ) より 1 階強い特異性である点に注意 **[事実の整理]**。$\kappa'$ 不連続 (G3 破れ) を G2 破れと分離して定量した超音速ノズル研究は**見つからなかった** (翼列の CIRCLE 系のみ、§3.2)。

### 2.5 Korte の実在気体 Sivells 設計 (AIAA 2000-0677) — 到達水準の相場

- **[事実]** 4 領域分割 (subsonic-throat / throat→conical / source flow / turning) で、亜音速壁は「スロート近傍は円弧相当・上流へ**曲率が指数減衰**する単一解析式」(紙面 p.3) — 亜音速側も接合を作らない。
- **[事実]** 遷音速解析用のダミー壁でさえ "preserves **continuity at the throat wall in radius, slope, and curvature**" を要求 (紙面 p.4)。超音速側は音速線に影響しない (双曲性) ため任意。
- **[事実]** 設計例 (M14, $R=3.0$) の CFD 検証で "slight **Mach number variations on the order of 0.02 Mach**"、帰属は "most likely due to the differences in the **transonic expansion solution** for the nozzle throat used in the design ... **characteristic of Sivells' design process**" (紙面 p.5)。→ 曲率連続を完全保証した古典チェーンでも、遷音速モデル誤差起因で 0.02 級のうねりが残るのが実測相場。
- **[事実]** CFD 帰還で壁を直接更新する系譜の正本: CAN-DO (AIAA 92-4009)、Korte & Hodge JSR 32(4) 1995、**Shope AEDC-TR-04-2 (2005) "Design Optimization of Hypersonic Test Facility Nozzle Contours Using Splined Corrections"** — 最後は「既存コンタ+スプライン補正を CFD で最適化」で、本ケースの帰還ループの直系先行例 (未入手、§10)。

---

## 3. Q2: 曲率・曲率微分を直接制御する形状表現

### 3.1 G2/G3 の定式化と本件での具体形

**[事実]** Gn 連続 = 再パラメトライズを許した Cn (β-constraints; Barsky & DeRose 1986, Farin CAGD 5th ed.)。平面曲線で G2 ⇔ 位置・接線・$\kappa$ 連続、**G3 ⇔ さらに $d\kappa/ds$ 連続**。

**本件への具体形 [事実の適用]**: 円弧は $\kappa(s)=1/R$ 一定なので $d\kappa/ds\equiv 0$。したがって **J での G3 化 = 設計壁スプライン側に「J で $d\kappa/ds=0$」を課すこと** (現行 `ModeFWall` の 5 次 B-spline は左端 $(r',r'')$ クランプ = G2 まで。$r'''$ 相当が未拘束)。5 次 B-spline は端条件を 1 本追加する自由度があり (左端 3 条件 + 右端 2 条件、または LSQ 化)、実装は小さい。ただし G3 にしても G4 ($\kappa''$) の跳びは残る。

### 3.2 κ(s) を設計変数にする表現 — 本命候補

**Korakianitis CIRCLE 法** (prescribed surface curvature distribution):

- 書誌: ASME J. Turbomach. 115(2) 1993; Computer-Aided Design 25(5) 1993; Applied Energy 89(1) 2012; ASME J. Turbomach. 135(4) 2013。
- **[事実、抄録+open access 検証研究]** 翼面の $\kappa(s)$ と**その勾配**を LE から TE まで連続に保つよう構成的に保証し、座標は $\kappa$ から積分で生成。曲率不連続が作る表面 Mach 分布の異常 (スパイク・こぶ) を除去、剥離泡縮小・損失低減を CFD/実験で実証 (Shen et al., J. Algorithms & Comput. Tech. 11(1) 2017, DOI 10.1177/1748301816665527 — 本文確認済)。
- **本件への移植形 [推測]**: ノズル壁を $\kappa(s)$ の区分スプライン (C1) で与え、$\theta(s)=\int\kappa\,ds$、$(x,r)=\int(\cos\theta,\sin\theta)\,ds$ で生成。**円弧は「$\kappa=$const 区間」として自然に埋め込まれ、接合の滑らかさは $\kappa(s)$ の関数クラスの選択だけで決まる** ($\kappa(s)$ を C1 にすれば壁は自動的に G3)。逆 MOC 点列へは $\kappa(s)$ 係数の LSQ フィットで追従。
- 系譜の語彙: クロソイド ($\kappa=as$、G2 止まり)、log-aesthetic curves (Miura; 超音速内部流適用例なし)、Class A Bézier (Farin CAGD 23(7) 2006 — Cao & Wang 2008 の反例と修正条件に注意)、MVC (Moreton & Séquin 1992: $\int(d\kappa/ds)^2ds$ 最小化 — 5 次要素と相性が良い fairing 汎関数)。

**曲率不連続→表面圧スパイクの実証 (翼列)**: Goodhand & Miller, ASME J. Turbomach. 133(2) 2011 — LE blend 点の曲率不連続が圧力スパイクを作り、除去でプロファイル損失低減・作動範囲拡大 **[事実、二次引用群]**。壁圧分布のこぶという現象の機序が本件と同型 **[推測]**。

### 3.3 CST (Kulfan AIAA 2007-62) — 精読結果

- **[事実]** $\zeta(\psi)=C^{N1}_{N2}(\psi)\,S(\psi)+\psi\zeta_T$、$S$ は Bernstein 展開。端点曲率 (LE 半径) は $S(0)=\sqrt{2R_{LE}/c}$ で**陽に制御** (p.5)、先頭係数のみが担う (p.10)。上下面 $Bu_0=Bl_0$ で前縁まわり曲率連続 (p.32)。
- **[事実 (不在の確認)]** **$\kappa'$ / G3 の明示制御は論文に一切登場しない**。"curvature-controlled CST" 系の後続専用論文も Web 調査で見つからず。
- **[事実]** 内部流路適用はダクト・ナセルインレット (スロートあり) の表現例まで (p.20–25)。次数は再現に BPO6–9、設計変数としては BPO4–6 で十分 (p.15)。
- **[推測]** CST が本件に効くとすれば理由は「制御できるから」ではなく「**単一 $C^\infty$ 関数なので接合が存在しないから**」。$r(\psi)$ は係数に線形なので、スロート内部点拘束 ($r=r_t$, $r'=0$, $r''=1/R$ 相当) は線形等式拘束で厳密に課せる (論文外の自然な拡張)。ただし任意点の $\kappa,\kappa'$ の狙い撃ちは非線形制約になり、κ(s) 表現より間接的。

### 3.4 補間 vs 近似 — B-spline の曲率振動

**[事実、CAGD 標準知見]** 点列を「通す」高次補間は曲率に波を作りやすい。定石は (a) LSQ 近似化 (自由度削減)、(b) fairing 汎関数 ($\int|r''|^2$ や $\int\kappa'^2$) の付加、(c) ノット配置の調整 (CAD 27(1) 1995 ほか)。**現行 `ModeFWall` は逆 MOC 点列の全点通過補間** (`make_interp_spline` k=5) なので、点列の離散雑音がそのまま $\kappa$ 波になる既知の構成。plan B2 の「最小二乗フィット」案は文献的に正当 **[事実の適用]**。ただし壁曲率の小波が本ケースのこぶの主因でないことは到達可能性の実測から明らか (§6)。

---

## 4. Q3: 「接合部を作らない」代替設計

- **CONTUR がまさにそれ** (§2.1): スロート壁点から下流の壁は全て MOC 出力で、固定円弧区間が超音速側に存在しない。「単一の滑らかな曲線で表す」のではなく「**壁を区分表現で持つこと自体をやめる**」のが古典の解。亜音速側も Korte の「曲率指数減衰の単一式」(§2.5) で接合レス。
- **starting line の実務**: 本ケースの「$M=1.05$ の縦線、壁足 = 円弧終端 = 可変壁始点」という構成は、調査した範囲の古典実務に**対応物がない** **[事実 (不在の確認)]**。古典は (i) CONTUR: throat characteristic (特性線) を初期線にし、旧 sonic line 方式の「壁の未計算 gap」を名指しで解消 (p.8–9)、(ii) Kliegel/Shelton: 壁 M=1.3/軸 M=1.1 の可変 M start line (§2.3)。縦線 starting line は特性線でないため厳密には Cauchy データの適合条件が異なり、壁足がそのまま可変壁始点になる本構成は「gap を可変壁で埋めている」構造に相当する **[推測]**。
- **教訓 (OSTI 1568032)**: minimum-length ノズル設計に収縮部を付けたら「設計法が湾曲音速線を考慮しないため衝撃波が立った」(p.20) — 遷音速整合を欠いた接合の失敗例。

---

## 5. Q4: 二階微分の「測定」側

- **[事実]** 古典設計法は $M''$ を「測る」のではなく「規定する」(多項式係数として。§2.1)。測定問題は CFD 帰還設計に固有。
- **[事実、標準手法]** 雑音下の 2 階微分推定の定石: **Savitzky–Golay** (局所多項式 LSQ の解析微分、`scipy.signal.savgol_filter(deriv=2)`; Anal. Chem. 1964)、**Tikhonov 正則化数値微分** (Knowles & Renka, EJDE Conf. Proc. 21; Chartrand, ISRN Appl. Math. 2011 は TV 系で不連続微分向き — 滑らかな $M''$ には Tikhonov/H² が適合)。窓幅/λ は L-curve で選ぶ。
- **本ケースの実践との照合 [事実]**: B10-a/B10-b の結論 (cell は atomicAdd 雑音床 ~6e-4 で $\epsilon/dx^2$ が $M''$ 信号を圧倒、node で 124× 改善し $M''=+0.19\pm0.03$) は、この標準論とちょうど整合する。「局所多項式フィット + 次数/窓幅感度の記録 + node データ」という現行手順 (B10 正本) は文献的にも適正で、追加すべきは高々 SG/Tikhonov による頑健化のみ。

---

## 6. 本ケースへの適用分析 — 表現差し替えの効く範囲

実測事実の再確認 (plan §9.2):

1. こぶ頂 ($x/r_t\approx1.43$) は**到達不能域内** ($x<x_{reach}=1.704$)。3 種の壁で軸 M が 4 桁一致 — 設計壁の自由度は届かない。
2. 波源ピークの C⁻ 壁足は $x_w=0.154$ = **凍結円弧上** ($x_j=0.2465$ の上流)。
3. J の $r,r',r''$ 連続化 (B2 相当) はこぶ不変 (0.1176→0.1152)。
4. 遷音速 Cauchy データの CFD 置換が最も効いた (58% 減)。

**帰結 [推測、ただし 1–4 は実測]**: こぶの支配源は「J の幾何連続性の階数」ではなく「**円弧区間そのものが作る遷音速〜低超音速の実流れ構造**と設計側モデルの残差」である。したがって:

- **J 下流の壁表現の差し替え (G3 クランプ・CST・κ(s)) を単体で行ってもこぶは消えない** — これは文献調査の前に判っていたべき帰結であり、本調査はそれを裏書きした (CONTUR の「接合は流れデータの問題」、Korte の 0.02 帰属)。
- 効きうる経路は 3 つ:
  - **(a) $R$ を上げる** — 円弧が作る 2 次元性 ($u_{wall}=1/(4R)$) と Sauer/CFD 残差 ($\propto\bar u\propto1/\sqrt R$、自前 R スイープで確認済) を源から減らす。文献の全会一致 (Cuffel 下限 2 / CONTUR 実例 5.5–6 / $R<12$ 補正 / 音速点 2 階微分符号 10.5)。副作用: ノズル長増 (加速度 $\propto1/\sqrt R$)、$x_{reach}$ 移動。既存 dv なので**表現・チェーン変更ゼロ**。
  - **(b) スロート下流円弧の廃止 (CONTUR 型)** — 可変壁 (または設計出力壁) の始点を J からスロート直後へ移し、遷音速アンカー (CFD) から throat characteristic を計算してそこから下流を一体生成。接合が消え、こぶ帯全域が設計到達域に入る。実装は `moc_inverse` の初期線の変更 + `ModeFWall` の区分再定義で、チェーン構造に踏み込む (実装コスト大)。既存 WIP `_CPlusMarch` が方向としては近い。
  - **(c) 遷音速アンカーの高次化 (Kliegel–Levine) + starting line の下流移動 (M=1.1–1.3)** — B6 残課題そのもの。古典実務 (Kliegel/Shelton) の実績構成。CFD アンカーが既にあるため増分は限定的だが、(i) ブートストラップ品質、(ii) CFD レス高速プレビュー、(iii) 低 R 対応で価値が残る。
- **§4 のジレンマ (始端 C2 整合 vs 出口一様性) は壁表現ではなく目標曲線の基底の問題**。自由 CP 3 個の Bernstein グローバル基底は台が広く ($x/r_t=1.38$ で遠方 CP2 が重み 11% — B7 実測)、始端拘束の変化を局所で吸収できない。対策は CP 増 (3→4–5、B10 の決定的射影で初期化) か、局所基底 (B-spline target) 化。EHVI への影響: 次元 3→5 は EHVI の実用域内だが評価回数は増える (次元あたり経験的に 1.5–2×) **[推測]**。

---

## 7. 採用候補の比較表

| 案 | 表現/変更 | 保証される連続性 | 追加 dv | 逆 MOC との相性 | 実装コスト | こぶへの効果見込み | リスク |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **A. R≥5 化** (表現不変) | 既存 dv `R` の値/下限変更 | 現行のまま (G2) | 0 | 変更なし | **極小** (config + 再設計 1 パス) | **高** (文献全会一致、源を減らす) | ノズル長増・MOO の探索域縮小; $x_{reach}$ 移動で target 再構築 |
| B. G3 クランプ | 5 次 B-spline 左端に $d\kappa/ds=0$ 追加 (または 7 次化/LSQ 化) | G3 | 0 | 変更なし | 小 | **低** (波源は円弧上・到達不能域) | 端条件過拘束で最初のノット区間が暴れる可能性 |
| C. κ(s) 直接表現 (CIRCLE 型) | 壁全体を $\kappa(s)$ C1 スプライン化、円弧= κ=const 区間、点列へ LSQ | 任意 (κ(s) のクラスで G3+) | 0 (係数は従属) | 点列フィットで追従可 | 中 (`ModeFWall` 差し替え + 積分幾何) | 単体では**低**、(b) 経路と併用で真価 | 弧長パラメータ化の数値整備; 検証系 (`validate`) 再構築 |
| D. CST 単一曲線 | スロート込み単一 $C^\infty$ 曲線 + 内部点線形拘束 | 全域 $C^\infty$ (接合消滅) | BPO 次第 (4–9) | 点列 LSQ で追従可 | 中 | 単体では**低**〜中 | κ の狙い撃ち不可; グローバル基底ゆえ局所修正が波及 (現行 Bernstein target と同じ病) |
| E. CONTUR 型接合廃止 | 可変壁始点をスロート直後へ、throat characteristic 初期線化 | 接合レス (壁は場の出力) | 0 | チェーン再編 (逆 MOC 初期線変更) | **大** | **高** (こぶ帯が設計到達域に入る) | 既存資産 (帰還マップ・taper・x_reach 系) の大改修; 検証やり直し |
| F. 目標基底の拡張 | 自由 CP 3→4–5 or B-spline target | (目標曲線の話) | +1–2 | 変更なし | 小 | こぶには無効 (§4 ジレンマ専用) | EHVI 評価回数増 (~1.5–2×/次元) |
| G. KL 高次遷音速解 | アンカー生成器の二次化 + start line M=1.1–1.3 | — | 0 | starting line 品質向上 | 中 (B6 残課題) | 中 (CFD アンカー導入済のため増分限定) | 実装検証の工数; CFD アンカーと役割重複 |

---

## 8. 推奨案と最小実験計画

**推奨: A (R=5 の単発 A/B) を判別実験として先行し、その結果で C+E (κ(s) 表現によるスロート下流の設計出力化) への投資可否を決める。B (G3 単体)・D (CST 単体) は非推奨** (波源が到達不能域にあるため空振りが濃厚)。F はジレンマ対策として独立に軽く実施可。

根拠: (i) 文献 4 系統 (Cuffel 下限 / CONTUR 実例・級数精度 / Sauer 摂動スケール / 自前 R スイープの $\propto1/\sqrt R$) が全て $R=2$ を低 R 帯と判定、(ii) A は表現・チェーン・dv 構成を一切変えず制約 §6 を満たす、(iii) A の結果自体が「こぶの $R$ 依存性」という物理情報になり、C/E/G の要否判定に使える。

### 最小実験計画 (A/B: R=2 vs R=5)

- **run**: `case/41.wind_tunnel_design/run_0037_r5_cold` (新連番、既存 run 再利用禁止)。現行チェーンで $R=5$ の仮壁を生成 (Sauer アンカー pass 1 相当 — CFD アンカーは R=2 場に凍結されているため使わない)、node Euler (`_config_euler_node`)、メッシュは node config で新規変換、品質 VERDICT 記録、cold 単段起動 (run_0031 と同一 CFD 条件)。対照 = `run_0031_baseline_cold`〜`run_0034_node_base` 系列 (R=2)。
- **測定** (全て既存ツール):
  1. こぶ振幅・位置: 意図非依存うねり指標 (軸 M の 5 次トレンド差、x/rt 0.55–2.5) — 対照は R=2 の 0.019–0.023 @1.43。
  2. $x_{reach,CFD}$ (壁面状態起点の C⁻ トレース) と、そこでの $(M,M',M'')$ (局所フィット感度込み)。
  3. 壁 $dM/dx$ の J 両側線形近似 (対照: 0.989→0.595 の 40% 跳び)。
  4. mdot 比・出口 $\varepsilon_M$・$|\theta|_{max}$ (`core_mode="traced"`)。
  5. `check_convergence.py` / `check_quasisteady.py` の VERDICT (報告に必須)。
- **判定基準**:
  - こぶ振幅が **~1/2 以下** ($\lesssim0.010$、$1/\sqrt{R}$ スケール比 0.63 より深い低下) → 円弧起因が支配と確定。**MOO に $R\ge5$ の下限制約** (または R を活かした探索) を正式採用し、C/E は保留。ジレンマ対策 F のみ実施。
  - こぶ振幅が **ほぼ不変** ($\ge0.015$) → 円弧の一次誤差以外 (starting line 壁端 1.1% 残差、縦線 starting line の構成問題) が支配 → **E (接合廃止・throat characteristic 化) を C の κ(s) 表現とセットで計画化** (新規 plan を起票)。G (KL) はその前段の安価な検証として実施。
  - 中間 → R 依存成分と構成依存成分の分離値が得られたとみなし、MOO での R 扱いと E の優先度をその比で判断。
- **併記事項**: R=5 では出口までの長さ・リップ位置が変わるため、出口指標は R=2 と直接比較せず「風洞要求 (<0.5%, <0.1–0.25°) に対するマージン」で評価する。

---

## 9. 本調査で判明した「今の理解が誤っている/曖昧だった点」

1. **「Sivells/CONTUR は曲率連続の非粘性コンタを要件にしている」は半分だけ正しい**。要件は事実だが、その保証手段は壁側の幾何拘束ではなく**軸分布の C2 整合 + 壁の MOC 一体生成**。`wall_modef.py` docstring がこの文言を「壁側 C2 クランプの根拠」として引くのはニュアンス違い (壁側 C2 は無意味ではないが、CONTUR の主張の写像ではない)。
2. **Sauer の $M''$ は「精度が悪い」のではなく「存在しない」**。軸上 $u$ は厳密に $x$ の 1 次 (式 11 は仮定)。$M''$ の 2.3 倍差は構造的欠落で、starting line を下流に置いても Sauer 内では改善しない。
3. **$R=2$ は文献上「小さい側の下限ちょうど」**。「妥当」とされる根拠 (Back–Massier–Gier) は壁静圧の一致であり、Cauchy データ用の $M',M''$ 精度は担保外。R スイープ (B2) の不整合 $\propto\bar u$ は Sauer の $u_{wall}=1/(4R)$ スケーリングそのもの。
4. **「風洞は $R_d$=5–10 が定石」の明文規則は見つからない**。裏付けは三角測量 (CONTUR 実例 5.5/6、$R<12$ 補正、音速点 2 階微分符号 10.5、Cuffel 下限 2) であり、「定石」として引用するときはこの 4 点を根拠に挙げるのが正確。
5. **$M=1.05$ 縦線 starting line + 壁足=可変壁始点という構成は古典に対応物がない**。古典は throat characteristic (CONTUR) か M=1.1–1.3 可変 M start line (Kliegel/Shelton)。
6. **mdot 比 0.994 は妥当と確認** (Sauer 積分 0.9934、Szaniszlo $R_N=2.29$ 実測と整合) — 疑う必要なし。
7. **B9→B10 の撤回 (「支配方程式切替による不可避な波」の否定) は文献とも整合**。Korte は残留 0.02 を「遷音速モデルと実流れの差」に帰属しており (=改善可能な設計データ問題)、「物理的に不可避」とは述べていない。同時に、**0.02 級はフル古典チェーンの実測相場**でもあり、そこから先を削るには円弧/遷音速整合そのものに手を入れる必要がある。
8. **Ogawa–Boyce の「円弧+Bézier, G1 接合」が成立するのは目的関数が推力 (積分量) だから**。軸 M 分布 (微分的量) を目的にする本ケースに同型パラメトライズの前例として引いてはいけない。
9. **接合特異点を円弧で丸めても波は消えず塗り広がる** (Genin–Stark: 丸め半径増で影響域 $L_i/R_{th}$ 1→1.5 に増加) — 「局所を滑らかにすれば波が消える」という直感への反例 (ただし G1 破れの丸めの話で、本ケースへの定量転用は不可)。

---

## 10. 主要書誌と入手可否

**ローカルにあり (精読済)**: Sivells AEDC-TR-78-63 (CONTUR)。Cuffel–Back–Massier AIAA J 7(7) 1969。Sauer NACA TM-1147。Szaniszlo NASA TN D-7848。Kulfan AIAA 2007-62。Korte AIAA 2000-0677。Korte NASA TP-3050 (※ノズル設計とは無関係と判明)。Ogawa–Boyce JPP 28(6) 2012。Genin–Stark DLR 2011。Arthur 1952 (無関係)。

**Web 公開 (未ダウンロード、入手可)**:

- Sivells AEDC-TR-78-63 の DTIC 版: apps.dtic.mil/sti/tr/pdf/ADA062944.pdf (ローカル PDF と同一内容)。F77/Python 移植: github.com/aldorona/contur, github.com/noahess/conturpy
- Vallabh & Skews 2017 (CSIR, open access): scielo.org.za (pid=S2309-89882017000100008)
- Shen et al. 2017 (CIRCLE 検証, open access): DOI 10.1177/1748301816665527
- Moreton 博士論文 (MVC): UCB CSD-93-732。Farin Class A: farinhansford.com/gerald/papers/classA.pdf
- Chartrand 2011 (正則化微分): DOI 10.5402/2011/164564。Knowles & Renka: EJDE Conf. Proc. 21

**有料/未入手 (優先度順)**:

1. **Shope, AEDC-TR-04-2 (2005)** "Design Optimization of Hypersonic Test Facility Nozzle Contours Using Splined Corrections" — 本ケースの帰還ループ直系先行例。DTIC 検索を推奨
2. Kliegel & Levine, AIAA J 7(7) 1969, 1375–1378 (DOI 10.2514/3.5355) — B6 残課題 (KL 実装) の一次資料
3. Hall, QJMAM 15(4) 1962, 487–508 (DOI 10.1093/qjmam/15.4.487)
4. Darwell & Badham, AIAA J 1(8) 1963 (DOI 10.2514/3.1965); Back & Cuffel, AIAA J 4(12) 1966
5. Sivells JSR 7(11) 1970 (DOI 10.2514/3.30160); Sivells AEDC-TR-56-11 (2D 曲率連続の原点, 1956)
6. Korte 系: AIAA 92-4009 (CAN-DO); Korte & Hodge JSR 32(4) 1995; Shope AIAA-2006-3665
7. Korakianitis CIRCLE 系: ASME J. Turbomach. 115(2) 1993 / 135(4) 2013; Goodhand & Miller ASME J. Turbomach. 133(2) 2011
8. Beckwith 系: NACA TN 2711 (1952, NTRS 公開のはず); Chen–Malik–Beckwith 1990 (NTRS 19910046726)

**見つからなかったもの (正直な報告)**: Rao の $R_d=0.382r_t$ の原典明記箇所。「AEDC-TR-71 系 Sivells ノズルレポート」(存在未確認 — 系譜は TR-56-11→JSR 1970→TR-75-118→TR-78-63 が正)。κ' を陽制御する CST 拡張。超音速ノズルで G3 効果を G2 と分離定量した研究。Pope & Goin の $R_d$ 推奨本文 (貸出制で未確認)。

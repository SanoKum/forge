# LES/DES における中心差分 (KEEP) スキームのスプリアス振動抑制 技術調査

<!-- ファイル名規約: convection-<short-slug>.md / 本文書は調査専用・コード変更なし -->

## メタ

- **area**: `convection`
- **status**: `draft`（**調査専用 / コード変更なし**）
- **related_docs**:
  - [`methods/convection/theory.md`](../../methods/convection/theory.md)
  - [`methods/convection/implementation.md`](../../methods/convection/implementation.md)
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md)
- **related_plans / surveys**:
  - [`turbulence-des-flux-survey.md`](turbulence-des-flux-survey.md) — **(a) RANS↔LES 中心↔風上ブレンド** (f_d 駆動 σ)。本調査が扱わない側面①
  - [`convection-pep-scheme-survey.md`](convection-pep-scheme-survey.md) — **(b) 接触面・多成分の圧力平衡保持 (PEP)**。本調査が扱わない側面②
  - [`turbulence-iddes-sst.md`](../../plans/active/turbulence-iddes-sst.md) — SST-DDES/IDDES 実装計画 (§4.8 で低散逸 flux を扱う)
  - [`discretization-median-dual.md`](../../plans/active/discretization-median-dual.md) — node/cell-vertex 両対応化 (本調査の結論が直結)
- **created**: `2026-06-29`
- **owner**: `CFD Dev`
- **調査実施**: deep-research ワークフロー（105 エージェント / 23 ソース取得 / 25 クレーム検証: 22 確認・3 棄却）

> **位置づけ**: 本文書は実装計画ではなく **技術調査レポート**。「KEEP のような中心差分 (運動エネルギー保存型 split-form) を LES/DES で安定に回すために、振動抑制として何をすべきか」を一次文献で裏取りしてまとめる。既存2サーベイ ((a) DES flux ブレンド / (b) PEP) が扱う側面は**意図的に範囲外**とし、その上位に立つ統括ノートとして「振動源の4分類」とそれぞれの対策を整理する。

---

## 0. エグゼクティブサマリ（結論先出し）

1. **KEEP の split convective form は「エイリアシング不安定」(振動源 #1) を散逸ゼロで抑える de-aliasing 機構である**が、その効力は**非粘性・非圧縮極限にスコープされ、圧縮性の線形安定は保証しない**。KEEP 単体で LES/DES が安定に回ると考えてはいけない (Kuya 2018 / Okumura-Kuya-Sawada 2023, 検証済)。
2. **圧縮性では運動エネルギー保存 (KEP) だけでは非線形安定に不十分**。一次保存量 (primary invariants) の大域保存も同時に要る。split form 族の中で **KGP (Kennedy-Gruber/Pirozzoli) 三重分割が全エネルギー方程式形式で最も頑健** (Coppola et al. 2019, 検証済)。forge は保存形 FV なので大域保存は満たすが、「どの split を使うか」は KGP が安全側。
3. **エントロピー保存 (EC) も圧力平衡保持 (PEP) も、中心スキームの局所線形安定を治さない**。EC two-point flux には固有の **anti-diffusion** 機構があり、密度波テストは CFL・時間積分法に依らず短時間で破綻。KEP+PEP 同時flux (Shima) でも正実部固有値が残る (Ranocha-Gassner 2021, 検証済)。→ **既存 PEP サーベイ (b) の成果は「接触面の圧力振動」を消すが「安定化」はしない**。両者を混同しないこと。
4. **したがって中心/split スキームには「明示的な散逸層」が必須**。これは2段構えにする:
   - **(i) de-aliasing 基盤の散逸**: KE 安定かつエントロピー整合な **matrix/entropy-variable 散逸** (Chandrashekar 2013, FV 構成・forge に移植可)。粗格子・通常の乱流で常時効かせる薄い散逸。
   - **(ii) 不連続用の局所衝撃捕獲**: エントロピー安定だけでは Gibbs 振動 (振動源 #3) は消えない。**低次 FV と高次作用素の convex blend** (Hennemann-Gassner の subcell FV/DG) か **localized artificial diffusivity (LAD)** をセンサ駆動で局所適用 (検証済: 「ES だけでは衝撃振動を抑えない・追加の shock-capturing が要る」)。
5. **forge への最重要含意**: 非構造で **KEEP 性と2次精度を保つには cell-vertex (= node-centered / median-dual) 離散化が必要** (Okumura-Kuya-Sawada 2023, 検証済)。**forge の median-dual 路線が KEEP の正しい置き場であり、pure cell-centered は KEEP に不利**。
6. **3つの棄却結論 (過強主張)**: 「split form は強安定/KEP を**保証**する」「split form は**散逸を足さずに**頑健性を上げる」は棄却 (0-3 / 1-2)。**split form 単体は安定保証ではない** → 上記 (i)(ii) の散逸層は省略不可。
7. **穴埋め (#2/#3/#4) の目玉**: (#4) **中心スキームの KE 保存は対流作用素の歪対称性と等価** (Verstappen-Veldman)、free-stream 保存は離散 GCL ($\sum_f A_f=0$)+同一再構成 Jacobian を要する (Kok/Sun-Deng)。**★ forge 直撃: フル FP32 でメトリック (体積・面積ベクトル・法線) を計算すると非 watertight 格子で free-stream が崩れ早期遷移するが、メトリックのみ FP64 で回復** (Karp et al. 2506.05150) — forge 既知の [[forge-freestream-nonorthogonal]] 弱点の**最小改修・最大効果の処方箋**。(#3) 衝撃捕獲は **directional LAD** (Kawai-Lele μ\*/β\*/κ\*, $f_\mathrm{sw}$=H(−∇·u)·Ducros, 係数 Cμ=0.002/Cβ=1.75/Cκ=0.01; スカラー LAD は高 AR で dt∝1/AR なので **directional 必須**)。(#2) KEEP は市松/エイリアシングを制御しないので LAD μ\* か非構造 differential filter を別途要す (具体は medium)。

---

## 1. 振動源の4分類 (本調査のフレーム)

中心スキームの「振動」は単一原因でなく、独立した発生源に分けて潰す。混ぜると過散逸で解像乱流を殺す。本調査は既存2サーベイが扱う ((a) RANS↔LES ブレンド・(b) PEP 接触面) を除く以下を扱う。

| # | 振動源 | メカニズム | 標準対策 (本調査の結論) | forge 現状 |
|---|--------|-----------|------------------------|-----------|
| **#1** | エイリアシング不安定 | 非線形積 (ρuu) が grid scale 超の波数を生成→折返し→KE/エントロピー非保存で発散 | **split/skew-symmetric 形** (KEEP/KGP)。ただし圧縮性では単独不十分→ **KE/entropy 整合の薄い散逸を併用** | KEEP 実装済 (`KEEP_SLAU_d`)。散逸は SLAU ブレンド |
| **#2** | 高波数エネルギー堆積 (2Δ odd-even) | 中心差分は最高波数に散逸ゼロ→2Δ 波が成長 | split-form の**暗黙 de-alias が主**。残れば**陽的低域通過フィルタ** (非構造は differential filter)。SVV/3-2則は構造格子前提 | 専用フィルタ無し。`duc` フロアの SLAU 残留散逸が代用 |
| **#3** | 不連続 (衝撃/接触) の Gibbs 振動 | 中心は不連続で必ずリンギング。**ES/split でも消えない** | **センサ駆動の局所散逸** (LAD: Kawai-Lele 等) or **subcell 低次FV blend** | `duc` blend で KEEP↔SLAU、Ducros (現状不活性) |
| **#4** | free-stream/GCL 非保存・odd-even デカップリング | 非構造/曲線座標で離散 metric 恒等式が破れ、低散逸化で桁落ちが顕在化 | **対称保存スキーム** (Kok 2009) ・離散 metric 整合・telescoping | 既知の非直交 free-stream 桁落ち ([[forge-freestream-nonorthogonal]])、`lowMachPrecond=2` で市松根治済 |

> **#1 と #3 を混同しない**: KEEP (split) は #1 のエイリアシング不安定を治すが、衝撃の Gibbs (#3) は治さない。逆に衝撃用散逸を全域にかけると #1 は安定だが解像乱流が死ぬ。**DES の本質は #3 の散逸を f_d で空間局在化し、LES 域は #1 の split + 薄い #1 散逸だけで安定化すること** (既存サーベイ (a) の結論と整合)。

---

## 2. #1 split-form / entropy-stable・EC 枠組みの比較 (検証済)

### 2.1 KEEP の機構と限界

- **機構** (Kuya-Totani-Kawai, *JCP* 375:823-853, 2018): 質量・運動量の対流項を **split convective form** に書き換え、離散化した質量・運動量式と「支配方程式間の解析的関係」が**全エネルギー流束を決定**する。これが運動エネ↔内部エネ交換を正しく解く鍵で、エントロピー保存にも効く。**追加散逸ゼロで頑健**。
- **限界 (検証時の釘刺し)**: この de-aliasing 効力は **非圧縮・非粘性極限にスコープ**され、**局所線形安定ではない**。「エントロピー保存に necessary」は原文の "important for" の言い過ぎ (EC flux は別経路で EC を達成)。

### 2.2 圧縮性では KEP だけでは安定でない — split 族の順位

- **Coppola, Capuano, Pirozzoli, de Luca (*JCP* 382:86-104, 2019)**: 「KEP だが一次変数を大域保存しない形式は非粘性で典型的に不安定」(非圧縮では KEP だけで十分なのと対照的)。反例 KG1/KG2 を非粘性 Taylor-Green で実証。
- **順位**: energy-preserving split は **2パラメータ族** (Feiereisen は δ=0, KGP は δ=1/4)。**KGP 三重分割が全エネルギー方程式形式で最も頑健**。Feiereisen と Coppola の C 形は「全エネルギー式 (全エンタルピー分割) かエントロピー式」とでのみ安定。
- **forge 含意**: forge は保存形 FV なので「一次変数の大域保存」は満たす (有利)。split の選択は **KGP 寄りが安全**。現 `KEEP_SLAU` の中心枝がどの分割形かを確認し、KGP 系でなければ移行を検討。
- **CAVEAT**: 構造格子・有限差分・倍精度の研究。**離散性質のヒューリスティックとして移植**でき、非構造 cell-centered FV / GPU / float32 での実証ではない。

### 2.3 EC / PEP は線形安定を治さない (Ranocha-Gassner)

- **Ranocha & Gassner (*CAMC* 2021, arXiv:2009.13139)** "Preventing Pressure Oscillations Does Not Fix Local Linear Stability...": 「PEP 性は EC flux の安定問題の治療ではない。単純な密度波で依然破綻し、CFL・時間積分法に鈍感」。「局所線形安定の欠如は **EC two-point flux に固有の anti-diffusion 機構**と強く関係」。「Shima の KEP+PEP flux のスペクトルも同様に**正実部固有値**を持つ」。
- **forge 含意**: **PEP サーベイ (b) の成果 (接触面圧力振動の除去) を「安定化」と読み替えてはいけない**。安定化は別途 §2.4 の散逸層が担う。CAVEAT: 高次 split-form DG での証明。**2次 median-dual FV で同じ anti-diffusion 強度が出るかは未確立** (→ open question)。

### 2.4 必要な散逸層の作り方 (Chandrashekar)

- **Chandrashekar (*Comm. Comput. Phys.* 14(5):1252-1286, 2013, arXiv:1209.4994)**: 中心 KEP flux は未決定自由度を持つが、**エントロピー整合を課すと全項が一意に決まる** (KEP+EC flux)。「中心 flux は衝撃や粗格子では散逸を足す必要がある。**KE 安定かつエントロピー条件を満たす scalar 人工散逸**を構成」。「**entropy-variable ベースの matrix (Roe 様) 散逸**は KE・エントロピー安定で、かつ original Roe と違いエントロピー違反解を生まない」。
- **forge 含意**: これは **FV 構成で密度ベース forge に直接移植可能**。現 `KEEP_SLAU` の散逸を「SLAU」から「entropy-variable matrix dissipation」に置き換える/補強するのが原理的に正しい方向。CAVEAT: 厳密 iff の運動量 flux 形 (`f_m=f_ρ·u_avg+p̃`) は **棄却 (1-2)** なので「算術平均速度形が KEP の必要十分」とは扱わない。

### 2.5 split/skew-symmetric の基盤 (Abe / Pirozzoli / Tadmor)

- **Abe et al. (*JCP* 353:193-227, 2018)**: split (skew-symmetric) 形の FR で、**primary conservation (PC) と KEP は別個の離散性質で別個の条件を要する**ことを多項式解析で厳密導出。
- **Pirozzoli (*JCP* 229:7180-7190, 2010)**: 任意次数の局所保存 split-convective 近似 (Ducros 2000 の一般化)。
- 基盤枠組み: **Fisher-Carpenter SBP + Tadmor two-point entropy-conservative flux**。entropy-stable DG (2024-2025) はこの上に乗る。
- **棄却 (重要)**: Pirozzoli の「split form が強安定を**保証**」(0-3)・「散逸なしで頑健性向上」(1-2) は過強。**split 単体は安定保証でない**。

---

## 3. #2 高波数エネルギー堆積 / フィルタ / SVV / de-aliasing — 一部検証 (穴埋め 2026-06-29)

> 「中心 split-form だけでは不十分=別途機構が要る」方向性は高信頼。ただし**非構造向けフィルタの具体的定式化・強度設定は medium 以下** (一次文献を fetch したが最終検証 25 件に残らず)。

- **検証済み (反証側)**: **KEEP の歪対称中心は KE を保存するが、それは高波数のエントロピー/エイリアシング堆積や odd-even 市松モードを制御しない** (split-form は FD でエイリアシング誤差をむしろ増やしうる: Reiss-Sesterhenn arXiv:1308.6672)。低散逸ほど市松が顕在化 (arXiv:2408.06821)。→ **これが forge が KEEP↔SLAU を `duc` ブレンドする理由そのもの**であり、#2 専用機構 (フィルタ/SVV/LAD μ\*) が要る根拠。
- **暗黙 de-alias が第一手**: split-form (KEEP/KGP) は de-aliasing に寄与するが**万能ではない** (上記)。**§4 の LAD 人工せん断粘性 μ\* が #2 を兼ねる**のが forge には現実的 (専用フィルタを足さずに高波数を減衰)。
- **非構造での陽的フィルタ (fetch 済・最終未検証)**: 構造格子の compact filter (Visbal-Gaitonde, Lele 6次)・SVV・3/2則は**構造格子/スペクトル前提**で非構造に直接乗らない。非構造の代替は **differential (Helmholtz 型) filter** — Bull & Jameson ([JBull *JCP* 2015](http://aero-comlab.stanford.edu/Papers/JBull_JCP_2015.pdf)) / Lund 系 (*JCP* 2015, S0021999115001898) が本命候補。ヘルムホルツ方程式を解く陰的フィルタなので GPU コストに注意。
- **forge 含意**: median-dual KEEP が回り始めたら、まず**陽的フィルタ無し**で 2Δ odd-even が出るかを Taylor-Green / backstep で確認。出れば **(優先) LAD μ\* → §2.4 matrix 散逸強化 → differential filter** の順で対処。SVV/compact filter は非構造 forge に非親和。

---

## 4. #3 不連続の Gibbs 振動 / localized artificial diffusivity (LAD) — 検証済 (穴埋め 2026-06-29)

> 穴埋め deep-research で一次文献を敵対的検証 (3-0 主)。LAD は forge の衝撃捕獲層の本命。

### 4.1 LAD の構成 (Cook-Cabot → Kawai-Lele 系)

中心スキームに**3+1種の局所人工拡散**をセンサ駆動で局所的に足す (Cook & Cabot 2004/2005 起源、Kawai & Lele *JCP* 227:9498-9526, 2008 が曲線・異方性格子へ、Kawai-Shankar-Lele *JCP* 229:1739-1762, 2010 が正準化):

- **人工せん断粘性 μ\***: 乱流の**未解像勾配**を減衰 (= #1/#2 の高波数堆積にも効く)
- **人工バルク粘性 β\***: **衝撃捕獲** (#3 本体)
- **人工熱拡散 κ\***: **接触不連続**捕獲
- **人工種拡散 D\***: 多成分での組成不連続 (forge の TP 多成分に該当)

### 4.2 係数の式構造 (検証済・正準値)

各係数は **〈4階微分 × グリッドスケール Δξ⁴ × センサ〉** (8階微分説は**棄却**、正しくは4階):

$$
\mu^* = C_\mu\,\rho\,\overline{\textstyle\sum_l \left|\frac{\partial^4 S}{\partial \xi_l^4}\right|\Delta\xi_l^4}\,\Delta^2,\quad
\beta^* = C_\beta\,\rho\,\overline{f_\mathrm{sw}\textstyle\sum_l \left|\frac{\partial^4 (\nabla\cdot u)}{\partial \xi_l^4}\right|\Delta\xi_l^4}\,\Delta^2,\quad
\kappa^* = C_\kappa\,\frac{\rho c_s}{T}\,\overline{\textstyle\sum_l \left|\frac{\partial^4 e}{\partial \xi_l^4}\right|\Delta\xi_l^4}\,\Delta
$$

- 定数: **$C_\mu=0.002$, $C_\beta=1.75$, $C_\kappa=0.01$** (Kawai-Shankar-Lele 2010)。$S$=ひずみ速度、$e$=内部エネルギー、$\Sigma_l$ は計算3方向、overbar は truncated-Gauss フィルタ。
- **$f_\mathrm{sw}$ (衝撃局在センサ) は $\beta^*$ にのみ付く**。$\mu^*/\kappa^*$ は4階微分+フィルタで**自己局在**する (滑らか域で4階微分が小)。

### 4.3 衝撃局在センサ $f_\mathrm{sw}$ (検証済)

$$
f_\mathrm{sw} = H(-\nabla\cdot u)\cdot\frac{(\nabla\cdot u)^2}{(\nabla\cdot u)^2 + (\nabla\times u)^2 + \epsilon},\qquad \epsilon\approx 10^{-32}
$$

- **負ダイレーション Heaviside × Ducros** の積。**圧縮 (衝撃) 域のみ発火、渦度支配の乱流域では作動しない** = 解像乱流を守る要。
- 正準形は $\min(\cdot,1)$ キャップと $\Omega=\max(|\nabla\times u|, 0.05c/\Delta)$ クリップを併用。Ducros 1999 起源、Hendrickson-Kartha-Candler 2018 が Heaviside を精緻化 (既存サーベイ (a) の「改良 Ducros」と同一系譜)。

### 4.4 異方性格子の罠と directional LAD (検証済・forge に直結)

- **スカラー LAD は高アスペクト比格子で過散逸**し、数値剛性で**安定 dt が 1/AR に律速**される。BL の高 AR セルを持つ forge には致命的。
- **Olson & Lele (*JCP* 246:207-220, 2013) の directional LAD** は係数を方向ごと ($\xi,\eta,\zeta$) に独立に加えて緩和。M=2 斜め SBLI で **dt 0.01→0.17 (14.3倍高速)、壁圧・Cf は同一**。
- **棄却 (1-2)**: 「LAD は強衝撃で reconstruction 法より急に発散する」は不採用 (LAD は強衝撃でも実用)。

### 4.5 乱流統計の保存 (検証済)

- LAD は中心高次と併用しても**解像 LES 乱流統計を保存**。M=1.3, $Re_\theta=1181$ 超音速 TBL で van Driest 平均・密度スケール Reynolds 応力が Pirozzoli DNS と整合 (細格子)。バッファ層の流れ方向 Reynolds 応力がわずかに過大 (粗格子で悪化)。
- ただし **LAD は微小な非零散逸を足す** (細格子 Cf は LAD off が良い) → **常時全域でなくセンサで局在**させるのが要。

### 4.6 forge 含意

- 衝撃用局所散逸は **「RANS↔LES の f_d ブレンド (a)」とも「#1 の de-alias 散逸」とも別レイヤ**で設計 (3レイヤ分離)。
- forge には Ducros センサ (`ducrosSensor_d`, 現状不活性) があり、$f_\mathrm{sw}$ はその系譜なので**再利用可能**。
- **directional LAD が必須** (forge は BL 高 AR セル前提)。非構造 median-dual で「方向別」をどう定義するか (面法線/エッジベース) は open question (§9)。4階微分の非構造 GPU 評価コストも要計測。
- 多成分なので $D^*$ (人工種拡散) も入れれば接触面の組成振動に効く (PEP サーベイ (b) と補完関係)。

---

## 5. #4 free-stream/GCL 保存・odd-even デカップリング — 検証済 (穴埋め 2026-06-29)

> **forge が最も必要とする分類**。穴埋め deep-research で敵対的検証 (3-0 主)。**forge 直撃の float32 知見が出た**。

### 5.1 KE 保存 ≡ 対流作用素の歪対称性 (Verstappen-Veldman)

- **Verstappen & Veldman (*JCP* 187:343-368, 2003)**: 中心スキームの KE 保存は離散対流作用素の**歪対称性 $C(u)+C(u)^\mathsf{T}=0$ と数学的に等価** (iff)。これを満たせば**任意格子で安定**かつ質量・運動量・KE (物理散逸 off 時) を保存。
- **罠**: 標準微分ステンシルは**非一様格子で係数行列に非零対角を生む**ため歪対称性を破る (実歪対称行列は対角=0 が必然)。→ KE 保存が壊れ、風上の暗黙減衰か陽的人工散逸が要り、**対流-物理散逸の微妙なバランス (乱流の本質) を損なう**。
- **高次精度と KE 保存は素朴には両立しない**: truncation error を物理空間で最小化しようとステンシル係数/補間をメッシュ依存にすると**離散積積分則 (discrete product rule) が破れ KE 保存が壊れる** (Kok; Morinishi/Vasilyev の「strict conservation か strict 4次の二択」)。→ **forge は「対称保存」を「高次化」より優先すべき**。
- CAVEAT: 原典は**非圧縮・スタッガード・構造格子**。collocated/非構造への一般化は **Trias et al. *JCP* 2014** が後続 (要精読)。

### 5.2 free-stream 保存の十分条件 (Kok / Sun-Deng)

- **Kok (*JCP* 228:6811-6832, 2009)** cell-centred 曲線座標 FV:
  - **(a) uniform-flow consistency (= free-stream 保存)** は**面積ベクトル和=0 の離散 GCL** $\sum_f A_f = 0$ を要する。体積/面積ベクトルを cell-vertex 座標から〈辺=直線・cell=trilinear 補間〉で計算すれば恒等式が成立。
  - **(b) 運動量と KE の両保存**は歪対称形=離散発散形の等価に依存し、それには**離散積積分則**が必要 = **メトリックを平均作用素に入れず、平均重みを 1/2 固定、メトリックは cell volume と面積ベクトルのみ経由**する設計でのみ成立。
  - 中心スキームは陽に加えない限り**数値散逸ゼロ** (散逸は SGS/粘性応力のみ)。
- **Sun/Deng et al. (*Comput. Fluids* 2016)**: free-stream 保存は FV (曲線座標) でも FD 同様に崩れ、崩れると重大誤差。十分条件: **(i) メトリック評価多項式の次数 ≤ 面流束 Gauss 求積の精度次数**、**(ii) 各 Gauss 点の Jacobian を流束/解と同一再構成アルゴリズムで評価** (さもないと離散メトリック項が相殺しない)。CAVEAT: 証明は2次元。
- CAVEAT: Kok の trilinear 構築は**六面体構造/曲線セル向け**。forge の**任意非構造/median-dual ポリヘドラ (可変 valence・非平面面) への直接適用は射程外** (§9)。

### 5.3 ★ float32 桁落ちと free-stream 喪失 (forge 直撃・検証済 3-0)

- **Karp et al. (arXiv:2506.05150, "Effects of lower floating-point precision on scale-resolving simulations of turbulence", §8.3)**: 高次圧縮性ソルバ SSDC を**フル FP32 でコンパイルすると free-stream 保存が崩れた**。FP32 でのメトリック幾何量 (**Jacobian・cell 体積・法線**) 計算が**非 watertight 格子**を生み、対流項経由で非物理アーティファクト・数値ノイズ増大・**早期乱流遷移**を起こした。
- **これらメトリック計算のみ FP64 に戻すと正常化**。
- **forge への含意 (決定的)**: forge の既知弱点 [[forge-freestream-nonorthogonal]] (float32 非直交 free-stream 崩壊) と**機構が完全一致**。**処方箋 = 混合精度: 本体は FP32 のまま、メトリック (体積・面積ベクトル・法線・wall_dist 等の幾何量) のみ FP64 で計算**。これは「全体を double にする」より遥かに安価で、過去の「double solve はフォールバック」議論 ([[mixed-precision-axisym-refuted]]) とも整合する (幾何だけ倍精度)。
- CAVEAT: コンパイル FP32 ビルド固有 (エミュレート FP32 は妥当に挙動)、SSDC スコープと原文が明記。ただし機構一致ゆえ含意は強い。

### 5.4 forge 含意のまとめ

- **LES 域で中心比率を上げる (=散逸を下げる) ほど free-stream 誤差が顕在化**するため、これは「中心化の前提条件」。
- 即効処方箋は **§5.3 のメトリック FP64 化** (最小改修・最大効果)。
- 構造的には **§5.1-5.2 の対称保存 + 離散 GCL ($\sum_f A_f=0$) + 同一再構成 Jacobian** を median-dual に移植 (Trias 2014 精読が前提)。
- odd-even/checkerboard は forge では `lowMachPrecond=2` で既に根治 ([[backstep-lowmach-checkerboard-precond2]]) だが、中心比率を上げると再燃しうるので #2 のフィルタ/LAD μ\* と合わせて監視。

---

## 6. forge への接ぎ木: 優先度付き推奨

検証済み結論から、「forge が KEEP を LES/DES で安定に回すために #1-#4 で何を足すか」を優先度順に。

| 優先度 | 施策 | 根拠 (検証済) | forge での具体 |
|--------|------|--------------|---------------|
| **P0** | **散逸を完全には消さない**。KEEP 単体は線形不安定なので、常時薄い KE/entropy 整合散逸を残す | Ranocha-Gassner (EC/PEP≠安定, anti-diffusion) / Coppola (KEP単独不十分) / split保証は棄却 | `KEEP_SLAU` の `duc` フロア (0.05) を**撤廃しない**。むしろ「LES 域でも 0 にしない」を原則化 |
| **P0** | **median-dual (node) を KEEP の主戦場にする** | Okumura-Kuya-Sawada 2023 (cell-vertex で KEEP 性+2次精度) | [[node-mode-periodic-and-backstep-status]] の median-dual を KEEP-LES の既定離散とする。pure cell-centered KEEP は KEEP 性が劣化する前提で扱う |
| **P0** | **★ メトリック幾何量 (体積・面積ベクトル・法線・wall_dist) のみ FP64 化** | Karp et al. 2506.05150 (フル FP32 メトリック→非 watertight→free-stream 喪失・早期遷移; メトリックのみ FP64 で回復, 検証 3-0) | forge 既知の [[forge-freestream-nonorthogonal]] と機構一致。**全体 double 化より遥かに安価な即効処方箋**。幾何セットアップを FP64 で計算し FP32 へ格納/転送。中心化を進める前提 |
| **P1** | **基盤散逸を entropy-variable matrix dissipation 化** | Chandrashekar 2013 (KE+entropy 安定, Roe のエントロピー違反を回避, FV 構成) | `KEEP_SLAU` の散逸枝を SLAU から/に加えて entropy-consistent matrix 散逸へ。defect-correction で陰解法と両立 ([[implicit-blockdplur-config]]) |
| **P1** | **3レイヤを分離設計**: ①RANS↔LES f_d ブレンド (既存(a)) / ②#1 de-alias 薄散逸 / ③#3 衝撃局所散逸 | ES だけでは衝撃振動を抑えない (3-0) | f_d (既存) ・matrix 散逸フロア・Ducros×Heaviside の LAD を**独立**に重ねる。limiter を衝撃センサ流用しない (既存(a)の結論) |
| **P1** | **衝撃捕獲を directional LAD で** (Kawai-Lele μ\*/β\*/κ\*+D\*) | Kawai-Shankar-Lele 2010 / Olson-Lele 2013 directional (検証 3-0); スカラー LAD は AR で dt∝1/AR | $\beta^*$=Cβ(1.75) + $f_\mathrm{sw}$=H(−∇·u)Ducros、$\mu^*$(Cμ=0.002) は #2 高波数も兼ねる。**directional 必須** (forge は BL 高 AR)。既存 `ducrosSensor_d` を流用。多成分は $D^*$ も |
| **P2** | **split は KGP 系を採る/確認** | Coppola 2019 (KGP 最頑健) | 現 `KEEP_SLAU` 中心枝の分割形を確認、非KGPなら移行検討。forge は保存形なので大域保存は OK |
| **P2** | **対称保存 + 離散 GCL ($\sum_f A_f=0$) を median-dual に** | Verstappen-Veldman 2003 / Kok 2009 / Sun-Deng 2016 (検証 3-0) | メトリックを平均作用素に入れず重み1/2固定・面積ベクトル和=0・同一再構成 Jacobian。非構造一般化は Trias 2014 精読が前提 |
| **P3** | **#2 陽的フィルタは三の矢** | KEEP は市松/エイリアシングを制御しない (検証); LAD μ\* が先 | LAD μ\* で足りなければ非構造 differential filter (Bull-Jameson 2015)。SVV/compact は非親和 |

**一言で**: forge の正しい道は「**median-dual 上の KEEP (KGP) + メトリック FP64 + 常時薄い entropy-consistent matrix 散逸 (de-alias) + f_d で局在化した directional LAD (衝撃/高波数)**」。KEEP も PEP も**安定化の代替にならない**ことを設計前提にする。**最小改修で最大効果は P0 のメトリック FP64 化** (既知の発散弱点に直結)。

---

## 7. 検証票

2回の deep-research。**第1回 (#1 中心, 25 claims: 22 確認/3 棄却)** + **第2回 穴埋め (#2/#3/#4, 25 claims: 23 確認/2 棄却)**。

**第1回 確認 (3-0 主)**:
- KEEP の split-form 機構と「解析的関係が energy flux を決定」(Kuya 2018 / Okumura-Kuya-Sawada 2023)
- 圧縮性で KEP 単独 ≠ 非線形安定、一次変数大域保存も要、KGP 最頑健 (Coppola 2019)
- EC/PEP ≠ 局所線形安定、EC two-point flux の anti-diffusion、Shima KEP+PEP も正実部固有値 (Ranocha-Gassner 2021)
- entropy 整合で flux 一意化 (KEP+EC)、中心は散逸要、entropy-variable matrix 散逸が KE+entropy 安定 (Chandrashekar 2013)
- ES/split だけでは衝撃 Gibbs を抑えない、convex 低次FV/高次 blend が要 (arXiv 2504.00173 / Hennemann-Gassner)
- 高次中心は本質的にエイリアシング不安定で安定化技術を要する、split が de-alias 手段 (Abe 2018 / Pirozzoli 2010 / Tadmor-Fisher-Carpenter)
- **非構造 KEEP-FV は cell-vertex (median-dual) が KEEP 性+2次精度に必要** (Okumura-Kuya-Sawada 2023)

**第2回 穴埋め 確認 (3-0 主)**:
- LAD μ\*/β\*/κ\*(+D\*) の式構造〈4階微分×Δξ⁴×センサ〉と正準係数 Cμ=0.002/Cβ=1.75/Cκ=0.01 (Kawai-Shankar-Lele 2010)
- β\* センサ $f_\mathrm{sw}$=H(−∇·u)·Ducros、衝撃のみ発火 (Springer chapter / arXiv 2511.18042)
- スカラー LAD は AR で dt∝1/AR、directional LAD で M=2 SBLI 14.3倍高速・Cf 同一 (Olson-Lele 2013)
- LAD は中心高次でも解像乱流統計を保存 (DNS Pirozzoli 整合)
- KE 保存 ≡ 対流作用素の歪対称性、非一様格子で標準ステンシルが破る (Verstappen-Veldman 2003)
- free-stream 保存は離散 GCL ($\sum_f A_f=0$)・離散積積分則・同一再構成 Jacobian を要する (Kok 2009 / Sun-Deng 2016)
- **★ フル FP32 メトリック→非 watertight→free-stream 喪失・早期遷移、メトリックのみ FP64 で回復 (Karp et al. 2506.05150)**

**棄却 (両回)**:
- 「運動量 flux `f_m=f_ρ·u_avg+p̃` **iff** KEP」(1-2) / 「split に中心差分で強安定を**保証**」(0-3) / 「split は**散逸なし**で頑健化」(1-2)
- 「LAD は8階微分ベース」(誤、正しくは4階) / 「LAD は強衝撃で reconstruction より急発散」(1-2 不採用)

**信頼度の留保**: #2 (非構造フィルタの具体的定式化・強度) は **medium 以下** — 「中心 split だけでは不十分」方向は高信頼だが、Bull-Jameson differential filter 等の一次文献は fetch したが最終検証 25 件に残らず。Verstappen/Kok は**非圧縮スタッガード/六面体構造**が原典で非構造一般化は外挿 (Trias 2014 要精読)。Karp は SSDC スコープ (機構一致ゆえ含意は強い)。

---

## 8. 一次文献 (検証済ソース)

| # | 文献 | 役割 |
|---|------|------|
| 1 | Kuya, Totani & Kawai, *JCP* 375:823-853 (2018) [S0021999118305916](https://www.sciencedirect.com/science/article/abs/pii/S0021999118305916) | KEEP 原典 (split convective form) |
| 2 | **Okumura, Kuya & Sawada, *JCP* 494:112521 (2023)** [10.1016/j.jcp.2023.112521](https://dl.acm.org/doi/10.1016/j.jcp.2023.112521) | **非構造 KEEP-FV (cell-vertex/median-dual 必要)・forge 直結** |
| 3 | Coppola, Capuano, Pirozzoli & de Luca, *JCP* 382:86-104 (2019) [S0021999119300245](https://www.sciencedirect.com/science/article/abs/pii/S0021999119300245) | split 族 2 パラメータ化・KGP 最頑健・KEP単独不十分 |
| 4 | Ranocha & Gassner, *CAMC* (2021) [arXiv:2009.13139](https://arxiv.org/pdf/2009.13139) | **EC/PEP ≠ 線形安定・anti-diffusion** |
| 5 | Chandrashekar, *CiCP* 14(5):1252-1286 (2013) [arXiv:1209.4994](https://arxiv.org/pdf/1209.4994) | KEP+EC flux・entropy-variable matrix 散逸 (FV) |
| 6 | Abe et al., *JCP* 353:193-227 (2018) [S0021999117307453](https://www.sciencedirect.com/science/article/pii/S0021999117307453) | PC と KEP は別個条件・skew-symmetric FR |
| 7 | Pirozzoli, *JCP* 229:7180-7190 (2010) [ResearchGate](https://www.researchgate.net/profile/Sergio-Pirozzoli/publication/222569549_Generalized_conservative_approximations_of_split_convective_derivative_operators) | 任意次数 split-convective (Ducros 一般化) |
| 8 | entropy-stable collocation DG, *Comp.&Fluids* (2025) [arXiv:2504.00173](https://arxiv.org/pdf/2504.00173) | ES+subcell shock-capturing が必要 |
| 9 | Glaubitz/Ranocha (2024) [arXiv:2406.14557](https://arxiv.org/pdf/2406.14557) | ES baseline は anti-alias 安定化を judicious に要する |
| 10 | **Kawai & Lele, *JCP* 227:9498-9526 (2008)** [S0021999108003641](https://www.sciencedirect.com/science/article/abs/pii/S0021999108003641) | LAD 曲線・異方性格子版 (μ\*/β\*/κ\*) |
| 11 | Kawai, Shankar & Lele, *JCP* 229:1739-1762 (2010) [S0021999109006160](https://www.sciencedirect.com/science/article/abs/pii/S0021999109006160) | LAD 正準係数・乱流統計評価 |
| 12 | **Olson & Lele, *JCP* 246:207-220 (2013)** [S0021999113002040](https://www.sciencedirect.com/science/article/abs/pii/S0021999113002040) | **directional LAD (高 AR で dt 緩和)** |
| 13 | Kumar & Vadlamani (ISTAM PA0234) [Springer 10.1007/978-981-97-0418-7_9](https://link.springer.com/chapter/10.1007/978-981-97-0418-7_9) | LAD 式・係数・センサの逐語ソース |
| 14 | Cook & Cabot, *JCP* (2004) [S0021999104004000](https://www.sciencedirect.com/science/article/abs/pii/S0021999104004000) | LAD 原型 (人工バルク粘性 hyperviscosity) |
| 15 | **Verstappen & Veldman, *JCP* 187:343-368 (2003)** [PDF](https://pure.rug.nl/ws/files/3005906/2003JCompPhysVerstappen.pdf) | KE 保存 ≡ 歪対称性 (対称保存) |
| 16 | **Kok, *JCP* 228:6811-6832 (2009)** [NLR PDF](https://reports.nlr.nl/server/api/core/bitstreams/2652c350-e8a1-472b-afb6-d26ddbcde12a/content) | 曲線 FV の free-stream/KE 両保存・離散 GCL |
| 17 | Sun/Deng et al., *Comput. Fluids* (2016) [S0045793016300184](https://www.sciencedirect.com/science/article/abs/pii/S0045793016300184) | FV free-stream 保存の十分条件 (metric/Jacobian) |
| 18 | **Karp et al. (2025)** [arXiv:2506.05150](https://arxiv.org/pdf/2506.05150) | **★ FP32 メトリック→free-stream 喪失、FP64 メトリックで回復** |
| 19 | Bull & Jameson, *JCP* (2015) [PDF](http://aero-comlab.stanford.edu/Papers/JBull_JCP_2015.pdf) | 非構造 differential filter (#2, fetch 済・最終未検証) |
| 20 | Reiss & Sesterhenn (2014) [arXiv:1308.6672](https://arxiv.org/pdf/1308.6672) | split-form が FD でエイリアシングを増やしうる (#2) |

---

## 9. open questions (穴埋め後に残る = 実装前に詰める)

1. **#4 非構造一般化**: Verstappen-Veldman/Kok の対称保存・離散 GCL ($\sum_f A_f=0$・メトリックを平均作用素に入れない・重み1/2固定) を、forge の **median-dual ポリヘドラ (可変 valence・非平面面)** へどう移植するか。Kok の trilinear 構築が使えない場合の離散メトリック恒等式の満たし方 → **Trias et al. *JCP* 2014 (collocated 一般化) の精読が次の一手**。
2. **#4 混合精度の実証**: Karp の「メトリックのみ FP64・本体 FP32」が forge の median-dual FV で有効か。符号付き和/Kahan 和・面積ベクトル telescoping でどこまで桁落ちを緩和できるか。**過去の double-solve 議論 ([[mixed-precision-axisym-refuted]]) を「幾何だけ FP64」で再評価**。
3. **#3 directional LAD の非構造定義**: Olson-Lele の方向別 ($\xi,\eta,\zeta$) を、構造化方向を持たない median-dual で**どう定義するか** (面法線方向ベース? エッジベース?)。4階微分作用素の非構造 GPU 評価コスト・精度、JST との優劣。
4. **#2 非構造フィルタの要否**: KEEP の暗黙 de-alias + LAD μ\* で 2Δ odd-even が足りるか、Bull-Jameson differential filter (ヘルムホルツ型・陰的・GPU コスト) が要るか。**Taylor-Green / backstep で実測**。
5. **Ranocha-Gassner の anti-diffusion は 2次 median-dual FV で出るか**: 高次 DG 現象か、forge の実離散でも効くか。**median-dual KEEP が実際にどれだけの (どんな形の) 散逸を要するか**を実測。既存 `KEEP_SLAU` ブレンドで足りるかの実証が起点。

---

## 変更ログ

- `2026-06-29` — 初稿 (調査専用)。第1回 deep-research (105 エージェント / 23 ソース / 25 クレーム検証: 22 確認・3 棄却) + forge コンテキストを統合。既存2サーベイ ((a) DES flux / (b) PEP) の上位に立つ統括として「振動源4分類」を定義。核心結論: KEEP/split は de-alias を与えるが圧縮性の安定は保証しない (Coppola/Ranocha-Gassner)、entropy-consistent matrix 散逸 (Chandrashekar) + 局所衝撃捕獲の2段散逸が必須、非構造では cell-vertex/median-dual が KEEP の正しい置き場 (Okumura-Kuya-Sawada)。
- `2026-06-29` — **穴埋め**: 第2回 deep-research (102 エージェント / 20 ソース / 25 クレーム検証: 23 確認・2 棄却) で #2/#3/#4 を grounding。§4 (LAD: μ\*/β\*/κ\* 式・係数・$f_\mathrm{sw}$・directional LAD)・§5 (free-stream: 歪対称性=KE保存・離散 GCL・**★Karp の FP32 メトリック→free-stream 喪失/FP64 で回復**)・§3 (#2 は KEEP 限界の反証側で grounding、フィルタ具体は medium) を確定。§6 推奨に **P0「メトリックのみ FP64 化」** (forge 既知発散弱点の即効処方箋) と **P1「directional LAD」** を追加。残る open question は非構造一般化 (Trias 2014)・混合精度実証・directional の非構造定義・#2 要否・anti-diffusion の 2次 FV 実測。

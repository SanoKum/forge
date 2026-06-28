# 極超音速燃焼風洞・飛翔体における凝縮液滴/雨滴の挙動 — 技術調査

<!-- ファイル名規約: <short-slug>.md (調査メモ) -->

## メタ

- **area**: `condensation` / `multiphase` (横断)
- **status**: `draft` (**調査専用 / コード変更なし**)
- **related_plans**:
  - [`condensation-nonequilibrium.md`](../../plans/active/condensation-nonequilibrium.md) — 非平衡凝縮モデルの実装計画 (本調査のテーマ1と直結)
- **related_docs**:
  - [`methods/thermophysics.md`](../../methods/thermophysics.md) — 多成分 TP gas 基盤
- **created**: `2026-06-28`
- **owner**: (未定)
- **手法**: deep-research ハーネス 2 回 (前段=素過程, 後段=end-to-end 深掘り)。各回 5 アングル × 22 ソース取得 × ~90 主張抽出 → 上位 25 を 3 票方式で敵対的検証。

> **位置づけ**: 本文書は実装計画ではなく **技術調査レポート**。極超音速燃焼風洞での凝縮液滴の挙動 (テーマ1) と、降雨中飛翔時の雨滴の燃焼器到達 (テーマ2) について、一次文献・支配無次元数・数値手法を棚卸しする。forge への展開を判断する場合は [`condensation-nonequilibrium.md`](../../plans/active/condensation-nonequilibrium.md) 等の plan 側で別途設計すること。

### 凡例 (引用の検証ステータス)

| 記号 | 意味 |
| --- | --- |
| ◎ | 本調査の敵対的検証 (3 票方式) で 3-0 確証 |
| ○ | 検証バッチで 2-0 / 2-1 など過半数で確証 (一部留保) |
| △ | 情報源は取得したが今回の検証バッチで独立確証されていない (要追検証) |
| ✗ | 検証で棄却された主張 (誤りとして本文に記録) |

---

## 0. エグゼクティブサマリ (結論先出し)

- **テーマ1 (風洞内凝縮)**: ノズルの急膨張で **N2 (空気主成分)・H2O・燃焼生成物が非平衡均質凝縮**を起こす。平衡飽和線を**大きく超える過冷却** (supercooling) を経て核生成・液滴成長が進み、生成物は**潜熱放出で試験気体の熱力学状態 (設計流れ条件) を変える** ◎。供試体まわりの衝撃層での蒸発は**滞留時間 (=供試体/プローブ寸法・衝撃離脱距離) 依存**で、薄い衝撃層では蒸発が間に合わず凝縮が残存する ◎。
- **テーマ1 の運命 (蒸発 vs 成長)**: 衝撃通過直後の水滴は昇温により **「凝縮成長」ではなく「蒸発」側**に推移する ◎。ただし弱衝撃 ($M_{sf,0}=1.3$) では蒸発による径減少はわずか (~0.2 μm)、強衝撃でブレイクアップが起きると最終径は初期径によらず **~4 μm に収束** ◎。
- **テーマ1 の汚染**: 燃焼加熱型 (vitiated) 風洞の **H2O/CO2 気相汚染**は着火遅れ・燃焼特性・推力を実質的に変え、scramjet 試験データを汚染する ◎。FAA が規制するほど water ingestion は重大 ◎。
- **テーマ2 (降雨中飛翔)**: 雨滴は弓状/斜め衝撃波通過後に **Weber 数** (遷移 We は **Ohnesorge 数**の関数) で支配される空力ブレイクアップ (bag〜sheet/wave-crest stripping) を起こす ◎。鈍頭体への到達/運動量伝達は **「ブレイクアップ時間 vs 機体到達までの利用可能時間の競合」**で決まる ◎。
- **数値手法**: 支配的なのは**一方向結合 Euler-Lagrange** + **KHRT/TAB** ブレイクアップ + 変形考慮の経験的抗力。最新は **Eulerian-Lagrangian (EL)** でブレイクアップ+蒸発を同時解析 (2D double-wedge / 3D cone-cylinder-flare, TAMU) ◎。
- **🔴 最大の空白 (あなたの2つの問いの核心)**: 「凝縮液滴/雨滴が**実機 scramjet 供試体の衝撃波系を通過して燃焼器に到達するか**」の **end-to-end 定量評価**を直接扱った一次文献は、2 回の深掘りでも確立できなかった。素過程 (凝縮・蒸発・ブレイクアップ) は堅固だが、それらを燃焼器入口までつないだ統合研究は **open** (§5)。
- **🔴 高マッハ域の限界**: shock-droplet 相互作用の確立データは概ね **M≲5-6**。真の極超音速 (M≫6) は TAMU 等が **Mach 5-30 の normal shock** を数値で埋めつつあるが、実験検証は乏しく、古典相関 (Ranger & Nicholls) の高マッハ外挿は**棄却** ✗。

---

## テーマ1: 風洞内で凝縮した液滴の供試体まわりの挙動

### 1-1. ノズル膨張での非平衡凝縮の物理

**空気主成分 (N2) の均質凝縮** ◎ — 極超音速ノズルの急膨張で純窒素 N2 が均質凝縮する。M=10/14/18 の純窒素流で凝縮の発現・消失が実験観測され (過冷却 ~22-25 K、凝縮質量分率 12-14%)、半経験モデル (簡易化 DFT でナノスケール N2 液滴をモデル化、自由エネルギー障壁・表面張力・臨界液滴密度を核生成・成長モデルに供給) が確立。軸対称 Mach10 ノズルの凝縮開始静温で公表値と**最大 3.5 K の偏差**で一致。
〈出典〉NASA NTRS [19900050882](https://ntrs.nasa.gov/citations/19900050882) / AIAA SciTech 2020-0381 [Homogeneous Condensation of Nitrogen in Hypersonic Wind Tunnels](https://www.researchgate.net/publication/338400653_Homogeneous_Condensation_of_Nitrogen_in_Hypersonic_Wind_Tunnels_A_Semi-Empirical_Model)

**非平衡過冷却が本質** ◎ — 蒸気相は平衡飽和線を**超えて大幅に過冷却**してから核生成する (平衡飽和点では始まらない)。過冷却の到達度は**ノズル形状とシーディング (汚染物) に依存**。contoured ノズルや低膨張率 (膨張率パラメータ $\dot{P}<1000$) では過冷却が著しく制限され**不均質核生成が支配的**になる。均質 (spontaneous) 核生成理論は $\dot{P}>1000$ で実験と一致、小 $\dot{P}$ で乖離。
〈出典〉NTRS 19900050882 / VKI Longshot [AIAA 2014-1153 (Grossir)](https://dipot.ulb.ac.be/dspace/bitstream/2013/208925/3/AIAA_2014_1153_Grossir.pdf)

**水蒸気の凝縮衝撃波** ◎ — H2O の非平衡凝縮はノズルのスロート〜末広部で**凝縮衝撃波 (condensation shock)** として顕在化。Schlieren 高速可視化で凝縮衝撃波・弱擾乱波の位置・角度から発現位置を実測できる。
〈出典〉Buttsworth et al. 2020, impulse tunnel TUSQ (N2/H2O 混合気) [ScienceDirect S030193222030584X](https://www.sciencedirect.com/science/article/abs/pii/S030193222030584X)

**理論的基盤 (古典核生成理論 CNT)** ◎ — Wyslouzil & Wölk のレビューが基盤文献。等温核生成速度の標準測定法 (膨張チャンバ、nucleation-pulse、超音速ノズル) と CNT 理論対実験の比較 (Wölk-Strey 経験補正含む) を扱う。
〈出典〉[J. Chem. Phys. 145, 211702 (2016)](https://pubs.aip.org/aip/jcp/article/145/21/211702/196311/Overview-Homogeneous-nucleation-from-the-vapor)

> **歴史的基準文献**: **Wegener & Mack「Condensation in Supersonic and Hypersonic Wind Tunnels」(1958, Advances in Applied Mechanics)** が古典レビュー (検索で最重要一次文献として浮上、paywall のため claim 抽出には未使用 △)。本テーマの出発点。

### 1-2. 燃焼加熱型 (combustion-heated) 風洞での凝縮

炭化水素燃焼を熱源とする燃焼加熱型風洞では、**燃焼生成物の水蒸気がノズル膨張で凝縮**して二相流を生じ、設計流れ条件 (ノズル出口の熱力学環境) を変える ◎。相変化と圧縮性流れの連成のため見積もりは困難で、縮約モデルは 2 段階に分ける: (1) 凝縮なしの実在気体ノズル流で**過飽和**出口状態を得る、(2) ノズル末端で過飽和→飽和へ**不連続ジャンプ**する。
〈出典〉Lin et al., Shock Waves 28:321-333, 2018 [10.1007/s00193-017-0724-x](https://link.springer.com/article/10.1007/s00193-017-0724-x)

### 1-3. 衝撃波背後で蒸発するか残存するか (供試体まわり)

凝縮物 (液滴/氷粒子) が衝撃波背後で蒸発するか残存するかは、**衝撃層内の滞留時間 (=供試体/プローブ寸法・衝撃離脱距離) に依存**する ◎:

- **大型プローブ**は離脱距離が大きく高温衝撃層での**滞留時間が長い → 凝縮物が蒸発**し、ピトー圧に影響しない (25.4 mm 半球プローブのピトー圧/淀み点熱流束は凝縮に完全に不感)。
- **細径プローブ**は薄い衝撃層で**蒸発が間に合わず凝縮を検出**できる。
- 凝縮は潜熱放出で静温を上げ、密度ほぼ不変ゆえ静圧も同程度上がる → **静温・静圧が凝縮検出に好適**。

〈出典〉VKI Longshot AIAA 2014-1153 (M=14, N2)、McBride & Sherman 参照

**衝撃通過後の運命 = 蒸発側** ◎ (深掘りで確認) — 衝撃通過直後の水滴は「凝縮成長」ではなく「蒸発」側に推移する。弱衝撃 ($M_{sf,0}=1.3$) では $d_{d,0}=$30/50 μm はブレイクアップせず、蒸発による径減少は ~0.2 μm にとどまる。$d_{d,0}=$70/90 μm では衝撃通過直後にブレイクアップし、最終径は初期径によらず **~4 μm に収束**。
〈出典〉Huang/Zhu/Davy [arXiv 2103.10576](https://arxiv.org/pdf/2103.10576) (逐語確認、ただし 1 次元衝撃管・弱衝撃・TAB モデル)

**極低温風洞 N2 凝縮の翼まわり影響** ◎ — 局所ガス温度が飽和線を越えかつ臨界 (過冷却/Wilson) 点より低温になって初めて、流れ特性 (圧力・温度・マッハ数) を顕著に変える量の凝縮物が生じる。凝縮影響は**圧力係数 $C_p$・抗力係数 $C_d$** として解析。3D 化 (ONERA M6 翼、NASA Langley TCT 条件) も Euler-Euler+古典核生成+液滴成長モデルで実施。
〈出典〉Cryogenics 2020 (NACA0012-64) [S0011227520301673](https://www.sciencedirect.com/science/article/abs/pii/S0011227520301673) / Eng. Appl. CFD 2024 (ONERA M6) [10.1080/19942060.2024.2387055](https://www.tandfonline.com/doi/full/10.1080/19942060.2024.2387055)

> ⚠️ **レジーム注記 (重要)**: これら **cryogenic wind tunnel は遷音速 (M≈0.8 級) であって極超音速ではない**。ONERA M6 翼は遷音速 CFD の標準検証ケース ($M\approx0.84$)。クライオ風洞の目的は高 Reynolds 数確保のためにガス (N2) を ~110-150 K に冷やすことで、凝縮の駆動は「**あらかじめ飽和すれすれに冷やした気体を、翼上面の弱い局所加速で飽和線を割らせる**」もの。極超音速風洞の「**大膨張で深い過冷却を作る**」とは**運転レジームも凝縮トリガも別物**である。したがって本 2 文献は **3D N2 凝縮を CNT+液滴成長+Euler-Euler で解く「数値手法リファレンス」としては有効**だが、**極超音速供試体まわり凝縮の直接証拠ではない**。M≈0.8 では強い衝撃層がなく、§1-3 冒頭の「衝撃層滞留時間に依存した蒸発 vs 残存」の論点はそのままは当てはまらない (関連は弱い遷音速衝撃のみ)。

> **限界 (negative finding)**: これら N2 凝縮研究は影響を**翼まわりの空力 ($C_p$, $C_d$)** として評価しており、凝縮物が**供試体内部 (インテーク〜燃焼器) へ残留して燃焼・点火・PIV/圧力/熱流束計測を汚染する end-to-end 経路は一次文献で確立できなかった** (§5 open)。

➡️ **あなたの問い (テーマ1) への含意**: 方向性としては「**衝撃背後の昇温で蒸発側**」「**大スケール供試体の衝撃層では蒸発が進みやすい**」。逆に薄い斜め衝撃層・狭隘低温域では残存し得る。ただし**燃焼器内残留・計測汚染の直接データは未確認**であり、定量判定は今後の課題。

### 1-4. vitiation 汚染としての試験データ汚染

地上 scramjet 試験で流入空気を水素/炭化水素燃焼で予熱すると、実飛行に存在しない **H2O・CO2 汚染**が加わり、燃焼過程と燃料の着火特性を実質的に変える ◎:
- AIAA 2021-4177: 水素燃焼予熱は実飛行に無い水蒸気汚染を加え燃焼に substantial effect (7.4%/18.3% H2O 流入で推力 4.3%/15.2% 低下)。
- Liang et al. (Acta Mech Sinica 2014): vitiated air heater 由来の H2O/CO2 汚染が RP-3 ケロシンの**着火遅れ**を地上試験と実飛行で変える (衝撃波管実測、800-1500 K, 0.3 MPa, H2O 10-20%/CO2 10%)。

〈出典〉[AIAA 2021-4177](https://arc.aiaa.org/doi/10.2514/6.2021-4177) / [Acta Mech Sinica 10.1007/s10409-014-0014-0](https://link.springer.com/article/10.1007/s10409-014-0014-0)

> **注意**: vitiation は**気相 H2O 汚染**であり、液滴取り込み (テーマ2) とは物理が別。切り分けが必要。

---

## テーマ2: 降雨中飛翔時の雨滴が燃焼器に到達するか

### 2-1. 対象の径スケールとリスク

雨滴・雲滴・燃料噴霧は高速機体に対しエロージョンおよび water ingestion のリスク ◎。径範囲: **雲滴 1-100 μm、雨滴 100-10000 μm、scramjet 燃料噴霧 10-200 μm**。極超音速飛翔体は上層雲の sub-micron 〜 海面近くの mm 雨滴に遭遇する ◎。
〈出典〉AIAA 2021-0751 [Evaluation of Droplet Aerodynamic Breakup Models](https://www.researchgate.net/publication/348240253_Evaluation_of_Droplet_Aerodynamic_Breakup_Models_in_Supersonic_and_Hypersonic_Flows) / FTaC 2024 [10.1007/s10494-024-00581-z](https://link.springer.com/article/10.1007/s10494-024-00581-z) / [TAMU FMECL](https://fmecl.engr.tamu.edu/research/hypersonic-droplet-breakup/)

### 2-2. 空力ブレイクアップの支配則 (Weber / Ohnesorge)

ブレイクアップは **Weber 数** (空力 vs 界面張力) が支配し、レジームは**遷移 Weber 数**で区分される。その**遷移 We は Ohnesorge 数の関数** (粘性が高いほど遷移 We が増大) ◎。モードは bag / bag-and-stamen / sheet stripping / wave crest stripping。入射衝撃 Mach 数が上がる (We 大) ほど**ブレイクアップ完了距離と最終液滴径がともに減少** (径20-80 μm、$M_0$=1.3-4.0、初期 We 10.0-4758.3 で網羅)。
〈出典〉AIAA 2021-0751 / Huang/Zhu/Davy FTaC 2024

**ブレイクアップ機構 (多段階連成)** ◎ — 液滴扁平化 → 表面不安定性・穿孔の形成 → 変形一次液滴のリガメントからの二次液滴放出。駆動は強い加速・せん断・表面張力に起因する流体力学的不安定性 (**Rayleigh-Taylor / Kelvin-Helmholtz**)。Weber 数で RTP mode (We<100) / SIE mode (We>1000) に大別。Ranger & Nicholls の無次元時間 $\tau$ で規格化すると Mach 2 と Mach 3 でブレイクアップ過程はほぼ同速度 (ただし Mach 2-3 限定、M≫6 一般化は未確立)。
〈出典〉Ullman, Bielawski & Raman, PRF [arXiv 2504.09007](https://arxiv.org/abs/2504.09007) (逐語確認) / TAMU FMECL

### 2-3. 鈍頭体への到達を決める一次機構

鈍頭体先端への液滴衝突・運動量伝達は、**ブレイクアップ時間と、液滴が弓状衝撃波通過後に機体へ到達するまでの利用可能時間との競合**で支配され、スケーリング則でよく特徴づけられる ◎。これが「雨滴がインテーク方向にどれだけ慣性を保って到達するか」を決める一次機構。
〈出典〉Briney & Balachandar, Phys. Fluids 35, 016103, 2023 (Mach 2/3/6 鈍頭体, one-way Euler-Lagrange + stochastic TAB) [リンク](https://pubs.aip.org/aip/pof/article/35/1/016103/2868366/Euler-Lagrange-stochastic-modeling-of-droplet)

> 🔴 **棄却された直感** ✗: 「弓状衝撃波で液滴がブレイクアップ→微小化→慣性低下で機体衝突を回避」という単純描像は **3 票棄権で棄却 (支持なし)**。正しくは上記の「ブレイクアップ時間 vs 到達時間の競合」が支配機構。large な雨滴はブレイクアップ完了前に機体/インテークに到達し得る。

### 2-4. 真の極超音速 (M≫6) の shock-droplet 相互作用

水滴の変形・ブレイクアップは **Mach 5 まではよく研究**されているが、それを超える高マッハ域のデータは乏しく、これが現在の研究動機 ◎。TAMU 等は **Mach 5-30 の normal shock** に対する液滴ブレイクアップを数値シミュレーションで取得し、極超音速域 (M≫6) を埋めつつある。弓状衝撃波構造は複雑で、液滴は衝撃波・膨張波の系を通過し非定常な加速履歴を経験する (供試体/インテーク前縁での運命に直結)。
〈出典〉AIAA 2025-1500/1501 [2025-1500](https://arc.aiaa.org/doi/abs/10.2514/6.2025-1500) / [TAMU FMECL](https://fmecl.engr.tamu.edu/research/hypersonic-droplet-breakup/)

> ⚠️ **留意**: M5-30 は **normal-shock マッハ数** (自由流飛行マッハ数ではない) かつ**数値のみ・予備的 (early results)**。古典相関 (**Ranger & Nicholls**) の Mach 30 までの外挿妥当性は **0-3 で棄却** ✗。

### 2-5. water ingestion の燃焼への定量影響

**規制になるほど重大** ◎ — エンジンコアへの water ingestion は熱力学サイクルの変化で圧縮機サージ・推力損失・フレームアウトを引き起こす。機序は「燃焼器内での液体水の蒸発→火炎温度低下→化学反応速度低下→完全燃焼阻害→燃焼効率・安定性の低下」で、燃焼器は **sub-idle 運転時に最も失火しやすい**。FAA は飛行安全上の脅威と認め **Advisory Circular 33.78-1** を発行。
〈出典〉[FAA AC 33.78-1](https://www.faa.gov/documentLibrary/media/Advisory_Circular/AC_33_78-1.pdf) / [Federal Register 00-3702](https://www.federalregister.gov/documents/2000/02/16/00-3702/advisory-circular-turbine-engine-power-loss-and-instability-in-extreme-conditions-of-rain-and-hail)

**水/燃料比 vs 燃焼温度の定量** △ — WFR=0.5 で液体水を噴射すると燃焼器出口 (火炎) 温度が、メタン+空気で ~39-43℃、ディーゼル+空気で 39℃ 低下。ディーゼルエマルジョン予混合では 52.8℃、メタン+蒸気予混合では 19.1℃。
〈出典〉MSc thesis [DiVA 1793877](https://www.diva-portal.org/smash/get/diva2:1793877/FULLTEXT01.pdf) (STAR-CCM+/Cantera、査読論文より証拠強度弱い)

**液膜形成** ◎ — 高い水取り込み条件下では IPS (inlet particle separator) に入った雨滴は不可避的に壁面に**水膜 (water film)** を形成する。超音速/極超音速インテークへの water ingestion 経路で stream tube への水分布・液膜形成が起こることの傍証。
〈出典〉[SAGE 09544100231162425](https://journals.sagepub.com/doi/10.1177/09544100231162425)

> **注意**: water ingestion の燃焼定量影響の一次根拠は主に**航空ターボエンジン (FAA, IPS) とガスタービン/内燃の水噴射**であり、scramjet/ramjet 燃焼器そのものの定量データではない。機序は確立だが scramjet への外挿は推論を含む。臨界 water/air 比などの定量閾値はケース依存 (「水 5% で熱効率 41.2% 低下」の主張は **1-2 で否決** ✗)。

---

## §4. 数値手法・支配無次元数の整理 (forge への含意)

| 対象 | 支配無次元数 | モデル/相関 | 数値手法 |
|---|---|---|---|
| ノズル非平衡凝縮 | 過冷却度・過飽和度、Knudsen | 古典核生成理論 CNT (自由エネルギー障壁)、Hertz-Knudsen 液滴成長、Wölk-Strey 補正 | Euler-Euler / Euler-Lagrange-Euler (気相-液滴-液膜)、2 段階縮約 (過飽和→飽和ジャンプ) |
| 液滴ブレイクアップ | **Weber** (遷移 We = f(**Ohnesorge**))、Mach、$\tau$ (Ranger-Nicholls 無次元時間) | KHRT、TAB (stochastic)、変形考慮の経験的抗力 | 一方向結合 Euler-Lagrange 粒子追跡 |
| 衝撃-液滴/蒸発の同時解析 | We, Oh, Re, Stokes、蒸発 Damköhler | RTP (We<100) / SIE (We>1000)、RT/KH 不安定性 | **Eulerian-Lagrangian (EL)** (2D double-wedge / 3D cone-cylinder-flare) |

**数値手法の到達点** ◎ — 超音速/極超音速の液滴問題で支配的なのは**一方向結合 Euler-Lagrange**。ブレイクアップは KHRT / stochastic TAB、抗力は液滴変形を取り込む経験的抗力。KHRT+変形抗力は物理的に妥当だが、**高 Mach・高 We 域の較正には追加実験データが必要** (系統的検証が不足)。最新は EL でブレイクアップ+蒸発を**同時解析** (TAMU)。
〈出典〉AIAA 2021-0751 / Phys. Fluids 35, 016103 / AIAA 2025-1501

主要数値手法ソース: [Modified Euler-Lagrange-Euler for condensing droplets/films (Int. J. Heat Mass Transfer)](https://www.sciencedirect.com/science/article/pii/S0017931022010067)

---

## §5. 未解決の問い (あなたの2つの問いへの正直な現状)

2 回の深掘りでも、**素過程は堅固だが燃焼器入口までの end-to-end 評価は文献から合成できなかった**:

1. **【テーマ1+2 の核心 / end-to-end】** 風洞凝縮物または降雨/SLD 雨滴が、**実機 scramjet 供試体の弓状衝撃→斜め衝撃→インテーク圧縮場**を通過する一連の経路で、ブレイクアップ+空力加熱蒸発が**燃焼器到達『前』に完了するか?** 蒸発完了距離 vs 実機インテーク長さの定量比較を行った一次文献は確立できなかった。canonical 形状の EL 解析と 1 次元衝撃管はあるが、実機インテーク統合評価は **open**。
2. **【テーマ1 後半 / 計測汚染】** 風洞 test-section の凝縮物 (N2/H2O/CO2 氷粒子・液滴) が供試体内部 (インテーク〜燃焼室) に残留して燃焼・点火・PIV/圧力/熱流束計測を**実際に汚染した『事例』**を報告した一次文献は見つからなかった。
3. **【テーマ2 / 燃焼定量】** scramjet/ramjet 燃焼器『そのもの』への water ingestion が失火・着火限界・推力・燃焼効率に与える定量データ (航空ターボエンジン/内燃機関ではなく超音速燃焼器固有の water/air 比 vs 性能低下相関)。機序は確立だが scramjet 固有の定量閾値は未確立。
4. **【高マッハ外挿】** $\tau$ 規格化でブレイクアップ速度が衝撃強度に依らずほぼ一定という知見は **Mach 2-3 のみ**確認。真の極超音速域 (M≫6、高温・高 We・蒸発 Damköhler 数が支配的になり得る域) でこの規格化や古典相関が成立するかは実験検証付きで未確立。
5. **【極超音速の供試体/翼まわり3D 凝縮】** 3D の N2 凝縮 CFD は存在するが**遷音速クライオ風洞 (ONERA M6, M≈0.8) のもの** (§1-3 レジーム注記)。**極超音速 (M≫6) の供試体/翼まわりで凝縮物の3D 場・蒸発残存を直接扱った一次研究は確立できなかった**。手法 (CNT+成長+Euler-Euler) は流用可だがレジーム外挿が必要。

---

## §6. 棄却された主張 (誤りとして記録)

| vote | 棄却された主張 |
|---|---|
| 0-3 ✗ | 弓状衝撃波で液滴がブレイクアップ→微小化→慣性低下で機体衝突を回避するのが雨滴到達の一次機構 (正しくは「ブレイクアップ時間 vs 到達時間の競合」) |
| 0-3 ✗ | 極低温風洞で注入 N2 液滴が全温 200 K 以下では蒸発せず流れ場に残留して粒子追跡可能 (全温閾値以下残留の直接実証は不成立) |
| 0-3 ✗ | Ranger & Nicholls モデルが Mach 30 までの極超音速域でブレイクアップ時間を妥当推定 (高マッハ外挿は未確立) |
| 1-2 ✗ | jet-A1 への水エマルジョンで水 5% 超で熱効率が約 41.2% 低下 (定量閾値はケース依存、確証されず) |
| 0-3 ✗ | 燃焼器出口温度が WFR 0→1 で単調低下する一般的定量関係 (条件依存で一般化不可) |

---

## §7. 主要一次文献リスト

**テーマ1 (凝縮)**
- NASA NTRS 19900050882 — M=10/14/18 純 N2 凝縮の実験+理論 ◎
- AIAA SciTech 2020-0381 — N2 均質凝縮 半経験モデル (簡易 DFT) ◎
- ScienceDirect S030193222030584X — Buttsworth et al. 2020, 水蒸気凝縮衝撃波 (TUSQ) ◎
- J. Chem. Phys. 145, 211702 (2016) — Wyslouzil & Wölk, CNT レビュー ◎
- Shock Waves 28:321-333 (2018) — Lin et al., 燃焼加熱風洞の凝縮縮約モデル ◎
- VKI Longshot AIAA 2014-1153 (Grossir) — 凝縮検出・プローブ滞留時間依存 ◎
- Cryogenics 2020 (S0011227520301673) / Eng. Appl. CFD 2024 — N2 凝縮 翼まわり影響 ◎ (**遷音速クライオ風洞 M≈0.8、極超音速ではない**。手法リファレンスとして有効、§1-3 レジーム注記参照)
- AIAA 2021-4177 / Acta Mech Sinica 10.1007/s10409-014-0014-0 — vitiation 汚染 ◎
- (歴史的基準) Wegener & Mack 1958 △

**テーマ2 (液滴ブレイクアップ/ingestion)**
- Phys. Fluids 35, 016103 (2023) — Briney & Balachandar, Euler-Lagrange stochastic ◎
- AIAA 2021-0751 — Droplet aerodynamic breakup models 評価 ◎
- FTaC 2024 10.1007/s10494-024-00581-z (arXiv 2103.10576) — Huang/Zhu/Davy, shock-droplet ◎
- arXiv 2504.09007 (PRF) — Ullman, Bielawski & Raman, ブレイクアップ多段階機構 ◎
- AIAA 2025-1500 / 2025-1501 + TAMU FMECL — 極超音速 (M5-30) EL ブレイクアップ+蒸発 ◎
- Int. J. Heat Mass Transfer S0017931022010067 — Euler-Lagrange-Euler 凝縮 △
- FAA AC 33.78-1 / Federal Register 00-3702 — rain/hail エンジン失火規制 ◎
- DiVA 1793877 — WFR vs 燃焼温度 (MSc thesis) △
- SAGE 09544100231162425 — IPS 水膜形成 ◎

---

## §8. forge への含意・次アクション候補

- **テーマ1 は [`condensation-nonequilibrium.md`](../../plans/active/condensation-nonequilibrium.md) と直結**: 非平衡凝縮 (CNT + Hertz-Knudsen 成長) は本調査で確立された標準路線であり、forge の NASA 多項式 TP gas と整合させやすい。Euler-Lagrange-Euler (気相-液滴-液膜) が数値手法の主流。
- **「蒸発 vs 成長」の判定**: forge で供試体まわりの凝縮物の運命を扱うなら、衝撃層滞留時間と蒸発時間スケール (Hertz-Knudsen) の比をローカルに評価する後処理が定性判定に有効。
- **end-to-end は研究空白 = 新規性の余地**: §5 の 4 問は一次文献が薄く、forge (多成分 TP + 非平衡凝縮 + 将来の Lagrangian 粒子) で実機インテーク統合解析を行えば学術的寄与になり得る。着手時は plan を起こすこと。
- **高マッハは外挿に注意**: M≫6 のブレイクアップ/蒸発モデルは較正データ不足。古典相関の流用は検証付きで。

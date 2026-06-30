# ピントル/ニードル弁式 超音速ノズル推力調整機構 — 文献・特許まとめ

> **出典・関連物の所在**: 本書は調査メモ本体 (`notes/investigations/`)。**出典 PDF・生成スクリプト `build_pptx.py`・図 `.figs/`・生ログ JSON (`research_papers_raw.json` / `research_patents_raw.json`) は [`papers/pintle_nozzle/`](../../papers/pintle_nozzle/) にある**。以下の本文で `` `*.pdf` `` 等を同名参照している箇所は、このディレクトリ基準で読むこと。成果物 `pintle_nozzle_literature_review.pptx` は本ファイルと同じ `notes/investigations/` に置いている。

飛翔体のサイドスラスタ／DACS (Divert and Attitude Control System) で用いられる、超音速ノズルのスロート (throat / choke point) にバルブ（**pintle / needle / plug**）を突き刺して**等価スロート面積 (equivalent throat area)** を変え推力を調整する機構に関する調査メモ。

- 呼称: **pintle nozzle** / **needle-valve thrust modulation** / **throttleable solid rocket motor** / **plug-controlled supersonic nozzle**。「チョークバルブ」は学術用語ではない。
- 調査日: 2026-06-21
- 観点: (1) 空力性能、(2) 振動・動特性 (FIV/FSI)、(3) 固体粒子 (Al2O3) エロージョン
- 調査手法: deep-research ハーネス（fan-out web 検索 → 一次出典 fetch → 敵対的3票検証 → 統合）。論文調査・特許調査の2本を実施。
- 検証元の生データ: 下記「## 調査の生ログ」参照。

> **重要な留意**: 多くの ScienceDirect/AIAA/AIP 原典が HTTP 403 で全文取得不可のため、検証は複数独立検索による抄録の逐語突合に依存（逐語一致は確認済み）。本ディレクトリに PDF があるのは入手済みの4本のみ。残りは DOI/URL で参照。

---

## 技術調査報告 PowerPoint（成果物）

- **`pintle_nozzle_literature_review.pptx`** — 収録 PDF 5 件を全文精読して作成した文献レビュー（16:9, 全11枚）。
  - 構成: 表紙 / 全体まとめ2枚（成熟度マップ・重要知見と研究不足領域）/ 各文献1枚×5 / 横断比較3枚。
  - 各文献スライドは「背景・目的／方法／主要結果・考察／結論／本調査における知見・示唆」を色分けし、右に代表図（出典・Fig番号・頁を明記）を配置。著者主張（本文）と調査者評価（知見＝橙色）を区別。
  - 生成スクリプト: **`build_pptx.py`**（python-pptx で pptx を生成、同一内容を PIL でQAレンダリングしオーバーフロー検出）。図は `.figs/` の切り出し図（pptx に埋め込み済みで自己完結）。

## このディレクトリの PDF（入手済み）

| ファイル | 対応文献 | 観点 |
| --- | --- | --- |
| `numerical study of the dynamic characteristics of pintle nozzles for variable thrust.pdf` | **Heo, Jeong & Sung, *J. Propulsion & Power* 31(1), 2015 (DOI 10.2514/1.B35257)** — 主要文献#1。スライディングメッシュ、動特性・ヒステリシス | (1)(2) |
| `Pan_Chen_Ye_2017_pintle_tip_shape_specific_impulse_AtlantisPress.pdf` | **Pan, Chen, Ye, Atlantis Press vol.74 (2017)** — 主要文献#7。先端形状5種、密度ベースRANS、Isp効率（本調査で新規DL、OA） | (1) |
| `Transient flow characteristics and performance of a ... solid rocket motor wwith a pintle valve.pdf` | 移動ピントルSRMの過渡流れ・性能 (等価スロート面積／動メッシュCFD) | (1)(2) |
| `modeling and dynamic characteristics analysis on solid attitude control motor using pintle thrusters.pdf` | ピントルスラスタ式ACM の5サブシステム結合モデル（cf. *Aerospace Sci. Tech.* 106, 2020） | (1)(2) |
| `numerical investigation of the shock interaction effect on the lateral jet controlled missile.pdf` | 側方ジェット制御ミサイルの衝撃干渉（DACS隣接・側噴流） | (1) |

未入手の主要文献は「## 主要文献リスト」に識別子付きで掲載。**Elsevier(ScienceDirect)/AIP/MDPI は Cloudflare/Akamai の bot 保護で curl・WebFetch とも HTTP 403/Access Denied となり自動DL不可**（機関アクセスや手動DLが必要）。本調査で自動取得できたのは bot 保護のない Atlantis Press(Pan 2017) のみ。

---

## エグゼクティブサマリ

| 領域 | 成熟度 | 状況 |
| --- | --- | --- |
| 空力（推力制御則・衝撃干渉・Isp損失） | ★★★ | 複数の査読一次文献で unanimous。CFD手法（動メッシュ2D軸対称RANS）も確立 |
| 動的荷重・圧力ヒステリシス | ★★☆ | ピントル固有・実験検証あり。ただし準静的荷重で構造振動ではない |
| **FIV/自励振動・FSI** | ★☆☆ | **ピントル固有の直接研究は未確認**。制御弁からの類推のみ |
| **ピントル先端の粒子エロージョン** | ★☆☆ | 固定ノズルのモデルはあるが**ピントル可動部固有は空白** |
| 特許（2010以降） | ★☆☆ | 確証は ADD 二重ピントル(2019)のみ。中国系・Raytheon・三菱重工は確証に至らず |

**密度ベースCFDソルバ開発者としての寄与余地**:
1. ピントル先端衝撃波の**非定常FSI/自励振動**（高解像度非定常ソルバ＝LES/DESの出番）
2. **二相流（ガス-粒子）× ピントル可動 × 衝撃干渉 × FSI の統合シミュレーション**（推力調整サイクル中の摩耗・荷重・性能の同時評価）

---

## 1. 空力性能

### 支配機構: 等価スロート面積制御
ピントル位置・形状は**等価スロート面積**の変化を通じて燃焼室圧 $P_c$・推進薬燃焼速度・推力を制御。固体モータは内部弾道 $P_c \propto (A_b/A_t)^{1/(1-n)}$（$n$=圧力指数 <1）を介す。
- ピントル**挿入** → スロート縮小 → $P_c$・燃焼速度・推力**増**／**引抜** → 逆
- Sun et al., *Chinese J. Aeronautics* 33(12), 2020 ／ *Aerospace Sci. Tech.* 106, 106130 (2020)。［high, 3-0］

### CFD手法は確立
**ダイナミック／スライディングメッシュ**を用いた2D軸対称RANS（RNG k-ε等）が、cold-flow検証で**燃焼室圧を相対誤差2%未満**で予測。物理スロート面積の解析予測が数値結果とよく一致。
- Sun et al. 2020 ／ Heo, Jeong & Sung, *J. Propulsion & Power* 31(1), 2015 (DOI 10.2514/1.B35257)。［high, 3-0］

### ピントル先端衝撃波と Isp 劣化（CFD的に最重要）
ピントル挿入で**先端に必ず衝撃波**が生じ、`lip shock` / `trailing shock` / `oblique shock` / `shock train` を形成。**流れの運動量を散逸させ Isp を一般に劣化**。構造はストローク・先端形状に依存。
- Lee, Park & Yoon, *Aerospace Sci. Tech.* 26(1), 2012/13（Roe FDS、k-ε/k-ω SST/SA 3モデル比較＋実験）。基礎必読。［high, 3-0, 4件］
- 先端形状5種比較: 標準94.81%効率/365.52 N に対しピントル付加で**約1.5pt低下**、ピントル中は **small-arc 最良**(93.30%/361.47 N)。Pan, Chen, Ye, Atlantis Press vol.74 (2017)、密度ベースRANS。［medium — 単一未検証CFD］

### 剥離・衝撃干渉と内挿型の弱点
**内挿型 (internal)** はストローク増大に対し**推力調整比 <2%** なのに**空力荷重比 ~22%**。ストローク **4–5 mm** で衝撃波がノズル壁に当たり**剥離（peeling/separation bubble）**。
- *J. Korean Soc. Propulsion Eng.* (internal pintle thruster)。［high, 3-0］構成固有で一般化不可。

---

## 2. 振動・動特性 — ピントル固有のFIV研究は手薄

### 確認できた事実: 圧力ヒステリシスと動的荷重
ピントル運動が $P_c$ 変動 → **ヒステリシス** → アクチュエータ荷重急増。設計条件次第で**定常時の約2.5倍**。挿入・引抜・往復で応答遅れと感度が異なる。
- Ha & Kim, *Acta Astronautica* 201, 2022（cold-flow＋実射検証）／ AIAA 2015。［high, 3-0］
- **準静的な動的荷重であり構造振動（自励振動）ではない**点に注意。

### システムモデル化
ピントル式ACM/DACS = 電磁・ピントル運動・スロットル・チャンバ圧充填・推力の**5サブシステム結合系**（非線形状態方程式）。
- *Aerospace Sci. Tech.* 106, 106130 (2020)。［high, 3-0］

### FIV/FSI は「類推」依存
スロットル域の**剥離・再循環・渦放出が弁体FIVの主因**という非線形FSI解析はあるが、**原典は原子力用パイロット式制御弁**でピントルノズルではない。
- Yang et al., *Actuators* 14(8):372, 2025。［medium, 2-1 — 類推扱い］

> **空白**: ピントル/ニードル自体の流体励起・自励振動（ロックイン、リミットサイクル、FSI連成不安定、バフェッティング周波数）を**超音速ノズル流れ場で直接扱った査読研究は未確認**。＝寄与余地。

---

## 3. 固体粒子エロージョン — モデルはあるがピントル固有は空白

### 確認できたエロージョンモデル
**アルミ等金属粒子**によるSRMノズル壁エロージョンは、**強度理論・エネルギー保存則・Hertz接触理論**で導いたエロージョン率モデルで定量化。**衝突角と粒子濃度が支配因子**。モデル予測は実験と一致。
- Zhou, Yin, Sun et al., *AIP Advances* 12(12):125106, 2022（Rocstar、JPLノズル）。［high, 3-0］

> **留意**: 対象は**固定 C-D ノズル**でピントル固有ではない。**ピントル先端・ニードル表面への粒子衝突エロージョンを直接定量した研究は未確認**。可動部・先端よどみ点・衝撃干渉域での局所衝突角/濃度集中というピントル固有の摩耗は未開拓。
> 「二相流が convergent部のエロージョン・熱流束を増す」「Al含有率38%で機械的エロージョンが熱化学エロージョンと同オーダー」等は票割れ/未確証で除外（仮説としては妥当、要追加検証）。

---

## 4. 特許動向（振動抑制・エロージョン抑制・空力性能向上）

主要出願人: **Northrop Grumman（旧TRW／旧Alliant Techsystems）**、**Aerojet Rocketdyne（旧Aerojet General／旧Atlantic Research）**、**韓国ADD**、**日本IHI Aerospace**。

### (1) 振動抑制
| 特許 | 出願人 | 要点 |
| --- | --- | --- |
| **GB2283537A / US5456425A** (1993-95) | Aerojet General | 多ピントル**比例制御**。bang-bang PWM を連続比例制御に置換し機体ジッタ除去。**総スロート面積一定**で $P_c$ 平衡化、推力応答をfill/dump動特性から分離 |
| **US3948042A** (1968出願) | US Navy | 閉ループサーボに位置FBで **"stiffness"項**を導入し振動減衰。**制御系振動であり空力的自励振動ではない** |

### (2) エロージョン抑制
| 特許 | 出願人 | 要点 |
| --- | --- | --- |
| **US6330793B1** (1999) | Atlantic Research/Aerojet | **幾何形状**で対処。中間域を頂点状に隆起させ**燃焼粒子をスロートを逸れる軌道に乗せる**。固定ノズルだが思想を応用可 |
| **US3712063A** (1970) | Bell Aerospace | **Columbium合金芯＋Grafoil＋Carbitex積層**＋washer溝の**fuel-richガスフィルム冷却**。冷却ピントルの基礎 |

### (3) 空力性能向上
| 特許 | 出願人 | 要点 |
| --- | --- | --- |
| **US7849695B1 / US8215097B2** (priority 2001, granted 2010) | Alliant Techsystems/Northrop Grumman | **load-balanced pintle**。切頭二角錐ヘッド(60°/52.5°)で連続可変スロート、ブリードで荷重相殺、作動荷重 **約900→100 lbf (~89%減)** |
| **US10738740B1** (priority 2019) | **韓国ADD** | **内外二重ピントル**を単一共通ギア駆動、内外で異なるギア比で**流量レンジ拡大**。**確証された唯一の post-2010 特許** |
| **US8016211B2** (priority 2007) | Aerojet | **外部リングアクチュエータ**で駆動系をホットガス外に配置 |
| **US7565797B2 / US6591603B2** (2001-04) | GHKN／旧TRW | 可動/固定プラグの**高度補償で推力係数最大化** |
| **JP2002227723A** (2001) | **IHI Aerospace** | ノズル全体を機軸方向に往復動させ固定ピントルを出し入れ |

---

## 主要文献リスト（一次・査読、優先入手推奨）

識別子は本調査で実際に到達した一次URL/DOIに基づく（PII = ScienceDirect の article ID）。手で推測したDOIは誤りを避けるため記載しない。

| # | 文献 | 識別子 | 本ディレクトリ |
| --- | --- | --- | --- |
| 1 | Heo, Jeong & Sung, *J. Propulsion & Power* 31(1), 2015 — スライディングメッシュ、動特性 | DOI 10.2514/1.B35257 | **✓ 有** (`numerical study of the dynamic characteristics...pdf`) |
| 2 | Sun et al., *Chinese J. Aeronautics* 33(12), 2020 — 等価スロート面積、2%精度検証 | ScienceDirect PII S1000936120302260 | 未（bot保護で403） |
| 3 | Lee, Park & Yoon, *Aerospace Sci. Tech.* 26(1), 2012/13 — 衝撃構造・3乱流モデル比較 | ScienceDirect PII S1270963812000879 | 未（bot保護で403） |
| 4 | Ha & Kim, *Acta Astronautica* 201, 2022 — 圧力ヒステリシス・動的荷重2.5倍 | ScienceDirect PII S0094576522004969 | 未（bot保護で403） |
| 5 | Zhou et al., *AIP Advances* 12(12):125106, 2022 — 粒子エロージョンモデル | pubs.aip.org/aip/adv/article/12/12/125106 | 未（Cloudflareで403） |
| 6 | *Aerospace Sci. Tech.* 106, 106130 (2020) — ACM 5サブシステム結合・needle valve | ScienceDirect PII S1270963820308129 | 未（bot保護で403／類似PDF有） |
| 7 | Pan, Chen, Ye, Atlantis Press vol.74 (2017) — 先端形状5種、密度ベースRANS | atlantis-press article 25880156 | **✓ 有** (本調査でDL) |
| 8 | Yang et al., *Actuators* 14(8):372, 2025 — 制御弁FSI（類推） | DOI 10.3390/act14080372 | 未（Akamaiで403） |

### 特許（Google Patents/USPTO）
GB2283537A, US5456425A, US3948042A, US6330793B1, US3712063A, US7849695B1, US8215097B2, US10738740B1, US8016211B2, US7565797B2, US6591603B2, JP2002227723A

---

## 残された問い（追加調査候補）

1. ピントル/ニードル自体の流体励起・自励振動（ロックイン、リミットサイクル、FSI連成不安定）を超音速ノズル流れ場で直接扱った査読研究は存在するか。
2. Al2O3 等のピントル先端・ニードル表面への衝突エロージョン（局所衝突角・濃度分布、可動部固有摩耗）を直接定量したCFD/実験はあるか。
3. LES/DES でピントル先端衝撃波の非定常性・剥離バブル（圧力変動スペクトル、バフェッティング周波数）を捉えた最新（2020以降）研究の有無と密度ベースソルバでの再現性。
4. 二相流 × ピントル可動 × 衝撃干渉 × FSI を同時連成した統合シミュレーション事例。
5. 中国（航天科技/航天科工系, CN系）・Raytheon・Lockheed Martin・三菱重工・防衛装備庁の2010年以降の特許（本調査で未確証）。

## 調査の生ログ

検証済み主張・反証・スコア等の完全な生データ:
- 論文調査: `research_papers_raw.json`
- 特許調査: `research_patents_raw.json`

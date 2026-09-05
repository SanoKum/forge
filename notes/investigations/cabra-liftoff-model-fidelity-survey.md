# Cabra H₂/N₂ 浮き上がり火炎 — 文献でのリフトオフ再現水準とモデル要件 (2026-09-05)

forge (case/48) の反応 ON 計算が実験 H/d≈10 に対し**リップ付着**となり PaSR `C_mix` にも非感応である件について、「どのくらい合うべきなのか」「合わせるにはどのモデルが要るのか」を文献調査した。結論: **平均反応速度 (laminar-chemistry) クロージャの付着火炎は文献で立証済みの想定内失敗モード**であり、ソルバのバグを示すものではない。H/d≈10 とその温度応答を説得的に再現した実績が最も強いのは着火遅れ統計を運ぶクロージャ (transported PDF / CMC) である (これ以外の説得的な再現例は本調査では見つからなかった — 「のみ」と断定するには網羅性が足りない点は codex レビューの指摘どおり)。

## 1. 実験基準と「合うべき水準」

- 基準実験: Cabra et al., Proc. Combust. Inst. 29 (2002) / NASA CR-2002-212082。ジェット H₂ 25 %/N₂ 75 %, d=4.57 mm, 107 m/s, 305 K; vitiated coflow 1045 K (φ=0.25 希薄 H₂/air 燃焼ガス, X_O₂ 0.147, X_H₂O 0.099), ξ_st=0.47。**H/d≈10** (OH 発光基準, 数値側は Y_OH≈2×10⁻⁴ 等値面が標準)。
- **コフロー温度超敏感性 (最重要)**: Cao, Pope, Masri, Combust. Flame 142 (2005): 「T_c が 10 K (約 1 %) 下がると H/D は 2 倍になり得る」。熱電対の絶対精度は ~30 K (3 %)。同一条件のはずの Sandia と Sydney のリグで、同じ H/D≈10 を出すコフロー温度が 1045 K vs ~1022 K と **23 K ずれる**。Gordon et al., Combust. Theory Modelling 11 (2007): Rayleigh 実測 1058 K vs 公称 1045 K。
- **機構敏感性**: Li 機構は Mueller より早く着火し H が joint-PDF で 2–3 d、composition-PDF で 6 d 変わる (Cao 2005; Gordon 2007)。Benim & Pfeiffelmann, Energies 13 (2020): **GRI-3.0 の H₂ サブセットはこの火炎に使用不可** (EDC で全 T_c 吹き消え、PDF で H/d 55–75)。機構の影響は mixing model や C_φ より大きい (Cao 2005 の明示的結論)。
- **査読が受け入れる比較作法**: 絶対 H を公称 T_c で当てるのではなく、**T_c を実験不確かさ ±30 K の範囲で調整して リフトオフ高さを合わせた上で構造を比較**するのが標準 (Cao/Pope/Masri 2005 が明示: 1045→1031–1040 K を使用; CH₄ 系 LES-PaSR 2025 も 1350→1320 K)。より強いテストは **H vs T_c の応答曲線の傾き**。

つまり「バチバチに合う」は公称条件の単発比較では**実験自体が保証していない**。±(3–5) d / T_c チューニング込みが文献の実勢。

## 2. モデル別の再現実績

| クロージャ | Cabra H₂/N₂ での実績 | 出典 |
| --- | --- | --- |
| **平均速度 (no-TCI) / global+EDM/EDC** | **1010–1050 K の全域でリム付着** (global 機構)。着火遅れを表現できず早着火が構造的 | Benim 2020; De & Dongre arXiv:2101.09198 (steady flamelet で付着, DJHC) |
| **EDC** | H/d 8.5 (k-ε) / 5 (RSM) と点は近いが、**T_c 応答が浅く ~1020–1030 K 以下で突然吹き消え** = 応答曲線が誤り。DJHC でも既定定数は早着火 | Cabra 2002; Myhrvold CST 2006; Benim 2020; De FTC 2011 |
| **PaSR** | RANS-PaSR の公刊例なし。CH₄ Cabra の LES-PaSR (2025) は時間スケール定義の選択に強依存。**この火炎は chemistry-controlled で C_φ 1.5–2.5 の変化にほぼ非感応** (=C_mix ノブが効かない forge の観察と一致) | HAL 2025; Cao 2005 |
| **Transported PDF (RANS)** | 最有力。composition PDF+MC+k-ε (FLUENT) で **H/d 13.5 @1045 K** (+35 %) かつ応答曲線の傾き正しい; joint PDF は T_c 1031–1040 K で H/D=10 を再現。mixing model 間の差は <3 d (**EMST 必須ではない**) | Gordon 2007; Cao 2005 |
| **CMC (RANS 1st order)** | 実験・joint-PDF と良好一致 (radially-averaged, 実質 x×η 空間の追加 1D 解法)。条件付き速度/C_Φ の選択は低 T_c で効く | Patwardhan PCI 32 (2009); arXiv:2102.04854 |
| **LES-CMC** | FV LES-CMC (arXiv:2011.14083) が **T_c 無調整で H/d≈10** (基部 8.7–11.5 d 変動)。記録上最もクリーンな無調整一致 | 同左; Stanković/Merci は ±5 d |
| **Flamelet/FGM (steady)** | H/d≈5 で T_c にほぼ非感応 (着火物理なし)、DJHC では付着。**不可** | Benim 2020; De & Dongre |
| **UFPV (LES)** | CH₄ Cabra で ±10 %; ただし T_c 感度は実験の 1/4。揺らぎ項を切ると早着火に戻る | Ihme & See CF 157 (2010) |
| **quasi-laminar LES** | H/d≈10 を再現した公刊例なし (不在であり否定結果ではない)。解像混合だけでは早着火 (Ihme の C″=0 テスト) | — |

## 3. 付着のメカニズム (なぜ mean-rate では絶対に早着火するか)

- 安定化は**自着火** (Gordon 2007 の HO₂ 前駆蓄積・対流-反応バランス; Yoo, Sankaran, Chen JFM 640 (2009) DNS: 基部の平均流速は層流火炎速度の ~10 倍で伝播説を棄却)。
- 着火核は**最反応性混合分率 ξ_MR≈0.05 (超希薄) かつ低スカラー散逸のポケット**でのみ発生 (Mastorakos PECS 35 (2009))。coflow が H₂ クロスオーバー温度超のため、可着火な希薄混合気はリップ直近まで存在する。セル平均 (Y,T) で Arrhenius を評価すると PDF の裾と χ 履歴が消え、そこが即着火する → 火炎基部が上流へ這い、付着。
- 同じ理由で**混合速度リミッタ (PaSR/EDC) の定数調整では直らない** (欠けているのは混合律速でなく着火遅れ統計)。forge の C_mix 非感応はこの文献知見と正確に一致。

## 4. forge への含意

1. **リップ付着はソルバ無罪の証拠パターン**。laminar-chemistry + PaSR 非感応 + 付着、の 3 点セットは文献の失敗モードと完全一致。
2. **「バチバチに合わせる」なら**: 最安の実績ルートは (a) **RANS 1st-order CMC** (x×η の追加解法; Patwardhan 2009) または (b) **Lagrangian composition PDF** (MC/IEM 混合で十分, ISAT 併用; Gordon 2007)。EDC/FGM/MEPDF は実績が悪く非推奨。いずれも数セッション規模の実装。
3. **合格判定の作法**: 公称 1045 K の単発一致ではなく、T_c を ±30 K 内で振った **H vs T_c 応答曲線** (Wu 2003 / Gordon データ, Cao 2005 Fig.1 / Benim Figs.8–13) で評価する。
4. **機構**: 現行 Jachimowski (1988) の 1000–1100 K 帯の着火遅れは文献比較群 (Mueller / Li 2004 / Keromnes 2013) で未評価。TCI 実装の前に **Cantera で ξ_MR 近傍条件の着火遅れを Li/Mueller と比較**する (安価; H が数 d 動く要因)。GRI-3.0 サブセットは使用不可。
5. **下流検証で閉じる場合**: 付着でも z/d 26 半径 T は実験整合済み (case/48)。加熱器用途 (出口状態予測) には「付着火炎は既知限界」と明記して下流量で閉じるのは文献的に防御可能 — ただし「リフトオフ高さを予測した」とは主張しない。

## 出典 (主要・全文確認済み)

- [Cabra 2002 NASA CR](https://ntrs.nasa.gov/citations/20030014612) / [Cao, Pope, Masri CF 2005](https://tcg.mae.cornell.edu/pubs/Cao_PM_CF_2005.pdf) / [Gordon et al. CTM 2007](https://tcg.mae.cornell.edu/pubs/Gordon_MPG_CTM_07.pdf)
- [Benim & Pfeiffelmann Energies 2020](https://www.mdpi.com/1996-1073/13/1/152) / [RANS-CMC 感度 arXiv:2102.04854](https://arxiv.org/abs/2102.04854) / [FV LES-CMC arXiv:2011.14083](https://arxiv.org/pdf/2011.14083)
- [De et al. DJHC-EDC (arXiv:2101.08764)](https://arxiv.org/abs/2101.08764) / [De & Dongre MILD-TCI (arXiv:2101.09198)](https://arxiv.org/pdf/2101.09198) / [Ihme & See CF 2010](https://web.stanford.edu/group/ihmegroup/cgi-bin/MatthiasIhme/wp-content/papercite-data/pdf/ihme_see_cf2010.pdf)
- [LES-PaSR Cabra CH₄ 2025 (HAL)](https://hal.science/hal-05386199v1/document) / [Mastorakos PECS 2009](https://www.sciencedirect.com/science/article/abs/pii/S0360128508000415) / [Yoo, Sankaran, Chen JFM 2009](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/threedimensional-direct-numerical-simulation-of-a-turbulent-lifted-hydrogen-jet-flame-in-heated-coflow-flame-stabilization-and-structure/F90F81BFEB5324713BC67D3B8A762579)
- 抄録のみ (paywall): Myhrvold CST 2006 (EDC 詳細), Patwardhan PCI 2009 (CMC の H/d 数値), Jones & Navarro-Martinez CF 2007 (stochastic fields)

関連: [chemistry-finite-rate-h2-survey.md](chemistry-finite-rate-h2-survey.md), [../sessions/chemistry-cabra-handover-2026-09-05.md](../sessions/chemistry-cabra-handover-2026-09-05.md)

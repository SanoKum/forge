# DES / WMLES 最新動向サーベイ（超音速ノズル・ピントルバルブ適用向け）

## メタ

- **area**: `turbulence`
- **status**: `draft`（調査専用・コード変更なし）
- **related_docs**:
  - `docs/turbulence/theory.md`
  - `docs/turbulence/implementation.md`
- **related_plans**:
  - [`turbulence-iddes-sst.md`](../../design/active/turbulence-iddes-sst.md) — 本サーベイを受けた IDDES 実装計画
- **created**: `2026-06-21`
- **owner**: `CFD Dev`
- **調査実施**: deep-research ワークフロー（94 エージェント, 21 ソース取得, 22 クレーム検証: 11 確認 / 11 棄却）

---

## 背景・目的

超音速ノズル内部流れ（ピントルバルブ・チョーク流れ）に向けた非定常解析手法の選定。
URANS から DDES/IDDES/WMLES のどれを採用すべきかを 2020–2025 年文献で検証した。

---

## 確認された主要知見（3-vote 検証通過）

### 1. DDES の設計意図（信頼度：高、3-0 票）

**出典**: Spalart et al. (2006) *Theoretical and Computational Fluid Dynamics* 20, 181–195  
**補強**: Chaouat (2017) レビュー `https://link.springer.com/article/10.1007/s10494-017-9828-8`

> "DES is very sensitive to the grid-size. In particular, the gray area where the model varies
> from URANS to LES mode may be problematic unless the separation is abrupt and fixed by the
> geometry."

DES97 はグリッドサイズに対して非常に敏感で、グレーゾーン問題が起きる。
DDES はこれを「**resistant to ambiguous grid densities**」と銘打って、
f_d シールド関数（剥離を検知して RANS→LES スイッチを遅延）で解決した。
2024 年の文献（arXiv:2511.08257, ScienceDirect DDES-AGAC 2024）でも Kelvin–Helmholtz 遅延として現役の問題として確認。

### 2. SBLI には IDDES-SST + WENO が現在の主流（信頼度：高、2-1 票）

**出典**: Troshin & Bakhne (2024) *Mathematical Models and Computer Simulations* 16, 100–111.  
`https://link.springer.com/article/10.1134/S2070048224010113`

衝撃波-境界層干渉（SBLI）の数値研究で **IDDES-SST + WENO（上流・対称両スキーム）**
を2グリッド密度で検証。評価指標：剥離長さ・壁面圧力分布・摩擦係数・1点圧力 PDF（2箇所）。
2024年時点で SBLI に対する IDDES-SST の実績として最も直近かつ包括的な一次文献。

### 3. DES のログ層ミスマッチ（信頼度：中、2-1 票）

**出典**: Gopalan, Heinz & Stöllinger, *JCP* 2013. `https://www.uwyo.edu/heinz/jcp13.pdf`

付着 BL の RANS-LES 界面で DES は ~15.5% の摩擦係数誤差（LUM は ~9.8%）。
方向性は最新文献で一致。ただし Re_tau 値の独立検証は画像 PDF のため困難。

### 4. WMLES 平衡壁モデルは高 Re で失敗（信頼度：高、3-0 票）

**出典**: Iyer & Malik, ICCFD11 2022. `https://ntrs.nasa.gov/api/citations/20220009649/`

Gaussian bump（Mach 0.2）：
- Re_L = 1×10⁶ → WMLES が剥離を正確に捉える
- Re_L = 2×10⁶ → **実験では剥離あるのに WMLES が剥離を予測できない**

Vreman 係数 c=0.07（標準）では NG、c=0.025 で回復（Agrawal CTR 2022）。
**平衡壁モデル + 標準 SGS 係数の組み合わせは高 Re WMLES では危険。**

### 5. ML 壁モデルは亜音速限定（信頼度：高、3-0 票）

**出典**: Davidson, arXiv:2410.17767 (Oct 2024, updated Mar 2025)

IDDES データで訓練した ML 壁関数は 5 つの亜音速ケースで検証済み。
**超音速・衝撃波を含む流れへの適用は一切未検証。** 超音速ノズルへの転用は不可。

### 6. 「万能な最良手法」は存在しない（信頼度：高、2-1 票）

**出典**: Chaouat (2017) レビュー + 2024–2025 年複数論文

> "the most appropriate method for a particular application will depend on the expectations
> of the engineer and the computational resources the user is prepared to expend on the problem."

2024 SBLI 論文・AIAA 2025-0891（ZDES/HRLMEC が IDDES を凌駕するケースあり）・
2025 超音速 shock-BL 論文いずれも「ケースによる」結論。

---

## 棄却された俗説（検証 0-3 票）

| 俗説 | 棄却理由 |
|------|---------|
| IDDES は剥離長を 15–28% 過小評価する（常に） | 根拠なし・過度な一般化 |
| LES コストは Re^1.76 でスケールする | 独立検証不可 |
| RANS コストは Re によらず一定 | 1-2 票のみ |
| RANS は本質的に非定常を再現できない | 0-3 票で完全棄却（URANS は再現可） |
| 平衡壁モデルは逆圧力勾配で常に不十分 | 1-2 票のみ、フロー依存 |

---

## 未解決問題（ピントルバルブ直結）

調査した 21 ソースの中に**超音速ノズル内部流れ・ピントルバルブを直接扱った
DES/WMLES 研究は皆無**だった。以下がオープンクエスチョン：

1. チョーク流れ + 斜め衝撃波列 + ピントル面再循環の組み合わせで DDES/IDDES はどう振る舞うか
2. 超音速領域での WMLES の格子収束基準（y+ とスパン方向セル数の最小値）は？
3. ML 壁モデルを圧縮性・衝撃波含む流れに一般化できるか？
4. 強い 3D 非定常衝撃パターン + 再循環でも IDDES-SST の WMLES モードは安定か？

---

## 手法選定マトリクス

```
ケース                              推奨手法
──────────────────────────────────────────────────────────
超音速ノズル（主に付着 BL）         SST-DDES + automatic wall treatment (y+~30–80)
ピントルバルブ（幾何固定剥離）       SST-DDES（剥離が幾何固定なら有利）
                                    → 再循環支配なら IDDES
衝撃波-BL干渉をきちんと捉えたい      IDDES-SST + WENO（Troshin 2024 実績）
避けるべき組み合わせ                 WMLES + 平衡壁モデル + 高 Re（失敗確認）
                                    ML 壁モデル + 超音速（未検証）
```

---

## Forge への実装コスト観点

- **SST-DDES が最小コスト**: 既存 SST に f_d シールド関数 + length scale 修正のみ。
- **SST-IDDES**: さらに f_B・f_e の追加が必要。WMLES モードは解像乱流流入がない限り不活性。
- 実装計画は [`turbulence-iddes-sst.md`](../../design/active/turbulence-iddes-sst.md) を参照。

---

## 引用ソース（一次・二次）

| 品質 | 出典 |
|------|------|
| 一次 | Troshin & Bakhne (2024), Springer: IDDES-SST for SBLI |
| 一次 | Iyer & Malik (2022), NASA NTRS: WMLES Gaussian bump failure at high Re |
| 一次 | Davidson (2024), arXiv:2410.17767: ML wall function (subsonic only) |
| 一次 | Gopalan et al. (2013), JCP: DES log-layer mismatch quantification |
| 二次 | Chaouat (2017), Flow Turbulence Combust: hybrid RANS-LES review |
| 原典 | Spalart et al. (2006), TCFD 20: DDES 論文 |
| 原典 | Shur et al. (2008): IDDES 論文 |

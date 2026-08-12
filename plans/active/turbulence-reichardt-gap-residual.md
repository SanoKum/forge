# 壁解像プロファイルと Reichardt 相関の 5% 差 — 残る候補の切り分け

## メタ

- **area**: `turbulence / boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5(e)
- **related_plans**:
  - [`turbulence-node-wf-representative-point.md`](../accepted/turbulence-node-wf-representative-point.md) (前段: 候補 A/B を棄却し SU2 で forge 固有性を除外)
- **created**: `2026-08-13`
- **owner**: `sano`

## 1. 事実 (確定済み)

case/26 平板の**壁解像**解は、$y^+\approx30$ で Reichardt 対数則より約 **5% 低い**:

| ソース | 収束 | 実測/Reichardt |
| --- | --- | --- |
| forge 壁解像 (`run_0064` 等) | 全残差 `NOT CONVERGED` / $\theta$,$C_f$ は `STEADY` | 0.9492–0.9504 |
| SU2 壁解像 (`run_0049`) | 全残差収束 (rms −5.0〜−12.6) | 0.9477–0.9490 |

**両者 0.1–0.2% で一致**。独立した 2 コード (一方は全残差収束) が同じ差を示す。

### 既に除外した候補

| 候補 | 結果 |
| --- | --- |
| forge 固有の離散誤差・物性定義・壁応力定義 | **除外** (SU2 が同値) |
| 自由流乱流条件 ($\mu_t/\mu$ 4 桁) | **除外** (準定常量に検出可能な変化なし) |
| 補間手法 (IDW vs 要素補間) | **除外** (差 0.115%) |

## 2. 残る候補

1. **有限 $Re$** — この $Re_\theta$ 帯 ($\approx$1500–2500) では対数則の漸近形が
   まだ成立しきっていない可能性。Reichardt 係数 ($\kappa$=0.41, $B$ 相当 7.8) は
   高 $Re$ フィット。
2. **入口/前縁の扱い** — 前縁からの発達距離、入口 BL 厚さ、
   前縁特異点近傍の履歴が $y^+\approx30$ のプロファイルに残っている可能性。
3. **Reichardt 相関そのものの適合限界** — Reichardt 式は buffer 層をつなぐ内挿式であり、
   $y^+$=30 は遷移域の端。SST の解と系統差があっても不思議ではない。

## 3. 切り分けの方針 (案)

- **1 の検証**: 同一セットアップで $Re$ を上げ (または下流station で $Re_\theta$ を伸ばし)、
  実測/Reichardt が 1 へ近づくかを見る。近づけば有限 $Re$ 効果。
  上表の SU2 は $x$ 増加 ($Re_\theta$ 増加) で 0.9490→0.9477 と**わずかに離れる**ので、
  この方向は一見不利 — まず素直に確認する。
- **2 の検証**: 入口位置・前縁扱いを変えた run で同量を比較。
- **3 の検証**: 他の壁法則 (Spalding / Nichols–Nelson / 純対数則) との比を並べ、
  差が Reichardt 固有か全法則共通かを見る。**共通なら相関側でなく解側**。

## 4. 判定

この差は**壁関数の欠損とは別物**であることが確定している (壁解像でも出る)。
したがって本 plan は**参照基準の精度問題**であり、node 壁関数の改良とは切り離して扱う。

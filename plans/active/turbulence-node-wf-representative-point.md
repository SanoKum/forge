# node 壁関数の代表点 (Normal_Neighbor) 診断 — 3 要因の分離

## メタ

- **area**: `turbulence / boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5(e)
- **related_plans**:
  - [`turbulence-node-wf-omega-source.md`](turbulence-node-wf-omega-source.md) (親: $\omega$ 側の切り分けはここで完了。E3 は不採用)
- **created**: `2026-08-13`
- **owner**: `sano`

## 1. なぜここを見るか

親 plan で $\omega$ 側の候補は出し切った:

| 介入 | case/26 $C_f/C_{f,\mathrm{KS}}$ | case/40 $\tau_w$/y+1 基準 |
| --- | --- | --- |
| 現行 | 0.7615 | 0.945 |
| $\omega$ pin (Dirichlet) | 0.8594 | 1.237 (棄却) |
| **E3 ($P_\omega$ 源項整合)** | **0.8445** | **≈1.061 (不合格)** |
| 目標 | 0.94±0.02 | ≈1.0 |

**両ケースが約 11–12% の同方向応答**を示すのに、**基準に対する符号が逆** (case/26 は不足、
case/40 は過剰) である。→ **単一ノブでは閉じない**。$u_\tau$ を決める**入口側** =
代表点 (Normal_Neighbor) の選び方と、そこで読む速度を疑う。

## 2. ★ 前提: 「第一内層速度が低い = 代表点選択の誤り」と即断しない

低い速度は**低い壁応力と自己整合した解の結果**でもある
(親 plan で確認: 課した $\tau_w$ と BL 発達が互いに整合)。
**壁解像場を同じ物理位置へ補間して初めて**、サンプリングの問題か場の問題かを区別できる。
したがって §3 の 3 要因を**分離してから**アルゴリズムに手を入れる。
**いきなり選択アルゴリズムを変更しない。**

## 3. 3 要因の分離 (この順)

### 3.1 代表点選択の幾何

各壁点について出力する:

- 選ばれた `irep`
- $y=(\mathbf{x}_I-\mathbf{x}_W)\cdot(-\hat{\mathbf n})$ (壁関数が使う法線距離)
- 実距離 $|\mathbf{x}_I-\mathbf{x}_W|$
- `best_cos` (内向き法線との cos。現行の選択基準)
- 接線方向オフセット $|(\mathbf{x}_I-\mathbf{x}_W) - y\,(-\hat{\mathbf n})|$
- `wall_dist[irep]`$/y$ (壁距離場との整合)

**現行は「内向き法線との cos 最大」だけで選んでいる**ので、曲面やコーナーでは
「法線上で最適なサンプル点」とは限らない。case/26 (平板) と case/40 (曲面ベル) で比較する。

### 3.2 代表点の速度場

- $U_t(\texttt{irep})$ (壁接線速度の大きさ)
- $\rho,\nu$ (代表点の値)
- 解いた $u_\tau$
- 壁法則残差 ($U_t/u_\tau - u^+(y^+)$ の収束後の値)
- 感度 $du_\tau/dU_t$

case/26 と case/40 で比較し、**どちらのケースで代表点速度が基準からどうずれているか**を見る。

### 3.3 壁解像基準を**同じ物理位置**でサンプルする (核心)

壁解像 run を代表点と同じ物理距離 $y_{\rm rep}$ へ法線補間し

$$U_{t,\rm coarse}(y_{\rm rep})\quad\text{vs}\quad U_{t,\rm ref}(y_{\rm rep})$$

を比較する。これで次の 3 つが分離できる:

| 一致/不一致 | 意味 |
| --- | --- |
| $U_{t,\rm coarse}\ne U_{t,\rm ref}$ | **粗メッシュの速度場自体が違う** (サンプリングの問題ではない) |
| $U_{t,\rm coarse}=U_{t,\rm ref}$ だが $u_\tau$ が違う | **Reichardt 則に基準速度を入れても真の $u_\tau$ を再現しない** (壁法則側の問題) |
| 代表点の $y$ が法線上の期待位置とずれる | **隣接ノードの選び方が悪い** (§3.1 と連動) |

基準は case/26 = `run_0050` (node 壁解像) と `run_0049` (SU2 壁解像)、
case/40 = y+1 low-Re 基準 run。

## 4. 検証と判定

- 診断は既存と同じ env ゲート方式 (OFF で確保も出力もしない、既定経路ビット不変)。
- 判定量は親 plan §7.4 の 3 点セットを引き継ぐ ($C_f/C_{f,\mathrm{KS}}$ は
  [`case/26.flat_plate_sst/tools/cf_retheta_analysis.py`](../../case/26.flat_plate_sst/tools/cf_retheta_analysis.py) で再生成)。
- **case/40 の y+1 基準との比較を、case/26 と同一定義で行う正式ツール化**も本 plan の範囲
  (現状は $\tau_w$ 比からの推定で、基準 run との直接再測定になっていない)。

## 5. 完了条件

- [ ] §3.1 幾何診断の出力と case/26 / case/40 比較
- [ ] §3.2 代表点速度・$u_\tau$・壁法則残差・感度の比較
- [ ] §3.3 壁解像基準の同一位置サンプリングによる 3 分離
- [ ] case/40 の y+1 基準比較を同一定義で正式ツール化
- [ ] 結果に基づく修正方針の決定 (アルゴリズム変更はこの後)
- [ ] `methods/turbulence/` 更新 → `status: done` で `accepted/` へ移動、[`plans/README.md`](../README.md) 同期

## 6. 変更ログ

- `2026-08-13` — 起票。親 plan で $\omega$ 側 (pin / limiter / $P_\omega$ 源項) を出し切り、
  **両ケースが同方向に約 11–12% 応答するのに基準に対する符号が逆**という結果を受けて、
  $u_\tau$ の入口側 (代表点) へ調査を移す。**「第一内層速度が低い = 代表点選択の誤り」と
  即断しない**ことを §2 に明記。

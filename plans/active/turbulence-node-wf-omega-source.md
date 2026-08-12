# node 壁関数 DOF の $\omega$ 生産を壁法則整合形にする (3 分離実験 → 恒久修正)

## メタ

- **area**: `turbulence / boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5 (automatic wall treatment)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §3.7
- **related_plans**:
  - [`../accepted/turbulence-node-wall-function-coverage.md`](../accepted/turbulence-node-wall-function-coverage.md) (第一内層への `wf_pk` 適用を入れた計画 = 本件の前段)
  - [`../accepted/turbulence-sst-su2-taw-coupling.md`](../accepted/turbulence-sst-su2-taw-coupling.md) §10c (`nodeOmegaWfDirichlet` を case/40 で棄却した記録)
- **created**: `2026-08-12`
- **owner**: `sano`

## 1. 目的 (現象)

node × SST 壁関数の $\tau_w$ が cell より系統的に低い (case/26 平板で $C_f$/壁解像真値 =
**node 0.843 / cell 1.027**)。切り分けの結果、**第一内層ノードの $\omega$ が壁ノードのピン値より
高く跳ねる**ことが分かった (壁 1.146e5 → 第一内層 **2.645e5** → 第二 9.635e4)。
高い $\omega$ が $\mu_t=\rho k/\omega$ を下げ、同時に散逸 $\beta^*\rho k\omega$ を上げて $k$ も下げる。

本 plan は、**主因を分離して確定し**、恒久修正を入れることを目的とする。

## 2. 確定している事実 (case/26, x=0.6, 詳細は [case README](../../case/26.flat_plate_sst/README.md))

- **壁関数 $P_k$ の積分量は $\tau_w$ をほぼ支配しない**。`nodeKwfDirichlet: 0` で実効生産が
  0.311→0.977 (**3.1 倍**) になっても $u_\tau$ は +0.55% (`run_0045`)。
  → 「生産規約 (cell factor 2.000 vs node 1.564) が真因」は**反証済み**。
- **$\omega$ スパイクは局所源項平衡そのもの**: 解像 $S$=1.0750e5, $F_1$=0.621 →
  $\sqrt{\alpha/\beta}S$=2.754e5 に対し実測 $\omega$=2.645e5 (比 0.960)。
- **$P_\omega$ は `wf_pk` と結合していない**: [`ransSource_d.cu:178`](../../solver_density_cuda/cuda_forge/ransSource_d.cu)
  は $P_\omega=\alpha\rho S^2$。
- **未解像な離散ひずみが過大**: 壁法則整合 $S_{\rm wf}=u_\tau^2 g/\nu$ に対し
  解像 $S$ は y+27 で **2.09 倍** ($S^2$ で **4.36 倍**)。
- **`nodeOmegaWfDirichlet` の A/B は交絡している**: このスイッチは $\omega$ ピンと
  **SST shear limiter $\max(a_1\omega,SF_2)$ の迂回**を同時に行う
  ([`turbulent_viscosity_d.cu:222`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu))。
  $C_f$ 0.843→0.957 は両者の合成効果。

**未分離**: 効いているのが (i) $\omega$ 状態か (ii) strain limiter 迂回か (iii) 両方か。

## 3. 3 分離実験 (先にこれをやる)

いずれも case/26 `run_0042` と同一 IC/設定、20000 step、x=0.6 の $u_\tau$/$C_f$/$\omega$/$\mu_t$ を比較。
**$P_\omega$/$D_\omega$/交差拡散/$\omega$ 輸送の各収支を診断出力に足してから**行う。

| # | 内容 | 実装 | 判定 |
| --- | --- | --- | --- |
| E1 | **$\omega$ ピンのみ**、limiter は従来式のまま | `turbulent_viscosity_d.cu` の `wfPin` 分岐を外し、ピンは残す | $C_f$ が 0.957 のどれだけを再現するか |
| E2 | **limiter bypass のみ**、$\omega$ は free | `wfPin` 相当の条件を `wf_pk>=0` に付け替え、$\omega$ はピンしない | 同上 |
| E3 | **ピンも bypass もせず、第一内層の $P_\omega$ だけ壁法則整合形へ置換** | 下記 §4 | **本命** |

E1+E2 の和が E3 の効果とどう並ぶかで機構が確定する。

## 4. 恒久修正候補 (E3)

壁関数 DOF (`wf_pk>=0`) では、$P_\omega$ の入力ひずみを解像 $S$ でなく**壁法則整合**にする:

$$S_{\rm wf}=\frac{u_\tau^2}{\nu}g,\qquad
\nu_{t,\rm wf}=\nu\left(\frac1g-1\right),\qquad
P_{\omega,\rm wf}=\alpha\rho S_{\rm wf}^2=\alpha\frac{P_{k,\rm wf}}{\nu_{t,\rm wf}}$$

($g=du^+/dy^+$ は既存の Reichardt 微分 `wallLaw_reichardt_duplus_dyp`、$P_{k,\rm wf}$ は既存 `wf_pk`)。
これは **$P_\omega$ を limited $P_k$ と結合させる SST-2003 の形**でもあり、$k$ 側の壁関数置換と
整合する (現行は $k$ だけ壁関数値・$\omega$ は解像 $S$ という**不整合**な組合せになっている)。

- **要確認 (Codex 指摘, 未検証)**: SST-2003 では $P_\omega$ は limited $P_k$ と結合すべきで
  $\alpha\rho S^2$ は erratum、という指摘があった。NASA Turbulence Modeling Resource の原文で
  裏を取ってから本文に反映する (**現時点では未確認事項として扱う**)。もし一般に正しいなら
  cell/内部場も含む別 plan に切り出す (本 plan は壁関数 DOF 限定に留める)。
- cell モードは既定でビット不変にすること (cell は第一セルで $\omega$ ピン済みのため本修正の
  対象 DOF が存在しない、という形にできるか要確認)。

## 5. 検証

- **case/26 平板** (本 plan の主戦場): $C_f$/壁解像真値が cell (1.027) に近づくこと。
  y+ 掃引 (`run_0042`/`0040`/`0043` = y+27/52/102) で **y+ 非依存が保たれる**こと。
- **case/40 ベルノズル**: $\tau_w$/y+1 真値が現行 0.945 から**過大化しない**こと
  (`nodeOmegaWfDirichlet` は 1.237 に過大化して棄却された。E3 が同じ轍を踏まないことが採否の分かれ目)。
- **OFF 回帰**: 既定でビット不変 (cell 全ケース + node 非壁関数ケース)。
- **判定**: 残差はこの case で構造プラトーなので、**派生量 ($u_\tau$) の準定常**を
  スナップショット間ドリフトで確認する (現行 6 run は 0.001–0.075%)。

## 6. 完了条件

- [ ] $P_\omega$/$D_\omega$/交差拡散/輸送の壁近傍収支を診断出力に追加
- [ ] E1/E2/E3 の 3 分離実験で機構を確定
- [ ] E3 が case/26 を改善し case/40 を過大化させないことを確認
- [ ] `methods/turbulence/` の §6.5 に仕様追記、`methods/index.md` 同期
- [ ] `status: done` にして `plans/accepted/` へ移動、[`plans/README.md`](../README.md) 同期

## 7. 変更ログ

- `2026-08-12` — 起票。case/26 の切り分け (`run_0040`–`run_0046`) と 2 度のレビュー訂正
  (①`nodeKwfDirichlet` 既定 1 による $\Sigma P_kV$ の誤読 → 生産規約説を撤回、
  ②$P_\omega=\alpha\rho S^2$ は `wf_pk` 非結合・`nodeOmegaWfDirichlet` は limiter bypass と交絡)
  を反映した状態で作成。**主因は未分離**であることを明記。

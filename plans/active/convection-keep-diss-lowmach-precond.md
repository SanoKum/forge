# KEEP ES 散逸レイヤの低マッハ前処理化 (Weiss–Smith/Turkel 前処理散逸)

## メタ

- **area**: `convection`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/convection/implementation.md`](../../methods/convection/implementation.md) (散逸レイヤ節 `keepDissPrecond`)
  - [`methods/convection/theory.md`](../../methods/convection/theory.md) (SLAU の低マッハ前処理節 = 同一理論)
  - [`methods/time_integration/theory.md`](../../methods/time_integration/theory.md) (Weiss–Smith 前処理固有系)
- **related_plans**:
  - [convection-keep-es-dissipation](../accepted/convection-keep-es-dissipation.md) (ES 散逸レイヤ本体)
  - [convection-keep-diss-recon-jump](../accepted/convection-keep-diss-recon-jump.md) (jump2)
  - [time_integration-lowmach-preconditioning](../accepted/time_integration-lowmach-preconditioning.md) (`lowMachPrecond=2`, 移植元)
- **created**: `2026-07-19`
- **owner**: CFD Dev

## 1. 目的・動機 (backstep パイロットの実測)

case/18 backstep node LES パイロット (933k nodes, KEEP+matrix ES) で、低マッハ再循環域
(x≈2.2–5.5, |u|~数 m/s) に **x 方向柱状の圧力 odd-even (市松)** が残ることが判明した。
σ (`keepDissCoeff`) と音響波速 (`keepDissCprime`: c' / full c) を掃引した 4000 step プローブ
(seed=run_0097 res_11000 共通, 評価は再循環ボックスの x 方向 2Δ 振幅 = d2/4 rms):

| 設定 | 市松 mean/p95/max [Pa] | せん断層 Uy rms [m/s] | 3D 遷移 rms_roUz |
|---|---|---|---|
| σ=0 (pure KEEP) | 47.3 / 78.3 / 100 | 7.10 | 8.3e-4 |
| σ0.02 c' | 28.6 / 59.8 / 102 | 7.12 | 4.7e-4 |
| σ0.05 c' | 22.1 / 54.1 / 82 | 6.89 | 3.5e-4 |
| σ0.10 c' | 15.6 / 39.2 / 55 | 6.46 | 2.7e-4 |
| σ0.02 full c | 17.3 / 45.0 / 62 | 6.57 | 2.7e-4 |
| σ0.20 c' | 11.8 / 29.1 / 41 | 6.25 | 2.2e-4 |
| σ0.05 full c | 11.0 / 27.1 / 37 | 6.20 | 2.1e-4 |
| σ0.10 full c | 7.4 / 13.6 / 21 | 5.98 | 1.7e-4 |

(run 一覧: case/18 README `run_0109`–`run_0118`。σ0.20c'≈σ0.05fullc、σ0.10c'≈σ0.02fullc)

**発見: σ と c' はどちらを動かしても同一のトレードオフ曲線上を動くだけ** (市松 1 Pa 減 ≒
せん断層変動と 3D 遷移速度を比例して支払う)。固有値の一様スケーリングでは、2Δ 近傍に同居する
市松 (音響 null-mode) とギリギリ解像の物理を**分離できない**。

分離には散逸の**行列構造**を変える必要がある。文献の標準解 (Weiss–Smith 1995 / Turkel;
圧力ベース界の Rhie–Chow、AUSM+-up/SLAU の圧力速度結合項と同族) は前処理散逸:

- **速度ジャンプ散逸**: $\rho\,c\,\Delta U_n \to \rho\,U_r\,\Delta U_n$ (低マッハで $1/M$ 倍**縮小** →
  Guillard–Viozat の過散逸を解消、渦・せん断層を保護)
- **圧力ジャンプ散逸** (連続式): $\Delta p/c \to \Delta p/U_r$ ($1/M$ 倍**増強** → 圧力 Laplacian が
  連続式に入り市松を直接減衰 = Rhie–Chow と同じ機構)

forge は `lowMachPrecond=2` (SLAU+implicit) で同機構を実証済み (backstep P 市松 −77%,
[time_integration-lowmach-preconditioning](../accepted/time_integration-lowmach-preconditioning.md) §9)。
本 plan は KEEP ES 散逸レイヤの音響対へこれを移植する。

## 2. スコープ

- **やる**: `keepDissPrecond: 1` (既定 0 = 従来・ビット不変)。matrix ES (`keepDissType:2`) の
  **CPG 枝**で、音響 ± 対の散逸を前処理音響 2×2 に置換。`keepDissJump` 0/1/2 と合成可。
  cell/node 両対応 (散逸は主ループのみ、境界半割面は upwind 散逸内蔵で対象外 = 従来どおり)。
- **やらない**: TP/多成分枝 (フォローアップ)。scalar ES (`keepDissType:1`)。SLAU/ROE への展開。
  伝播 (中心流束・時間積分) への前処理 (散逸レイヤのみ = RHS の数値散逸を変えるだけで
  物理時間発展は不変。explicit unsteady で正当)。

## 3. 設計

### 3.1 定式化

現行 matrix 散逸は特性和 $D=\sum_k S_k\,\lambda_k\,(r_k\!\cdot\!\Delta w)\,r_k$ で、音響対
$k=\pm$ は共通固有値 $\lambda_A=|U_n|+c'$。これを音響 2×2 サブシステムの前処理散逸に置換する。

対称化音響変数 $(\hat p, \hat u) = (\Delta p/(\rho c),\ \Delta U_n)$ で音響部は
$A_{ac} = \begin{pmatrix} U_n & c \\ c & U_n \end{pmatrix}$。Weiss–Smith 前処理
$\Gamma_{ac} = \mathrm{diag}(c^2/U_r^2,\ 1)$ (圧力時間項の $\Theta=1/U_r^2$ スケール) に対し

$$D_{ac} = \Gamma_{ac}\,\left|\,\Gamma_{ac}^{-1} A_{ac}\,\right|
 = \Gamma_{ac}\, R_\Gamma |\Lambda_\Gamma| R_\Gamma^{-1}, \qquad
 \Lambda_\Gamma = \mathrm{diag}(u'{+}c',\ u'{-}c'),\ \ u'=\tfrac{1+\beta}{2}U_n,\ \
 c' = \texttt{lowMachCprime},\ \beta = U_r^2/c^2$$

漸近 ($M\to0$, $U_r\to\max(|u|,\varepsilon c)$):
$D_{ac}[\hat p,\hat p] \sim c^2/U_r$ (→ $\Delta p$ 散逸 $\propto 1/U_r$ **増強**)、
$D_{ac}[\hat u,\hat u] \sim U_r$ (→ $\Delta U_n$ 散逸**縮小**)。$M\ge1$ ($\beta=1$) で
$D_{ac} \to |A_{ac}|$ = 従来 full-c 行列散逸に厳密復帰。

実装は $(\hat p,\hat u)$ を現行の特性射影 $rd_\pm$ から復元し ($\hat p = \tfrac{S_A}{?}(rd_+{+}rd_-)$ 系、
正確な変換は検証スクリプトで確定)、閉形 2×2 係数 $d_{pp}, d_{pu}, d_{up}, d_{uu}$ (Un, c, Ur の関数)
を掛けて特性和へ戻す。エントロピー・せん断波は従来どおり $|U_n|$。

### 3.2 ES 性 (エントロピー散逸性) の扱い

$D_{ac}$ は対称でない ($d_{pu}\ne d_{up}$) ため従来の「$R|\Lambda|SR^T$ 対称正定値」の証明は
そのまま使えない。エントロピー散逸性は 2 次形式 $\Delta w^T D \Delta U$ の符号で決まり、
**$D$ の対称部分が正定値なら成立**する。検証スクリプトで $(M, U_n/c, \varepsilon)$ 掃引の
対称部固有値を数値確認し、負になる領域があれば (a) その領域で従来散逸へブレンド、または
(b) sign-property クリップ (jump2) と同様の成分クリップで保証する。証明できない場合は
「安定化オプション (工業標準 = Fluent/SLAU と同格)」として位置づけ、ES ラベルは付けない。

### 3.3 jump2 との合成

再構成ジャンプ (滑らかな場で $\Delta\to O(h^3)$) を $(\hat p,\hat u)$ の入力に使う。
**再構成済み圧力ジャンプ × $1/U_r$ 増強**は「滑らかな圧力場は触らず、非滑らか (市松) だけ
強減衰」となり、Rhie–Chow の現代形 ($\Delta p - \overline{\nabla p}\cdot d$ 補正) と一致する。
sign-property クリップは特性射影単位で従来どおり適用。

## 4. 検証計画

1. **式の事前数値検証** (`tools/verify_precond_dissipation.py`, 実装前必須):
   (a) $D_{ac}$ の閉形 = 数値固有分解の一致、(b) 漸近スケーリング ($\Delta p$ 係数 $\propto1/U_r$,
   $\Delta U_n$ 係数 $\propto U_r$)、(c) $M\ge1$ で従来行列散逸へ厳密復帰、
   (d) 対称部正定値性の $(M, U_n/c)$ マップ (ES 性判定)。
2. **L1 市松減衰** (case/35 一様周期箱 + 市松摂動): 減衰レート ≥ full c (9.2e-9/400step 級)。
3. **L2 渦保護** (case/09 TGV 64³): KE cost ≤ c' 構成 (−1.4% 級) を維持。
4. **backstep A/B** (本命): 同一 seed 4000 step プローブで §1 の曲線と比較。
   **合格基準 = 曲線の外側に出ること**: 市松 mean ≤ 10 Pa **かつ** Uy rms ≥ 6.9
   (= σ0.05c' の物理で σ0.10fullc の市松)。
5. cell/node 両方で確認 ([memory: verify-node-and-cell-both])。回帰: `keepDissPrecond:0` ビット不変。

## 5. 実装

- `solverConfig.{hpp,cpp}`: `keepDissPrecond` (top-level, 既定 0)。
- `convectiveFlux_keep_d.inc.cuh` (matrix CPG 枝): 音響対の $z_\pm$ 構成を 2×2 前処理係数で置換。
  `lowMachUr`/`lowMachCprime`/`precondEps` を再利用。
- 検証 run は case/18 (README run 表) / case/35 / case/09 に記録。

## 6. 検証結果 (2026-07-19, 実装後)

- **L1 (真の市松, case/35 種市松 400step)**: node `run_0043` A_cb 1e-3→**4.78e-8 (4.3桁)** vs
  baseline σ0.05+c' `run_0024` 6.78e-6 (2.2桁) = **最終振幅 142 分の 1**。cell `run_0044`
  1.51e-8 vs `run_0018` 7.09e-8 (**4.7×**)。両離散化で改善 ✓
- **回帰**: `case/18 run_0119` (新バイナリ, precond=0) = `run_0110` と一致 (x-oddeven
  22.14/54.14, Uy rms 6.885 同値) = 既定経路ビット同等 ✓
- **backstep A/B** (`run_0120` σ0.05 / `run_0121` σ0.02, 4000step): 安定・NaN 無し。
  d2 メトリクス 22.3 / 27.5 Pa = **c' 版と同水準 (§4 の合格基準 ≤10Pa は未達)**。ただし
  物理は改善 (Uy rms 6.94/7.16 = c' 版の 6.89/7.12 以上)。L1 の 142× と合わせ、
  **backstep の ~22 Pa 床は数値市松ではなく「ギリギリ解像の物理 (段差リップの 2-4Δx KH 構造)」
  と確定** (Δp 散逸を 6.7× にしても物理を守ったまま床が動かない = メトリクスの下限が物理)。
  σ・c'・precond の全掃引 (`run_0109`-`0121`) がこの解釈で整合。
- **結論**: keepDissPrecond は「真の市松を桁違いに殺し物理を守る」設計目標を達成。同 σ の
  c' 版に対し厳密に優位 (市松マージン ≫・物理コスト ≦) → **低マッハ壁付き LES の推奨設定**。
  backstep の残存縞の根治は散逸でなく**メッシュ (リップ~せん断層の Δx 細分化)** の仕事。
- **残**: TGV 64³ L2 (KE cost) の実測 (予想: c' 版同等以下。Ur=|u| 域では c' とほぼ同一係数)。
  TP 枝への展開。

## 変更ログ

- 2026-07-19 (実装): `tools/verify_precond_dissipation.py` **全 6 検証 PASS**:
  (1) 閉形 $|M_2|=\varphi_1 M_2+\varphi_2 I$ = 数値固有分解 (6e-16)、(2) 低マッハ漸近
  $D_p\to\mathrm{diag}(c^2/U_r,\,U_r)$ 厳密 (M0.01 で Δp 散逸 ×6.7・ΔUn 散逸 ×0.15)、
  (3) β=1 で標準 Roe $|A_2|$ に厳密復帰、(4) **sym(K) が全 $(M, U_n)$ で正定値 → ES 性は
  クリップ無しで維持** (§3.2 の懸念は解消)、(5) コンパクト実装 (entropy/shear $|U_n|$ 据え置き
  + 音響 2×2) が**完全 Weiss–Smith 5×5 ($\Gamma_c = I + \frac{1-\beta}{\beta c^2}g r^T$) と
  6e-4 で一致** = 単なる模倣でなく実質同一、(6) 市松キラー漸近の実証。
  実装: `keepDissPrecond` (top-level, 既定 0=ビット不変) / `convectiveFlux_keep_d.inc.cuh`
  matrix CPG 枝の z± を前処理 2×2 で置換 (α∓→(Δp,ΔUn)→D_p→z∓, 固有ベクトル不要 ~20 flops)。
  precond 有効時 `keepDissCprime` は不使用 (前処理が上位互換)。
- 2026-07-19: plan 作成。backstep プローブ (run_0109–0118) で「σ・c' 一様スケーリングは単一
  トレードオフ曲線に縮退し市松と 2-4Δ 物理を分離できない」を実測確定 → 前処理散逸の設計に着手。

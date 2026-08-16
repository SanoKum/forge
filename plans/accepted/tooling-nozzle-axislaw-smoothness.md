# axis-Mach: 軸 Mach 則の滑らかさ比較 (A: knot / B: monotone B-spline / C: 非負 dν/dx B-spline)

## メタ

- **area**: `tooling / optimization`
- **status**: `done`
- **related_docs**:
  - `methods/design/overview.md` (§ 軸 Mach law、本計画で比較候補と数学的条件を追記)
  - `design/forge_design/geometry/axis_law.py` (A: `KnotQuinticAxisLaw`)
  - `design/forge_design/geometry/wall_axismach.py` (`wall_qa`, `AxisMachCFDWall`)
  - `design/forge_design/geometry/moc_inverse.py` (`fill_levels`, `cplus_flux_wall`, `inverse_design`)
- **related_plans**:
  - 親: [tooling-nozzle-axismach-chain.md](../accepted/tooling-nozzle-axismach-chain.md)
  - 前提: [tooling-nozzle-axismach-knot-law.md](../accepted/tooling-nozzle-axismach-knot-law.md) (A6 knot 則、本計画の比較基準)
- **created**: `2026-08-16`
- **owner**: `Claude (autonomous)`

## 1. 目的

現行 `KnotQuinticAxisLaw` (A) は knot で $M,M',M''$ が連続 (C²) だが、$M'''$ が
knot で不連続にジャンプする (基準条件 R=3, $L_c$=45, $M_d$=6 で
$M'''(K^-)\approx0.118$, $M'''(K^+)\approx-0.00060$)。見た目にも加速則の切替が急に見える。
これが**実害** (壁の fairness・特性線網の質) かを定量化し、より滑らかな軸則
(B: 内部 C³/C⁴ の monotone B-spline / C: 非負 $d\nu/dx$ B-spline) が実際に改善するかを
同一 $R=3,\ L_c=45$ で比較する。壁を後処理で平滑化するのではなく、**軸則自体を滑らかにして
から同じ逆 MOC に通す**。

**注記 (θ_w<μ_w の位置づけ)**: これは一般的な物理上の壁角上限ではない。本計画の逆 MOC は
「軸から C⁻ を後退させて壁に届かせる」構成であり、C⁻ が軸と上壁を単調・一対一に結ぶための
**位相条件**として扱う (「Mach 角が壁角の物理的上限」とは説明しない — knot-law plan の記述は
この位相条件を指しており、物理法則の主張ではない)。本計画では `min(μ_w-θ_w)` を**この逆 MOC
位相に対する余裕の補助指標**として報告し、health gate の根拠は §7 の characteristic
crossing / orientation flip 検査そのものに置く。

## 2. スコープ

- **やる**: A (診断追加のみ、数値不変) / B (C³/C⁴ monotone B-spline, $M(x)$ 直接表現) /
  C (非負 $d\nu/dx$ B-spline, Prandtl–Meyer angle 表現) の実装、同一逆 MOC・同一解像度
  (600/1200/2400) での比較診断 (軸則の滑らかさ・壁 fairness・特性線網トポロジ)、
  hard gate を満たす上位 2 案の CFD 評価。
- **やらない**: E 側 $M''=0$ の単純削除案 (次数を1落とす代替案、事前確認で $L_c=45$ の
  解決にならないと判明済み — §4.4)。knot 数を増やす A の拡張。壁表現側 (LSQ B-spline 等,
  A14 で既に棄却) の再訪。

## 3. 関連 docs と前提

- 軸 M 則の現行仕様: `methods/design/overview.md` § 軸 Mach law (A: knot 則)。本計画で
  §比較候補 (B, C) の数学的条件を追記する。
- 逆 MOC の網構造: 同 § 壁の決め方 (`fill_levels`, `cplus_flux_wall`)。
- 基準条件: `case/42.isobutane_wt/problem_ib_m6_R3_Lc45_MK2.5.yaml` 相当
  (semi-perfect, Hall throat anchor, `start_line: throat_char`, `wall_mode: cplus`,
  `exit_mode: characteristic`, `axis_dx0: 0.03`)。

## 4. 設計方針

### 4.1 A: 現行基準 (診断のみ追加、数値不変)

`KnotQuinticAxisLaw.deriv(x, order=3)` は既存の多項式微分ループが `order` を一般的に
扱うため**そのまま呼べる** (コード変更不要)。追加する診断:

- $\max|M'''|$ (knot 前後それぞれの区間で)
- knot 左右の $M'''(K^-), M'''(K^+)$ とジャンプ $\Delta M'''$
- $J_{\rm axis}=\int (M''')^2\,dx$ (区間ごとの解析多項式を Gauss 求積で厳密積分)

### 4.2 B: C³/C⁴ monotone B-spline ($M(x)$ 直接表現)

次数 $k=5$ (quintic) の clamped B-spline、単純内部ノット (多重度 1) で**構成的に C⁴**
(次数 $k$・単純ノットは $C^{k-1}$)。制御点 $c_0,\dots,c_{n-1}$ を未知数とする。

**ハード拘束** (線形、`BSpline` 基底関数の $x_A,x_E$ での値・1階・2階微分行を通じて構成):
$$M(x_A)=M_A,\quad M'(x_A)=M'_A,\quad M''(x_A)=M''_A,\quad M(x_E)=M_d,\quad M'(x_E)=0$$
$M''(x_E)=0$ は既定でハード拘束とするが、過拘束の疑いがあるため **hard / soft
(2 次ペナルティ, 重み $w_E$) の両方を実装し数値で比較する** (§4.4, §6)。

**単調性 (構成的、密サンプル依存にしない)**: B-spline の**制御多角形が非減少**
($c_0\le c_1\le\cdots\le c_{n-1}$) なら曲線自身が非減少 (variation-diminishing 性質の
帰結、標準的な shape-preserving spline の定理)。線形不等式 $c_{i+1}-c_i\ge0$ として
QP に埋め込む。

**平滑化目的**: $J_{\rm axis}=\int_{x_A}^{x_E}(M'''(x))^2\,dx$ を最小化。$M'''(x)$ は
次数 $k-3=2$ の区分多項式 (下位次数 B-spline) なので、$J_{\rm axis}=c^\top H c$
($H$ は 3 階微分基底の Gram 行列、各ノット区間を 6 点 Gauss–Legendre で厳密積分) の
凸 2 次形式。等式拘束 + 単調不等式のもとで `scipy.optimize.minimize` (SLSQP,
フォールバック trust-constr) で解く凸 QP。

**内部ノット**: 本数 $n_{\rm interior}$ と配置 (スロート側にクラスタ、指数
$\alpha=0.5$ で $s_i=(i/(n+1))^\alpha$) をパラメータ化。$n_{\rm interior}$ を増やして
(a) $L_c=45$ で実行可能 (単調・拘束充足) になる最小本数、(b) $J_{\rm axis}$ が頭打ちする
本数、の両方を確認し、**制御点数を不必要に増やさない** (§6 感度表)。

### 4.3 C: 非負 $d\nu/dx$ B-spline (Prandtl–Meyer angle 則)

$\nu(x)=\nu(M(x))$ (Prandtl–Meyer 角)、$q(x)=d\nu/dx$ を次数 $k_q=4$ (quartic, 単純ノットで
構成的 C³) の B-spline で表現する。**非負性は制御係数で構成的に保証**する
(B-spline は非負基底関数の凸結合 [partition of unity] なので、制御点 $c_i\ge0$ なら
$q(x)\ge0$ が全域で厳密に成り立つ — 密サンプルの偶然に依存しない)。

**Hall アンカーの chain rule 変換**: 軸則 A/B が直接使う $(M_A,M'_A,M''_A)$ を
$(q(x_A), q'(x_A))$ に変換する。$d\nu/dM$・$d^2\nu/dM^2$ は
`moc_kernel.dnu_dM(M, gas, order)` で計算する — CPG は閉形式
$d\nu/dM=\sqrt{M^2-1}/[M(1+\frac{\gamma-1}{2}M^2)]$ (Newton 反復で既に使っている式と同一)
を解析微分、semi-perfect は内部 $(M,\nu)$ テーブルを 3 次スプラインで微分
(数値差分ではなくスプライン微分、既存の `gas.nu`/`gas.mach_of_nu` テーブルと整合)。

$$q(x_A)=\left.\frac{d\nu}{dM}\right|_{M_A} M'_A,\qquad
q'(x_A)=\left.\frac{d^2\nu}{dM^2}\right|_{M_A}(M'_A)^2+\left.\frac{d\nu}{dM}\right|_{M_A}M''_A$$

**ハード拘束**: $\nu(x_A)=\nu(M_A)$ (基底行 nu=0)、$q(x_A)$・$q'(x_A)$ (上式、基底行
nu=0,1)、$q(x_E)=0$。$q'(x_E)=0$ は既定ハード (B と同様 hard/soft 両方実装)。
**追加の等式拘束** (積分値、基底関数を Gauss 求積で積分したベクトルとの内積):
$$\int_{x_A}^{x_E} q(x)\,dx=\nu(M_d)-\nu(M_A)$$
(端点条件だけでは自動的に成り立たない独立条件 — 明示的に QP へ加える。事後にも
`BSpline.antiderivative()` で再検証する)。

**平滑化目的**: $J_\nu=\int (q''(x))^2\,dx$ を最小化 (B と同じ Gram 行列パターン、
nu\_order=2)。自由度 ($n_{\rm interior}$) と正則化 (拘束の hard/soft 重み) への感度を
B と同様に確認する。

**$\nu(x)\to M(x)$**: $q$ の B-spline 反導関数 (`BSpline.antiderivative()`, scipy 標準機能
で検証済み — quad との差 <1e-9) から $\nu(x)=\nu(M_A)+F(x)-F(x_A)$、
$M(x)=\nu^{-1}(\nu(x))=$ `gas.mach_of_nu(nu(x))`。$M(x)$ は解析的に厳密だが、
$\nu\to M$ 反転が既存実装 (`mach_of_nu`, 6000 点テーブルの線形補間) を経由するため、
$M'(x),M''(x),M'''(x)$ の診断は**細かい一様格子上の $M(x)$ サンプルに 5 次補間
B-spline を当てて微分**する (生の非一様点への数値差分ではなく、`AxisMachCFDWall` の
壁曲率診断と同じ「補間 spline で微分」方式。解像度 2001/4001/8001 で収束確認)。

### 4.3b 実装で見つかった「前倒れ (front-loading)」問題と spread 拘束

B・C とも、端点条件 + 単調性/非負性 + 平滑化目的**だけ**で QP を解くと、最適解が
$x_E$ よりずっと手前 (実測 M6/R3/$L_c$45 で $x\approx20$、$L_c$ の半分未満) で $M_d$ に
到達し、残りを $M'''=0$ (または $q''=0$) の平坦域として済ませる。**平坦域は目的関数を
下げるため、最適化にむしろ好まれる**。これを逆 MOC に通すと、単一 quintic が起こした
fold と同じ病理 (曲げがスロート近傍に圧縮される) が再現され、`AxisMachCFDWall.validate`
の spline リンギングで検出された。対策として、中間点 $x_A+f_x L_c$ ($f_x$=0.5 既定) での
到達量に上限 $M_A+f_M(M_d-M_A)$ ($f_M$=0.75 既定、B は $M$、C は $\nu$ で) を課す
**spread 不等式拘束**を追加した (knot 則の $M_K$ と同じ発想を不等式として一般化)。

### 4.3c 特性線網の「seam」(良性境界) の発見

初期値線は「軸目標配列 (index 0..$n_{ax}$-1)」+「スロート特性線配列 (index $n_{ax}$..)」の
L 字構成 (`inverse_design`)。この 2 つの配列種別が切り替わる添字境界・レベルが小さい
(充填の初期数十段) 領域は、**健全な設計 (A) でも三角形の符号反転が局在する** (2026-08-16
実測、A で 96 個、全て seam 内)。物理座標 (throat から 1 $r_t$ 以内) でも同じ領域に
留まることを確認し (index 窓だけでは解像度を上げると窓からはみ出す)、
`characteristic_topology_diagnostics` は index 窓と物理座標窓の両方で seam を除いた
「内部反転数」を報告する。B の残存反転もこの同じ seam 領域に局在することを確認し、
A・B とも解像度 (n_axis≥1200) で内部反転 0 に収束することを確認した。

### 4.4 E 側 $M''=0$ 単純削除案 (主案にしない理由、事前確認)

- A 側 3 条件 + E 側 $M,M'=0$ の単一 4 次式: M6/R3 の単調性上限 $L_c\approx14.18$
  ($L_c=45$ を解決しない)。
- $M''(E)$ を自由にした単一 5 次でも単調上限 $L_c\approx24.1$ ($L_c=45$ に届かない)。
- knot 後段を 4 次→3 次に落とすと knot の $M'''$ ジャンプはほぼ不変、E に $M''$
  ジャンプが増える (問題を下流へ移すだけ)。

これらは B/C の代替にはならないため、B/C の soft-$M''(E)$ 比較 (§4.2, §4.3) の範囲でのみ
再訪する。

## 5. 実装ステップ

1. `design/forge_design/geometry/moc_kernel.py`: `dnu_dM(M, g, order)` (CPG 解析式 /
   semi-perfect スプライン微分、両方を有限差分でユニットテスト検証)。
2. `design/forge_design/geometry/axis_law_bspline.py` (新規): B-spline 共通ヘルパ
   (`_basis_matrix`, `_gram_matrix`, `_clamped_knots`) + `MonotoneBSplineAxisLaw` (B) +
   `NonnegDnuBSplineAxisLaw` (C)。`QuinticHermiteAxisLaw`/`KnotQuinticAxisLaw` と同じ
   `__call__`/`deriv`/`gates` インターフェース。
3. `design/forge_design/geometry/moc_diagnostics.py` (新規): 軸則の滑らかさ診断
   (共通インターフェース経由、A/B/C 全対応) / 特性線網トポロジ診断 (`fill_levels` の
   三角形 signed area・orientation flip・退化・隣接間隔ゼロ化・点堆積・非隣接 C±
   polyline 交差) / 壁曲率診断 (`AxisMachCFDWall._spl` 由来、deriv 3 まで)。
4. `design/forge_design/geometry/wall_axismach.py`: `AxisMachCFDWall.r()` の `deriv` 上限
   0..2 → 0..3 (壁 spline は次数 5 なので deriv=3 は無変更で評価可能、境界分岐に
   deriv=3 の分岐を追加するのみ)。
5. `design/forge_design/evaluate/runner_axismach.py::design_chain`: `geometry.axis_law` に
   `bspline_M` (B) / `bspline_dnu` (C) を追加。
6. `case/42.isobutane_wt/compare_axislaw_ABC.py` (新規): A/B/C × 解像度 3 段の一括比較
   ドライバ (JSON + PNG 生成)。
7. `design/tests/run_axislaw_bspline_tests.py` (新規): B/C の端点条件・連続性・単調性・
   非負 $d\nu/dx$・$\nu\leftrightarrow M$ 往復・CPG/semi-perfect・$x>x_E$ 一様性・不正入力・
   解像度依存・回帰 (A 不変)。既存 `run_axislaw_tests.py` / `run_inverse_tests.py` /
   `run_axismach_wall_tests.py` / `run_hall_tests.py` / `run_gas_tests.py` は無変更で
   ALL PASS を維持する。

## 6. 検証

- **単体**: `run_axislaw_bspline_tests.py` ALL PASS (端点条件・単調性/非負性・chain rule・
  ν↔M 往復・CPG/semi-perfect 両対応・$x>x_E$ 一様性・不正入力・A 回帰)。既存 5 スイート
  (`run_axislaw_tests` / `run_inverse_tests` / `run_axismach_wall_tests` / `run_hall_tests`
  / `run_gas_tests`) ALL PASS、回帰なし。
- **設計比較 (R=3, $L_c$=45, $M_d$=6, semi-perfect Tt=1550K)**: n_axis∈{600,1200,2400} で
  A/B/C を実行 (`case/42.isobutane_wt/compare_axislaw_ABC.py` →
  `compare_axislaw_ABC.json`)。

  | 案 | hard gate (600/1200/2400) | θ_max [°] | max\|M'''\| | max\|dκ/ds\| | min(μ_w−θ_w) [°] |
  | --- | --- | --- | --- | --- | --- |
  | A: knot | FAIL/PASS/PASS (600 は粗メッシュの点堆積のみ) | 14.11 | 0.2991 | 9.49 | 1.777 |
  | B: monotone B-spline | FAIL/PASS/PASS (同上) | 23.47 | **0.0604** | 10.47 | **0.089** |
  | C: 非負 dν/dx | FAIL/FAIL/FAIL (内部反転 216→282→436、解像度で悪化) | 21.42 | 6.569 | 128.7 | 4.784 |

  C は spread 拘束後も特性線網の内部反転 (seam 除外済) が解像度を上げるほど増え、
  壁 spline リンギング (~1°) も解消しない → **hard gate 不合格、CFD 対象外**。
  A・B は n_axis≥1200 で内部反転 0・wall QA/CFDWall 合格 (600 は粗メッシュの点堆積の
  みが原因で 1200 以降解消、fold ではない)。**B は $M'''$ を 5 倍改善するが θ_max が
  悪化し margin が 0.089° まで縮む** — 軸則の滑らかさが壁 fairness に直結しないことを
  示す。
- **CFD (A・B、procedures/ 準拠、新規 run)**: `run_0051_ib_m6_R3_Lc45_MK2.5_tp_inletfix`
  (A, 36000 step) / `run_0062_ib_m6_R3_Lc45_bsplineM_tp` (B, 12000 step)。両 run とも
  メッシュ品質 PASS・NaN 0。`check_convergence.py` は両方とも入口 BC 由来の warm 床で
  NOT CONVERGED (申し送り済み、[[axismach-chain-status]])。設計区間の軸 M は 8k→12k で
  ~3.6e-5 に凍結 (相対比較に有効)。‖ΔM‖∞: A 0.042% vs B **0.183%** (12k 換算で B は
  A の約 4.4 倍)、ε_M: A 0.0046% vs B 0.023%、ε_θ: A 0.0034° vs B 0.0134°。
  → **軸則の滑らかさ改善は CFD 精度を改善しない**、θ_max/margin の悪化がそのまま表れる。

## 7. 影響範囲

- 新規: `design/forge_design/geometry/axis_law_bspline.py`, `moc_diagnostics.py`,
  `design/tests/run_axislaw_bspline_tests.py`,
  `case/42.isobutane_wt/compare_axislaw_ABC.py`,
  `case/42.isobutane_wt/problem_ib_m6_R3_Lc45_bsplineM.yaml`。
- 変更: `moc_kernel.py` (+`dnu_dM`, CPG/`GasCPG`/semi-perfect 3 分岐)、`wall_axismach.py`
  (`AxisMachCFDWall.r` の `deriv` 上限 0..2→0..3)、`runner_axismach.py`
  (`geometry.axis_law: bspline_M`/`bspline_dnu` 追加)。
- 既存の `axis_law: quintic` / `axis_law: knot` 経路・生産 M4.2/M6 設定は無変更
  (回帰テストで確認済み)。**生産は A (knot) を維持、B/C は不採用**。
- docs: `methods/design/overview.md` に比較候補の数学的条件と結果を追記、
  `methods/index.md` は変更なし。`case/42.isobutane_wt/README.md` の run 一覧に
  `run_0062` を追加。

## 8. 完了条件

- [x] A/B/C の実装・単体テスト ALL PASS
- [x] 3 解像度での比較 (JSON + PNG) 完了、hard gate 判定・ranking 確定
- [x] 上位 2 案 (A・B) の CFD 評価 (VERDICT 明記、C は hard gate 不合格のため対象外)
- [x] `case/42.isobutane_wt/README.md` の run 一覧・成果物パス同期
- [x] status done → accepted へ移動、`plans/README.md` 同期

## 変更ログ

- 2026-08-16: 起票・実装 (A 診断/B/C/moc_diagnostics)・3 解像度比較・CFD (A・B)・結果
  ページ ([リンク](https://claude.ai/code/artifact/cf8a4c74-8d1f-44f4-a256-652cc122d00d))
  まで完了。**結論: A (knot) を生産として維持、B/C は不採用** (B は margin 薄く CFD
  精度も劣化、C は hard gate 不合格)。status done → accepted。

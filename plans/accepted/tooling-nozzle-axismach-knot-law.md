# axis-Mach A6: knot 付き区分 C² 軸 Mach 則 (高マッハ・長ノズル対応)

## メタ

- **area**: `tooling / optimization`
- **status**: `done`  <!-- 2026-08-16 起票・同日完了。M6 イソブタン (Tt 1550 K) で単一 quintic が全点 REJECT になったのを受けて着手 -->
- **related_docs**:
  - `methods/design/overview.md` (axis-Mach チェーン § 軸 Mach law / § 逆 MOC 初期値線)
- **related_plans**:
  - 親: [tooling-nozzle-axismach-chain.md](../accepted/tooling-nozzle-axismach-chain.md) (§5.6 で A6 を「長いノズルが要るとき再訪」と保留)
  - 前提: [tooling-nozzle-axismach-throat-characteristic.md](../accepted/tooling-nozzle-axismach-throat-characteristic.md) (初期値線 = スロート特性線)、[tooling-nozzle-semiperfect-gas.md](../accepted/tooling-nozzle-semiperfect-gas.md) (M6 ケースのガス)
- **created**: `2026-08-16`
- **owner**: `Claude (autonomous)`

## 1. 目的

出口 M6 (面積比 86.6、$\nu_d = 95.9°$) の風洞ノズルを axis-Mach チェーンで設計できるようにする。
単一 quintic Hermite の軸 M 則は単調性から $L_c \le 3.6\,\Delta M/M'_A$ (M6, R=3 で 17.9 $r_t$) に
縛られ、この長さでは M6 の壁が物理的に成立しない (§4.1)。knot 1 個の区分 $C^2$ 則で $L_c$ を
30–60 $r_t$ に伸ばし、壁角を壁マッハ角以下に保った設計を得る。

## 2. スコープ

- **やる**: `KnotQuinticAxisLaw` (2 区分 C²、DOF $L_c$・$M_K$) / runner の `geometry.axis_law: knot`
  / 逆 MOC の軸点列スロート側細分 (`axis_dx0`) / M6 (Tt 1550 K) の設計スイープ + TP CFD 評価。
- **やらない**: 単一 quintic の生産設定 (M4.2 以下、`axis_law: quintic` 既定) の変更。knot 2 個以上。
  $M_K$ の最適化 (スタディで振るだけ)。粘性 δ\* 補正の M6 適用 (別途)。

## 3. 関連 docs と前提

- 軸 M 則の現行仕様: `methods/design/overview.md` § 軸 Mach law: 5次 Hermite。
- 逆 MOC のレベル充填 (`fill_levels`) と C⁺ 流束閉包 (`cplus_flux_wall`): 同 § 壁の決め方。
- 前提: 初期値線 = スロート特性線 (`start_line: throat_char`)、壁 = `wall_mode: cplus`、出口 = 終端特性線。

## 4. 設計方針

### 4.1 M6 が単一 quintic で成立しない機構 (2026-08-16 実測)

R∈{2,3,4}, $L_c$∈{max,12,9}, n_axis 500–2500 の全点が壁 QA (spline リンギング 5–13° / κ₀R 16% 乖離 /
壁半径非単調) で REJECT。解像度を上げても壁の Δθ=12° の折れ (x≈0.52) は不変。原因は次の連鎖:

1. 終端特性線 (E→F) は一様域の上流境界なので $\mu_d$ の直線: $x_F - x_E = r_F/\tan\mu_d$。
   M6 では $9.3/\tan 9.6° = 55\,r_t$ (M4.2 は 16)。したがって**壁の膨張区間は $x_E - 55$ 付近で
   終わる**必要があり、$x_E \approx 18$ では膨張が**スロートから 2 $r_t$ 以内**に押し込められる。
2. その結果 $\theta_{\max} = 27.7°$ が壁マッハ角 $\mu_w = 24.7°$ ($M_w$ 2.4) を超える。
   逆 MOC (軸から C⁻ を後退させて壁へ届かせる構成) では $\theta_w \ge \mu_w$ で C⁻ の傾き
   $\tan(\theta-\mu)$ が正になり、後退線が壁近傍で折り返して**極限線 (fold)** を作る:
   C⁺ 線上に点が堆積し (x=0.76, r=0.99 に 20 点、M 1.73→1.93 の多価)、閉包壁に折れが出る。
3. R を上げると窓が広がる (R=20 で $L_c$ 35) が θ_max≈μ_w とぎりぎりで、根本解でない。

つまり**必要なのは長い $L_c$** であり、単一 quintic はスロート勾配 $M'_A$ を保ったまま長くできない
(長い側で $M'<0$)。序盤の急膨張と後段の緩い膨張を分ける自由度 = knot。

### 4.2 KnotQuinticAxisLaw

区分 1 $[x_A, x_K]$: 5 次 Hermite、始端 $(M_A, M'_A, M''_A)$ (Hall アンカー) → 終端 $(M_K, s_K, 0)$。
区分 2 $[x_K, x_E]$: 4 次 $M = M_K + \Delta M_2 (2u - 2u^3 + u^4)$、$u = (x-x_K)/L_2$。
$M' = 2\Delta M_2 (1-u)^2 (1+2u)/L_2 \ge 0$、$M'' = -12 \Delta M_2 u(1-u)/L_2^2 \le 0$ で
**構成的に単調・凹**、終端 $M' = M'' = 0$。knot 勾配 $s_K = 2\Delta M_2 / L_2$。
$L_1$ は台形則 $L_1 = 2\Delta M_1/(M'_A + s_K)$ (区分 1 の $M'$ が線形に落ちる長さ) で決め、
$s_K$ が $L_1$ に依存するので 1 次元 bisection。DOF は $(L_c, M_K)$ の 2 個、$L_c$ に上限なし。
ゲートは quintic と同じ ($M' \ge 0$ hard / $M''$ 反転 ≤ 1 品質、knot 則は構成上 1 回)。

### 4.3 軸点列のスロート側細分 (`axis_dx0`)

初期値線 (スロート特性線) は自身が C⁻ なので、その点は C⁻ 担体として退化する。初期値線と壁の
間の場は**軸点から後退する C⁻ だけで**埋まり、スロート直後の壁点は最初の数本の C⁻ でしか
決まらない。M6 は場が $x_{end} \approx L_c + 127\,r_t$ に伸び、n_axis 1200 でも dx 0.15 で
壁の最初の 0.3 $r_t$ が円弧から ±2e-4 $r_t$・±0.2° ジッタし QA を落とす。
軸点列を初項 dx₀=0.03・公比 1.05 の等比で細分し、等間隔値に達したら等間隔に切り替える
(全域等比だと下流が粗くなり出口平坦部に 1e-5 のジッタ → 単調 QA 落ち)。
M4.2 生産 (n_axis 500 等間隔) は無変更 (`axis_dx0` 未指定)。指定すると κ₀R 1.034→1.016 と改善。

## 5. 実装ステップ

1. `design/forge_design/geometry/axis_law.py`: `KnotQuinticAxisLaw` (評価/deriv/gates/admissible_Lc_range)。
2. `design/forge_design/geometry/moc_inverse.py`: `_axis_grid` + `inverse_design(axis_dx0=)`。
3. `design/forge_design/evaluate/runner_axismach.py::design_chain`: `geometry.axis_law` / `M_knot` / `axis_dx0`、
   `Lc_mode: max` は knot で拒否 (窓開放)。戻り値に `axis_law` / `M_knot` / `x_K`。
4. `design/tests/run_axislaw_tests.py` §8 (knot 端点・C²・単調・窓・不正入力・`_axis_grid`)。
5. `case/42.isobutane_wt/problem_ib_m6*.yaml` + 設計スイープ (`study_design_m6.json`) + TP CFD。

## 6. 検証

- **単体**: `run_axislaw_tests.py` ALL PASS (knot 端点誤差 1e-9、区分 2 単調凹、窓 (0.01, 400) が
  quintic 上限 17.9 を超える)。回帰: `run_inverse_tests` / `run_axismach_wall_tests` / `run_hall_tests`
  / `run_gas_tests` ALL PASS。
- **設計スイープ (M6, Tt 1550)**: R∈{2,3,5} × $L_c$∈{30,40,50,60} × $M_K$∈{2,2.5,3} の 36 点**全合格**
  (`study_design_m6.json`)。θ_max 11–21°、κ₀R 0.87–1.00、全長 5.3–6.9 m。
- **CFD (TP)**: `case/42.isobutane_wt/run_0039_ib_m6_R3_Lc45_MK2.5_tp` (L_c 45) — ‖ΔM‖∞ **0.041% $M_d$**、
  出口 ε_M 0.030% / ε_θ 0.004°、M_exit 5.9999、品質 PASS、NaN 0。残差は TP 入口 BC の床で
  NOT CONVERGED (既知の申し送り)、設計区間の軸 M は 8k→12k で 2.3e-3 動く → 継続 run で確認。
- **判定基準**: 軸 M ゲート 0.5% $M_d$ 以内、壁 QA 合格、出口 ε_M < 0.1%。

## 7. 影響範囲

- `design/forge_design/geometry/axis_law.py`, `moc_inverse.py`, `evaluate/runner_axismach.py`, `tests/run_axislaw_tests.py`
- 既存ケース: 既定 (`axis_law: quintic`, `axis_dx0` なし) は無変更。
- docs: `methods/design/overview.md` § 軸 Mach law に knot 則と細分を追記、`methods/index.md` 変更なし。

## 8. 完了条件

- [x] knot 則実装 + テスト
- [x] M6 設計スイープ全合格
- [x] M6 TP CFD 1 点で軸 M ゲート内
- [x] M6 R/L_c/M_K/L_U の CFD スタディ (11 点) と結果ページ・README 同期
- [x] status done → accepted へ移動

## 変更ログ

- 2026-08-16: 起票・実装・単体/設計スイープ/CFD 1 点まで。
- 2026-08-16 (同日): **M6 CFD スタディ完了** (`case/42.isobutane_wt/run_0051` (L_c 45, 36k step) / `run_0052–0061`,
  `study_cfd_m6.json`, [結果ページ](https://claude.ai/code/artifact/4d988b47-6cea-41ae-8dba-8646c19d6cbe))。
  R3 × L_c∈{30,40,45,50,60}、R∈{2,5} @L_c50、M_K∈{2,3}、L_U∈{8,16}: **全 11 点 0.03–0.08% $M_d$**
  (ε_M ≤0.017%, ε_θ ≤0.007°)。L_c 30→60 で θ_max 19.7→11.2°・dM 0.061→0.035%・全長 5.3→6.9 m; R2 最良 0.029%;
  M_K は鈍感 (2/2.5/3 で 0.040%); L_U8 は悪化 (0.078%)、16 は 12 と同等。**推奨 R 2–3 / L_c 45–50 / M_K 2.5 / L_U 12
  (全長 ≈ 6.0–6.3 m)**。全 run 品質 PASS・NaN 0、軸 M は 8k→12k で 2–5e-5 に凍結 (基準点 36k まで 3e-5)。
  途中、TP `inlet_Pressure` BC の静圧参照ブレンドが M6 (M_in 0.011) で指数成長発散 (run_0040) → rf 0.5→0 に修正
  (`boundaryCond_d.cu`; M4.2 残差床 3e-2→9e-7, run_0050)。設計基準点は当初 L_c 45 で走っており (run_0039/0051)、
  真の L_c 50 基準は run_0061 (0.054%)。status done → accepted。

# ノズル設計 axis-Mach チェーン (①風洞): Hall 遷音速 + 5次 Hermite 軸 Mach + 逆 MOC + CFD-in-the-loop

## メタ

- **area**: `tooling / optimization`
- **status**: `draft`  <!-- 2026-08-15 起票。実装未着手 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (現在仕様。実装着手時に本チェーン節を追加する)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) (§4.7 のアンカー機構「遷音速系列 + kernel MOC で $x_k$ に $M,M',M''$ を固定」の具体化にあたる)
  - 兄弟: [`tooling-nozzle-phase3-windtunnel.md`](tooling-nozzle-phase3-windtunnel.md) (現行 B8 生産チェーン。**壊さず回帰対照として維持**)
  - 兄弟: [`tooling-nozzle-walldriven-chain.md`](tooling-nozzle-walldriven-chain.md) (wall-driven チェーン。**W5 以降は保留** — 本計画を主線とするユーザ判断 2026-08-15。W1–W3 の成果物は流用)
  - 出典調査: [`notes/investigations/nozzle-throat-curvature-shape-representation-survey.md`](../../notes/investigations/nozzle-throat-curvature-shape-representation-survey.md)
  - 引き継ぎ: [`notes/sessions/2026-08-15-handover-wind-tunnel-nozzle-design.md`](../../notes/sessions/2026-08-15-handover-wind-tunnel-nozzle-design.md)
- **created**: `2026-08-15`
- **owner**: `sano`

## 1. 目的

ユーザ方針文書 (2026-08-15、「軸中心 Mach 数分布を用いた軸対称超音速ノズル設計方針」) を受け、
軸中心 Mach 分布 $M_{\rm axis}(x)$ を**低自由度の一次設計変数**とする Sivells 型逆設計チェーン

$$\text{Hall 遷音速解} \rightarrow M_{\rm axis} \text{ (5次 Hermite, 自由度 } L_c\text{)} \rightarrow \text{逆 MOC} \rightarrow \text{Euler CFD} \rightarrow \text{characteristic 追跡でアンカー更新}$$

の反復を、**既存 mode F 資産 (`moc_inverse` / `cminus_cfd` / `feedback`) の差し替え強化**として実装する。
完成時には①風洞 $M_d=4$ 同一要求で B8 チェーンと出口一様性・軸 M うねり振幅を直接比較でき、
軸 Mach law の指定→実現→補正が閉ループで回る状態を得る。粘性排除厚補正 (A7) までを視野に入れるが、
本計画の主要スコープは inviscid ループの確立 (A1–A5) とする。

## 2. スコープ

- **やる**
  - 壁上流部の構成確定と組立: **入口直管 + U→T 5次 Hermite** (W1 資産
    `wall_walldriven.UpstreamThroatPoly` を流用。上流円弧は使わない — ユーザ指定 2026-08-15)
  - 5次 Hermite 軸 Mach law (閉形式係数、設計ゲート、最小 $L_c$ 探索) の実装
  - Hall (1962) 遷音速解の実装 (CONTUR 実装形、Kliegel–Levine $S=R+1$ 置換) と Sauer との差し替え可能化
  - 逆 MOC への配線 (starting line 状態値の Hall 化、$E \to F$ 閉包、壁 QA 指標)
  - Euler CFD 検証と CFD-in-the-loop アンカー更新 ($x_A$, $M_A$, $M'_A$, $M''_A$) の反復ドライバ
  - $L_c$ による長さ制御 (外側 root finding、optional 制約)
  - (条件付き) 区分 5次 $C^2$ spline への拡張
- **やらない**
  - wall-driven チェーン W5 以降の継続 (保留。W4 厳密化の 2 経路は当該 plan に記録済み)
  - real-gas EOS・最大壁角/曲率などの制約付き最適化 (親 plan の将来項目)
  - RANS δ\* 固定点反復の完遂 (A7 着手時に後続 plan を起票して分離)
  - ①以外の機種 (②矩形風洞・③ベル 等) への展開

## 3. ユーザ方針の評価と既存資産との照合

原方針はゼロからの新規実装ではなく、**現行 mode F チェーンの弱点 3 つ (B5 / B10-c / 縦 starting line) に
それぞれ閉じた解を与える差し替え**として読める。照合結果:

| 原方針の要素 | 既存資産 | ギャップと判断 |
| --- | --- | --- |
| 5次 Hermite 軸 Mach law (6 条件をちょうど消費、自由度は $L_c$ のみ) | `geometry/bezier.py::MachBezier` (自由 CP 3 個 + 端点条件) | **B10-c の確定 failure「3 自由 CP では始端 C² 整合と出口一様性を両立できない」の構造的解決**。Hermite は両端 6 条件を構成的に満たし、過不足の自由度が残らない。差し替え採用 |
| Hall 遷音速解から $M_A, M'_A, M''_A$ を取得 | `geometry/transonic.py::SauerThroat` (1次のみ)。Hall/KL は未実装 (調査 §2.3 に文献情報) | **B5 で「Sauer の Cauchy データが CFD と乖離する」ことは実測済み**。Sauer 1次には $M''$ が構造的に存在しない問題 (調査 §9-2) も Hall 高次化で解消。CONTUR レポート (`papers/nozzle_design/`, ローカル) を実装ソースにする |
| CFD 場での characteristic 再追跡 → $x_A$・アンカー更新 (§14) | `geometry/cminus_cfd.py` (`extract_x_reach` ほか、B10/W0 で検証済み。W0 では $x_{\rm reach}$ での $M,M',M''$ 抽出まで実施) | **そのまま使える**。有効性の実測根拠は B6 (CFD アンカー化だけで junction wave 58% 減)。原方針 §15 の平滑化フィット (生の 2 階差分禁止) だけ新規実装 |
| 逆 MOC (軸 M 指定 → 場 → 壁流線) | `geometry/moc_inverse.py::InverseMOC` (検証済み。スロート直後の wedge 未計算域は既知の制限) | 再利用。starting line の**幾何は当面据え置き、状態値のみ Sauer→Hall** に置換する最小差分で始める |
| $x_E \ne x_F$ の明示区別、terminal Mach line 近似 $x_F - x_E \approx r_F\sqrt{M_d^2-1}$ | 現行実装は E/F を明示区別していない ($x_{\rm axis\_end}$ 止まり) | 新規。target 軸配列を $x > x_E$ で $M = M_d$ 一定として F まで延長し、壁が $\theta \to 0$ に漸近することを確認する |
| $L_c$ を主要長さ制御パラメータとする外側ループ (§22) | なし (現行は $x_{\rm axis\_end}$ を直接与える) | 新規 (小規模ユーティリティ)。$x_F(L_c)$ の root finding |
| δ\* 補正 (§23–25) | `feedback` v1 の δ\* 経験式 (CFD-free) 実装済み | A7 で流用。RANS δ\* 抽出との固定点反復は後続 plan |

**wall-driven チェーンとの関係**: W4 (wave-cancellation 壁生成) は平面近似止まりで軸対称厳密化が未達
(実 CFD データで march 発散)。ユーザ判断により軸 Mach 指定路線を主線化し、W5 以降を保留する。
W1–W3 で得た資産 (段階起動 runner・cone 検証・`WallDrivenCFDWall`) と W0 の計測手法
(`analyze_r5_ab.py` の trend-residual うねり指標) は本チェーンの検証にそのまま流用する。

**円弧こぶ問題との関係**: W0 で「R=2 円弧が軸 M こぶの支配要因」と確定している。本チェーンでは
2 段で対処する: (i) **スロート上流の壁は円弧でなく入口直管 + U→T 5次 Hermite** (W1 資産流用、
ユーザ指定 2026-08-15) とし、$\rho_t$ は端点条件 $r''(T) = 1/\rho_t$ として入る。**スロート T
から下流は自由設計** — 円弧を幾何として引き継がない (Hall/Sauer の $\rho_t$ はスロートを骨接円
[osculating circle] で局所近似する漸近解の**入力パラメータ**にすぎず、実際の壁がその先も円弧を
描くことを要求しない。T での曲率一致 $r''(T)=1/\rho_t$ で接続条件は完結する — ユーザ訂正
2026-08-15)、(ii) 残るこぶ成分は **CFD 反復でアンカー $x_A, M_A, M'_A, M''_A$ を実測値に更新して
吸収する** (B6 がこの経路の有効性を半減として実測済み)。将来 $\rho_t$ (= $R \cdot r_t$) を動かす
判断は outer loop の話として分離する。

## 4. 記号と repo 内対応

| 原方針の記号 | 定義 | repo 内の対応物 |
| --- | --- | --- |
| $x_t$, $r_t$, $\rho_t$ | スロート位置・半径・壁曲率半径 | `probdef` の throat 諸元。$\rho_t = R \cdot r_t$ (campaign は R=2) |
| $x_A$ | 軸 Mach 指定開始点 (MOC 制御の始点) | 初回 = 遷音速 starting line の軸上位置 (Hall 場、$M_{\rm start}\approx1.05$。親 plan §4.7 のアンカー点に相当)。CFD 反復後 = $x_{\rm reach,CFD}$ (`cminus_cfd.extract_x_reach`、B10 の origin と同一) |
| $M_A, M'_A, M''_A$ | $x_A$ での接続条件 | 初回 = Hall + kernel MOC 場から評価。反復後 = CFD 軸データの平滑化フィットから評価 |
| $x_E$ | 軸上で初めて $M = M_d$ に到達する点 ($M'(x_E)=M''(x_E)=0$) | 現行 `x_axis_end` に相当 (ただし現行は E/F 未区別) |
| $x_F$ | 物理出口 (断面全体が $M_d$, $\theta=0$) | 新規に明示。第一近似 $x_F - x_E \approx r_F\sqrt{M_d^2-1}$, $r_F = r_t\sqrt{A_e/A_t}$ |
| $L_c = x_E - x_A$ | 軸 Mach 加速区間長 = **主要設計自由度** | 新規 |

## 5. 設計方針

### 5.1 軸 Mach law: 5次 Hermite (閉形式)

$s = (x - x_A)/L_c \in [0,1]$、$M(s) = \sum_{i=0}^{5} a_i s^i$。両端 6 条件
($M_A, M'_A, M''_A$ / $M_d, 0, 0$) から閉形式:

$$a_0 = M_A,\quad a_1 = L_c M'_A,\quad a_2 = \tfrac12 L_c^2 M''_A,\quad \Delta M = M_d - M_A$$
$$a_3 = 10\Delta M - 6a_1 - 3a_2,\quad a_4 = -15\Delta M + 8a_1 + 3a_2,\quad a_5 = 6\Delta M - 3a_1 - a_2$$

(端点条件の充足は代数的に検算済み。単体テストでも機械精度で確認する。)

**設計ゲート** (原方針 §8, §28):

- hard: $M'(x) \ge 0$ 全域 (⇒ $M_A \le M \le M_d$、overshoot ゼロが自動保証)
- 品質指標: $M''$ の符号反転回数 $\le 1$ ($M'$ の一峰性)
- 端点 $C^2$ 条件は構成的に満足 (テストのみ)

**最小 $L_c$ 探索**: $M'\ge 0$ が破れる境界を bisection で求め、その接続条件で許される
$L_{c,\min}$ を返す (原方針 §10)。壁角・曲率による実用下限は MOC 後の QA で別途評価。

### 5.2 Hall 遷音速解

Hall (1962) の軸対称遷音速系列を、CONTUR (Sivells AEDC-TR-78-63、ローカル PDF) の実装形
(Kliegel–Levine の $S = R+1$ 置換込み) で `geometry/transonic.py` に追加する。
API は `SauerThroat` 互換 (`mach` / `theta` / `x_axis_of_mach` / `starting_line`) +
新規 `axis_anchor(x)` → $(M, M', M'')$ (解析微分または同一級数からの評価。生差分は使わない)。

**検証**: (i) 1次項が Sauer に退化すること (係数比較 + $R \to \infty$ 数値一致)、
(ii) Cuffel–Back–Massier (AIAA J 1969、ローカル PDF) の小 R 実験データ照合、
(iii) 既存 run のスロート近傍 CFD 場との直接比較 (B5 と同じ手法。Sauer 比で乖離が縮むこと)。

### 5.3 壁の全体構成、逆 MOC への配線と E→F 閉包

- **壁の全体構成** (上流→下流):
  **入口直管** ($r = r_U$) → **U→T 5次 Hermite** (`wall_walldriven.UpstreamThroatPoly` 流用:
  直管と $C^2$ 接続、$r'(T) = 0$, $r''(T) = 1/\rho_t$、適合条件 $\mu \le 20$ 等の検証済みロジックごと使う)
  → **T (スロート) から下流は完全に自由** — 円弧・その他の幾何を一切挟まない (ユーザ指摘
  2026-08-15: 「スロート点以降は自由に設計、接続部のみ曲率保持」)。T での曲率一致
  $r''(T)=1/\rho_t$ が唯一の接続条件であり、それは U→T Hermite の端点条件と Hall/Sauer 遷音速解の
  パラメータ $\rho_t$ が共有することで自動的に満たされる。T 以降の壁は**逆 MOC の出力**そのもの
  (下記) であり、独立した幾何 DOF は持たない。
- **starting line と T→A 間の扱い**: Hall/Sauer の遷音速解は「スロートを骨接円 $\rho_t$ で局所近似する
  漸近解」であり、$\rho_t$ はモデルの入力パラメータであって、実際の壁がその先も円弧を描くことを
  要求しない。T のすぐ下流に置く transonic starting line ($M_{\rm start}\approx1.05$ 程度、
  `SauerThroat.starting_line` / `HallThroat.starting_line` 互換 API) は**壁形状ではなく場の
  コーシーデータ**であり、その線上の壁点位置・状態も遷音速解の場そのものから評価する
  (T に近いほど U→T Hermite の曲率に収束する近似で、独立した幾何設計ではない)。
  この starting line から下流は、既存の `InverseMOC.fill` (軸 Mach law 駆動) がそのまま
  壁流線を出力する — これが「T 以降は自由設計」の実体であり、別途「kernel 領域の壁マーチ」を
  設けない。初回は starting line の**幾何を当面据え置き、線上の状態値評価だけ Sauer → Hall に
  差し替える**最小差分から始める (`InverseMOC.fill` の位相を触らない)。縦線構成が古典に
  対応物を持たない問題 (調査 §9-5) は認識しつつ、starting line の characteristic 化は
  本計画では扱わない (乖離が CFD 反復で吸収しきれない場合の将来項目)。
- **wedge 制限**: `moc_inverse` 既知の未計算 wedge (convex-hull 補間で壁流線が横断) は
  据え置き、壁 QA でその区間の曲率振動を監視する。
- **E→F 閉包**: target 軸配列を $x_E$ より下流へ $M = M_d$ 一定で
  $x_E + r_F\sqrt{M_d^2-1} + \text{margin}$ まで延長して fill し、壁流線が terminal Mach line
  を横切って $\theta_w \to 0$ に漸近することを確認して物理出口 $x_F$ を確定する。
  下流端で $|\theta_w|$ が閾値を下回らなければ延長不足/場の不整合としてエラーにする。
- **壁 QA** (原方針 §8.3, §28): 壁半径単調増加・最大壁角・$\kappa_w$ の振動/ジャンプ・
  characteristic crossing の有無を数値化して返すユーティリティを設ける。

### 5.4 CFD-in-the-loop アンカー更新

反復 $k$ で:

1. Euler CFD (node、段階起動) を実行し、`metrics/extract.axis_mach` で軸 M を抽出
2. `cminus_cfd` でスロート域からの代表 C⁻ を追跡し $x_A^{(k+1)} = x_{\rm reach,CFD}$ を確定
3. 軸データに**平滑化フィット** (局所 3–5 次多項式 or 平滑化 spline) を当て、
   **同一 fitted 関数から** $M_A, M'_A, M''_A$ を評価する (原方針 §15。生の 2 階差分は禁止。
   フィット窓/次数への感度を必ず併記する — W0 の trend-residual 手法を流用)
4. Hermite を再構築 ($x_E$ 固定 or $L_c$ 固定は run の目的で選ぶ) → MOC → 次の CFD

収束判定 (原方針 §18): $\max_x |M_{\rm CFD,axis} - M_{\rm target}|$、overshoot
$(M_{\max}-M_d)/M_d$、出口 $\sigma_M$・$\theta_{\rm exit}$、MOC/CFD characteristic 位置差
($x_{\rm reach}$ の予実差)。数値ゲートは §7 に定める。

軸 M 分布**全体**の誤差 $e^{(k)}(x)$ を設計自由度へ射影する補正 (原方針 §17) は、
アンカー更新のみで収束しない場合の A6 (knot 自由度追加) とセットで扱う。

### 5.5 長さ制御

$x_F \approx x_A + L_c + r_F\sqrt{M_d^2-1}$ を初期 guess とし、全長目標がある場合は
MOC の実測 $x_F(L_c)$ に対する root finding を外側ループとして提供する (optional。
①風洞の現行要求には全長制約がないため、既定は $L_c$ 直接指定 + $L_{c,\min}$ 報告)。

### 5.6 段階拡張 (条件付き)

単一 quintic + アンカー更新で §7 の軸 M 誤差ゲートを満たせない場合のみ、内部 knot 1 個の
区分 5次 $C^2$ spline へ拡張し、$M(x_m)$ / $M'(x_m)$ / knot 位置を追加自由度とする
(原方針 §11)。単一高次多項式 (7次/9次) へは上げない。

## 6. 実装ステップ

ユーザ方針 §27 の Phase 1–7 を A1–A7 に対応させる (B 系列・W 系列との衝突回避のため A 系列と呼ぶ)。

1. **A0 — methods 更新**: `methods/design/overview.md` に「axis-Mach チェーン」節を追加
   (§5 の設計方針を現在仕様として記述)。`methods/index.md` 同期。**実装着手前に行う** (開発フロー準拠)。
2. **A1 — 軸 Mach law** [Phase 1–2]: `design/forge_design/geometry/axis_law.py` 新設 —
   `QuinticHermiteAxisLaw` (閉形式係数・$M/M'/M''$ 評価・設計ゲート・$L_{c,\min}$ 探索)。
   テスト `design/tests/run_axislaw_tests.py` (端点条件機械精度・ゲート判定・$L_{c,\min}$ の性質)。
3. **A2 — Hall 遷音速** [Phase 1]: `geometry/transonic.py` に `HallThroat` 追加 (§5.2)。
   テストは Sauer 退化・文献照合・CFD 照合の 3 系統。
4. **A3 — 壁組立 + 逆 MOC 配線 + E→F 閉包 + 壁 QA** [Phase 3]: `moc_inverse.py` の starting line
   状態評価の抽象化 (Sauer/Hall 差し替え可能に)、target 延長と $x_F$ 確定、壁 QA ユーティリティ
   (`geometry/wall_qa.py` 等)。全域壁 (直管 + `UpstreamThroatPoly` + T 以降は逆 MOC 壁流線 + 平行部。
   円弧やその他の中間幾何を挟まない) を CFD 用壁 CSV へ出力する組立クラスを追加
   (`WallDrivenCFDWall` を雛型に)。放射源流照合・既存 mode F 回帰 (Sauer 指定で従来結果不変) を確認。
5. **A4 — Euler CFD 検証** [Phase 4]: `case/41.wind_tunnel_design/problem_m4_axismach.yaml` 新設、
   既存 `runner_wt` 流用 (壁 CSV の供給経路のみ追加)。**段階起動必須** (W3 の教訓)、node、
   run 命名 `run_00NN_axismach_<slug>` (既存連番に継続)、README run 表同期。
6. **A5 — アンカー更新ループ** [Phase 5]: `feedback/` に反復ドライバ追加 (v2 `euler_loop` の
   凍結マップ帰還とは別系統として実装し、共通部品を流用)。平滑化フィット + 感度併記。
   §7 のゲートで収束を宣言し、B8 と同条件比較。
7. **A6 — 区分 quintic $C^2$** [Phase 6]: A5 が軸 M 誤差ゲート未達の場合のみ着手 (§5.6)。
8. **A7 — 粘性補正** [Phase 7]: v1 δ\* 相関で壁補正 → RANS で $\delta^*_{\rm CFD}$ 抽出 →
   固定点反復。**着手時に後続 plan を起票して分離** (本計画は inviscid ループ完成で done にできる)。

## 7. 検証

- **単体 / ビルド**: A1–A3 の各テスト (§6)。ソルバ本体のコード変更はないためビルドは不変。
  既存 `design/tests/` の全テストが回帰しないこと (特に `moc_inverse` 系)。
- **検証ケース**: `case/41.wind_tunnel_design` (①風洞 $M_d=4$, $P_t=1$ MPa, $T_t=800$ K,
  $r_t=10$ mm, R=2 — **B8 と同一条件**)。run 投入は AGENTS 準拠:
  `check_mesh_quality.py` / NaN チェック / `check_convergence.py` / `check_quasisteady.py` の
  VERDICT を結果報告に添付。README run 表を同期。
- **判定基準**:
  - **設計ゲート** (MOC 前): $M' \ge 0$、overshoot 0、$M''$ 符号反転 $\le 1$ (§5.1)
  - **壁 QA** (MOC 後): 壁半径単調増加・characteristic crossing なし・$\kappa_w$ にジャンプなし、
    $x_F$ で $|\theta_w|$ 閾値以下 (§5.3)
  - **A5 収束** (Euler 後): $\max_x|M_{\rm CFD,axis}-M_{\rm target}| \le 0.5\%\,M_d$
    (phase3 と同じゲート)、CFD overshoot $< 1\%$ (開発段階基準。最終設計では厳格化)、
    出口 $\sigma_M$・$|\theta_{\rm exit}|$ を B8 の値と併記
  - **同条件比較**: 軸 M うねり振幅 (W0 の trend-residual 指標) を B8 (masked 0.0198) /
    Korte ベンチマーク (0.02) と比較し、劣後しないこと
  - **アンカー予実**: $x_{\rm reach}$ の MOC 予測と CFD 実測の差が反復で縮小すること

## 8. 影響範囲

- **コード**: `design/forge_design/geometry/` (新規 `axis_law.py`・`wall_qa.py`・壁組立クラス、変更
  `transonic.py`・`moc_inverse.py`。`wall_walldriven.py` の `UpstreamThroatPoly` は流用のみで
  既存挙動は変えない)、`design/forge_design/feedback/` (反復ドライバ追加)。
  ソルバ本体 (`solver_density_cuda/`) は不変。
- **既存チェーン**: B8 (phase3) は回帰対照として不変。Sauer 指定の既存経路は挙動不変を回帰テストで保証。
  wall-driven チェーンは W5 以降保留 (当該 plan に追記済み)。
- **ケース**: `case/41.wind_tunnel_design` に problem YAML と run を追加。README run 表同期。
- **docs**: `methods/design/overview.md` (A0)、`methods/index.md`、`plans/README.md`。

## 9. 完了条件

- [ ] A0: `methods/design/overview.md` に axis-Mach チェーン節を追加し `methods/index.md` を同期
- [ ] A1–A3 実装 + 単体テスト ALL PASS + 既存 mode F 回帰不変
- [ ] A4: Euler CFD 検証 run が品質・収束・準定常ゲートを通過 (VERDICT 添付)
- [ ] A5: アンカー更新反復が §7 のゲートで収束し、B8 との同条件比較を記録
- [ ] (条件付き A6 は必要になった時点で本計画に追記)
- [ ] status を `done` にし §10 に変更ログを記録、`plans/accepted/` へ移動
- [ ] `plans/README.md` の一覧を同期 (移動元・移動先の両セクション)

## 10. 変更ログ

- `2026-08-15` — 初稿。ユーザ方針文書 (軸中心 Mach 分布による設計方針) を取り込み、
  既存 mode F 資産との照合 (§3) と A0–A7 の実装ステップを定義。B10-c の failure
  (3 自由 CP で両端 C² 両立不能) を 5次 Hermite で構造的に解消する位置づけを明記。
  ユーザ判断により wall-driven チェーン W5 以降を保留し、本計画を主線とする。
- `2026-08-15` — ユーザ指定を反映: スロート上流の壁表現を**入口直管 + U→T 5次 Hermite**
  (W1 資産 `UpstreamThroatPoly` 流用) に確定。上流円弧は不採用、下流円弧は starting line
  壁足までの stub に縮退 (§5.3 壁の全体構成)。
- `2026-08-15` — ユーザ訂正を反映: 前項で入れた「T→starting line 壁足の円弧 stub」を撤回。
  スロート T から下流は完全に自由設計であり、接続条件は T での曲率一致
  $r''(T)=1/\rho_t$ のみ (U→T Hermite の端点条件と Hall/Sauer の $\rho_t$ パラメータが共有)。
  T 直後の transonic starting line は壁形状ではなく Hall 場のコーシーデータであり、
  そこから下流は既存 `InverseMOC.fill` の出力がそのまま壁になる。独立した「kernel 領域の
  壁マーチ」も設けない (§5.3 を再構成)。

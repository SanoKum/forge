# axis-Mach: D 案 = 端点アンカー + 内部補間点 1 点の C⁴ 区分 5 次軸 Mach 則

## メタ

- **area**: `tooling / optimization`
- **status**: `done`
- **related_docs**:
  - `methods/design/overview.md` (§ 軸 Mach law → D の数学的定義)
  - `design/forge_design/geometry/axis_law_onepoint.py` (本計画で新規)
  - `design/forge_design/geometry/moc_diagnostics.py` (トポロジ余裕・壁 fairness の窓分割を追加)
- **related_plans**:
  - 親: [tooling-nozzle-axismach-chain.md](../accepted/tooling-nozzle-axismach-chain.md)
  - 前提: [tooling-nozzle-axismach-knot-law.md](../accepted/tooling-nozzle-axismach-knot-law.md) (A = 比較基準、生産)、
    [tooling-nozzle-axislaw-smoothness.md](../accepted/tooling-nozzle-axislaw-smoothness.md) (B/C の教訓)
- **created**: `2026-08-16`
- **owner**: `Claude (autonomous)`

## 1. 目的

A/B/C 比較の教訓は「軸 $M(x)$ の数学的滑らかさを最小化しても、MOC 後の壁 fairness や
特性線網トポロジは改善しない」。そこで**自由度を増やさず、上下流アンカーを完全固定したまま、
内部補間点 1 点 $P=(\xi_P,\eta_P)$ だけで膨張配分を変える** D 案を実装し、$L_c$ も外側の
設計変数として continuation で扱い、**逆 MOC 後の特性線網・壁形状で** A と比較する。
現行 A の生産設定は変更しない。D が A より明確に良い場合だけ採用を検討する。

## 2. スコープ

- **やる**: `OnePointC4AxisLaw` (Bernstein 基底、4×4 線形解)、`gates()`、runner の
  `geometry.axis_law: onepoint`、`moc_diagnostics` の拡張 (トポロジ余裕・壁 fairness の
  throat 窓/下流窓分割・seam 物理分布)、$L_c$=45 の (ξ,η) 探索 + $L_c$ continuation
  (60→30)、最終候補の 1200/2400 再評価、A との比較表・PNG。CFD は設計段階で D が A を
  上回った場合のみ。
- **やらない**: 単一 6 次式 (主案にしない、事前確認で単調不可能)。内部点 2 点以上。
  既存 `quintic`/`knot`/`bspline_M`/`bspline_dnu` の数値変更。

## 3. 関連 docs と前提

- D の数学的定義: `methods/design/overview.md` § D。
- 基準条件: `case/42.isobutane_wt` R=3, $M_d$=6, semi-perfect Tt=1550 K, Hall anchor,
  `start_line: throat_char`, `wall_mode: cplus`, `exit_mode: characteristic`, `axis_dx0: 0.03`。
- 診断: `moc_diagnostics.py` (A/B/C 比較で導入)。

## 4. 設計方針

### 4.1 曲線 (methods 参照)

2 区間 $[x_A,x_P]$, $[x_P,x_E]$ の 5 次 Bernstein。A 側 3 + E 側 3 + P の $C^0$〜$C^4$ 5 +
$M(x_P)=M_P$ 1 = 12 条件で一意。Bernstein では $b_0,b_1,b_2$ (A 側)・$c_3,c_4,c_5$ (E 側)・
$b_5=c_0=M_P$ が直接決まり、$(b_3,b_4,c_1,c_2)$ を $C^1$〜$C^4$ の 4×4 線形系で解く。
微分は Bernstein の差分公式 ($\frac{d}{dt}\sum b_iB_{i,n}=n\sum(b_{i+1}-b_i)B_{i,n-1}$) を
再帰適用し区間長で割る (物理座標べき乗基底に戻さない)。

### 4.2 探索と continuation

- $L_c$=45 で $\xi_P\in[0.04,0.20]$、$\eta_P\in[0.15,0.80]$ の格子 → $M'\ge0$ で MOC 前に
  篩う → n_axis=600 で MOC → hard gate + 指標 → 順位付け。基準候補は A の knot 相当
  ($\xi_P\approx0.083$, $M_P\approx3.0$)。
- continuation: $L_c$=60 の健全解 → その (ξ,η) を初期値に 50, 45, 40, 35, 30 と短縮、各 $L_c$
  で局所再探索。$\ell_P=\xi_P L_c$、$\ell_P/R$ も追跡し、どれが一定かを事後判断。
- 最終候補は n_axis 1200/2400 で再評価し、A との順位・トポロジ診断が収束することを確認。

### 4.3 評価 (ranking)

1. 特性線網トポロジの健全性と余裕 (内部 flip 0、最小 signed area 比、最小隣接 spacing)
2. 壁半径単調・wall QA
3. $\min(\mu_w-\theta_w)$ の余裕 (暫定目安 ≥1° — この逆 MOC 構成の設計余裕、物理的
   壁角上限ではない)
4. throat 窓 ($0\le x\le0.25\,r_t$) **外**の壁 fairness ($\max|d\kappa/ds|$, $\int(d\kappa/ds)^2ds$)
5. $\theta_{\max}$
6. 軸則の滑らかさ ($M'''$、tie breaker)
7. $L_c$

## 5. 実装ステップ

1. `design/forge_design/geometry/axis_law_onepoint.py` (新規): `OnePointC4AxisLaw`。
2. `moc_diagnostics.py`: `wall_curvature_diagnostics(x_split=)` で throat 窓/下流窓を分けて返す、
   `characteristic_topology_diagnostics` にトポロジ余裕 (最小 signed area 比・seam 内反転の
   物理座標分布) を追加。
3. `runner_axismach.py`: `axis_law: onepoint` (`onepoint_xi_P`, `onepoint_eta_P`)。
4. `case/42.isobutane_wt/compare_axislaw_D.py`: 探索・continuation・最終評価・A 比較・PNG。
5. `design/tests/run_axislaw_onepoint_tests.py`。

## 6. 検証

- **単体**: `run_axislaw_onepoint_tests.py` ALL PASS (端点 6 条件 1e-15、P 通過、P で 0〜4 階連続
  2e-13、条件数 ≤2e4、単調 gate が事前計算表と一致、不正入力、無次元↔物理、CPG/semi-perfect、
  有限差分検算)。既存 6 スイート ALL PASS。
- **単調可能領域** (`compare_axislaw_D.py`, 帯探索 0.0025 刻み): (ξ,η) 平面上の**細い尾根**
  (η 幅 0.01–0.07)。L_c=45: ξ 0.02–0.18 に帯あり (ℓ_P=3.74 で M_P∈[2.92,3.11] — 事前表と一致)、
  L_c=60: ξ=0.03–0.05 (ℓ_P 1.8–3.0) にしか帯がない、L_c を短くすると帯が広がる。粗い η 格子
  (9×14) では 8/126 しか当たらないため帯を先に見つけてから候補を置く方式にした。
- **A/D 比較 (n_axis 600 探索 → 1200/2400 最終)**: A・D とも全解像度で hard gate 合格 (n=600 の
  raw 点堆積は seam 内のみ)。順位・flip 判定は解像度で逆転しない。

  | L_c | A: θ_max / μ−θ / min sin角 | D 最良 (M''反転3): θ_max / μ−θ / min sin角 | D 単峰 (M''反転≤1): θ_max / μ−θ |
  | --- | --- | --- | --- |
  | 60 | 11.2° / 3.45° / 0.062 | 22.5° / −0.62° / 0.030 | (帯なし) |
  | 50 | 13.0° / 2.41° / 0.051 | 18.1° / 5.48° / 0.026 | 23.4° / −2.41° |
  | 45 | 14.1° / 1.78° / 0.041 | 18.1° / 5.05° / 0.020 | 21.5° / −0.78° |
  | 40 | 15.6° / 1.05° / 0.033 | 17.9° / 4.52° / 0.013 | 22.6° / −1.91° |
  | 35 | 17.5° / 0.18° / 0.027 | 19.9° / 2.63° / 0.009 | 24.1° / −3.05° |
  | 30 | 19.7° / −0.90° / 0.023 | 16.6° / 0.87° / 0.016 | 22.4° / −2.24° |

  (@2400。下流 max|dκ/ds| は A 0.8–1.4 / D 0.4–1.6 で同等、max|M'''| は D 0.09–0.31 vs A 0.29–0.31)
- **D 最良の正体**: μ−θ で A を上回る候補は全て M'' 符号反転 3 = **内部プラトー** (L_c45 で
  x≈12–16 に M'≈0.007、M 3.9 で膨張が止まり x≈30 で再加速)。μ−θ が大きいのは壁がその区間で
  局所的に直線化するため。A と同じ品質条件 (M'' 反転 ≤1) を課すと **D は全 L_c で μ−θ<0、
  θ_max 21–27°** で A に劣る (前倒れ = B/C と同じ病理)。
- **トポロジ余裕**: 内部三角形の最小 sin 角 (signed Jacobian 相当) は D が A の約半分
  (0.020 vs 0.041 @L_c45/2400)、D の最小位置は x≈10–12, r≈3.5 に解像度収束 (A は位置が動く)。
- **摂動ロバスト性 (L_c45, n=600)**: M''_A ×0.5/1.0/1.5 で A は 3/3 合格 (μ−θ 不変 1.78°)、
  D 最良 (ξ,η 固定) は ×1.5 で M'<0 (帯の端にいる)。帯中央に再フィットすると μ−θ 2.10/1.06/0.60°。
  M'_A ±10% は A・D とも κ₀R ゲート落ち (Hall アンカーと R の整合を壊す摂動なので参考)。
- **内部点位置の法則**: ξ_P・ℓ_P/r_t・ℓ_P/R のいずれも L_c に対し一定にならない (D 最良 ℓ_P
  2.4→2.5→3.6→4.8→7.7→1.8; 単峰 D は ℓ_P≈1.3–2.1 r_t ≈ 0.4–0.7 R に留まる = 前倒れ側)。
- **CFD**: D は「トポロジ余裕が A 以上」を満たさず、同一品質条件下で fairness/長さでも A を上回らない
  ため**実行しない** (設計棄却)。

## 7. 影響範囲

新規: `axis_law_onepoint.py`, `run_axislaw_onepoint_tests.py`, `compare_axislaw_D.py`。
変更: `moc_diagnostics.py` (追加のみ)、`runner_axismach.py` (選択肢追加)、docs。
既存 axis_law 経路は無変更。

## 8. 完了条件

- [x] D 実装 + 単体テスト ALL PASS、既存 6 スイート回帰なし
- [x] $L_c$=45 探索 + continuation + 最終 3 解像度評価、A 比較表・PNG
- [x] 採否判定: **D 不採用 (設計段階で棄却、CFD 実行なし)、生産は A を維持**
- [x] status done → accepted、`plans/README.md` 同期

## 変更ログ

- 2026-08-16: 起票・実装 (`OnePointC4AxisLaw`, `moc_diagnostics` 拡張 [throat/下流窓・sin 角余裕・
  seam 分布・壁内側限定・点堆積の seam 除外], runner `axis_law: onepoint`)・探索/continuation/
  3 解像度評価・[結果ページ](https://claude.ai/code/artifact/b3142a8a-0c48-47d5-b7bc-f4653cc939c0)。
  **結論: 1 点補間 + C⁴ 拘束では、前倒れ (単峰側) か内部プラトー (μ−θ を稼ぐ側) のどちらかしか
  表現できず、A の knot 則 (序盤急・後段緩の単峰、C²) を超えない。D 不採用。** status done → accepted。

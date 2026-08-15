# axis-Mach チェーン A8: 初期値線をスロート特性線にする (+ 逆 MOC 充填のベクトル化)

## メタ

- **area**: `tooling / optimization`
- **status**: `done`  <!-- 2026-08-15 起票・同日完了 (A/B 検証まで) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (「axis-Mach チェーン」節)
- **related_plans**:
  - 親: [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md)
    (A0–A5 で確立した本チェーン。本計画はその starting line 構成の差し替え)
  - 兄弟: [`../active/tooling-nozzle-phase3-windtunnel.md`](tooling-nozzle-phase3-windtunnel.md) (B8 生産構成、回帰対照)
- **created**: `2026-08-15`
- **owner**: `sano`

## 1. 目的

axis-Mach チェーンの逆 MOC は初期値線 (Cauchy データ線) に **$M_{\rm start}=1.05$ の縦線**を
使っている。これは古典 (CONTUR/Sivells) に対応物がなく、(a) 設計が軸に効き始める位置が
スロートより大きく下流にずれる、(b) 特性線でないため単位過程の担体割当てが線上で反転し
**スロート直後に未計算の楔**が残る、(c) $[T, x_A]$ の壁を「Hall 模型が仮定する骨接放物線」で
埋める必要がある (幾何 DOF ではないと整理したが、設計壁ではない区間が実在する) という 3 つの
歪みを持ち込んでいた。

初期値線を**スロート壁点 $(0, r_t)$ から軸へ下る C⁻ (throat characteristic)** に差し替え、
**壁全体をスロートから MOC 出力にする**。あわせて逆 MOC 三角充填のスカラーループを
ベクトル化し、生産解像度の設計時間を実用域に戻す。

## 2. スコープ

- **やる**
  - `HallThroat.throat_characteristic()`: Hall 場の中で $(0, r_t)$ 発の C⁻ を軸まで積分
  - `inverse_design(..., start_line="throat_char"|"vertical")` (既定は当面 `vertical` = 回帰保護)
  - `AxisMachCFDWall`: 壁流線がスロート始点のとき骨接放物線区間を持たない構成を許す
  - `runner_axismach`: `geometry.start_line` の配線 (アンカー $x_A$ = C⁻ 軸着地点)
  - `InverseMOC.fill_arrays`: フロント一括ベクトル充填 (スカラー版と同一計算)
  - 単体テスト (特性線の性質 / スカラー-ベクトル等価 / 壁始点) と CFD 検証 run
- **やらない**
  - starting line そのものの高次化 (Hall 場の精度は据え置き)
  - `_CPlusMarch` (質量流束閉包) の完成 — 別課題として親計画の future work に残す
  - A7 (粘性 δ\*)

## 3. 関連 docs と前提

- 現在仕様: [`methods/design/overview.md`](../../methods/design/overview.md)「axis-Mach チェーン」
- 親計画 §10 の実測: 縦線構成での MOC 質量流束リーク (n_axis=500 で 0.79%) と
  生産解像度 `n_axis_inv: 2000` (逆 MOC が 1 パス 6.8 分の 90.8% を占める)
- 楔の既知問題と、棄却された「向き非依存 fill」の記録は `moc_inverse.InverseMOC.fill` docstring

## 4. 設計方針

### 4.1 スロート特性線 (throat characteristic)

曲率のあるスロートでは幾何スロートの壁は既に超音速 (Hall, $R=2$, $\gamma=1.4$ で $M_w=1.129$)
なので、そこから右進特性線 $dr/dx = \tan(\theta-\mu)$ を Hall 場の中で軸まで下ろせる。
これが古典的超音速 MOC の初期値線であり、その軸着地点 $x_A$ が「壁の設計が軸に影響を
及ぼせる最初の位置」になる。$r$ 等間隔に $n_{\rm start}$ 点へリサンプルして
`starting_line` と同じ規約 (軸→壁) で返す。

縦線との違い (実測、$R=2$, $M_d=4$):

| | 縦線 ($M_{\rm start}=1.05$) | スロート特性線 |
| --- | --- | --- |
| 軸着地 $x_A$ | 0.248 | 0.537 |
| 壁足 | $(0.248,\,1.0153)$、$\theta$ は Sauer 線形化の過小評価を手当てして上書き | $(0,\,1)$ 厳密、$\theta=0$ 厳密 |
| 壁流線の始点 | $x_A>0$ → $[T,x_A]$ を骨接放物線で埋める | スロートそのもの (全域 MOC 出力) |

軸着地が下流にずれるのは**設計対象区間が狭まったのではない**: 縦線構成では壁足から
下ろした C⁻ の軸着地が $x=1.81$ で、$[0.248, 1.81]$ の軸目標は**実際には壁で実現できない**
帯 (到達不能域) だった。特性線構成ではこの帯が構造的に消え、$x_A$ 以降の軸目標がすべて
壁の自由度で到達可能になる。

### 4.2 充填のベクトル化

三角充填は「レベルごとに隣接ペアを単位過程に通す」構造で、1 レベル内のペアは互いに独立。
`interior_vec` (`moc_kernel`) で予測子-修正子・軸対称源項・係数平均をそのまま配列化し、
`fill_arrays` はレベル数 $n$ 回の numpy 呼び出しで済ませる。スカラー版と**同一計算**
(等価性を単体テストで機械精度検証)。

## 5. 実装ステップ

1. **A8-1** `transonic.HallThroat.throat_characteristic` (済)
2. **A8-2** `moc_inverse.inverse_design(start_line=...)` と `wall_axismach.AxisMachCFDWall`
   のスロート始点対応 (済)
3. **A8-3** `moc_kernel.interior_vec` / `InverseMOC.fill_arrays` (済)
4. **A8-4** `runner_axismach` への `start_line` 配線 + 単体テスト追加
5. **A8-5** 設計 A/B (縦線 vs 特性線、同一解像度) — リーク・壁 QA・$x_F$・時間
6. **A8-6** Euler CFD 検証 run (新連番) + アンカー更新 1–2 パス、親計画 §10 の値と比較
7. **A8-7** methods 更新 → commit / push、本計画を `accepted/` へ

## 6. 検証

- **単体**: 既存 `design/tests/run_{inverse,hall,axislaw,axismach_wall}_tests.py` ALL PASS
  (回帰) + 新規: 特性線が Hall 場で $dr/dx=\tan(\theta-\mu)$ を満たす・壁足が $(0,1)$・
  軸着地の $M$ が Hall 軸値と一致・`fill_arrays` ≡ スカラー `fill`
- **設計 A/B**: 同一 $R$/$M_d$/解像度で mdot リーク・$r_F$ 誤差・壁 QA・$\theta_{\max}$
- **CFD**: `case/41.wind_tunnel_design` 新連番 run。AGENTS 準拠
  (`check_mesh_quality.py` / NaN / `check_convergence.py` / `check_quasisteady.py` の
  VERDICT 添付、README run 表同期)
- **判定**: $\max_x|M_{\rm CFD,axis}-M_{\rm target}| \le 0.5\% M_d$ を維持しつつ、
  (i) 到達不能域が消える (軸目標の全域が壁で実現可能)、(ii) 軸うねりが親計画の
  0.0019 を悪化させない、(iii) 出口 $\varepsilon_M$/$\varepsilon_\theta$ が悪化しない

## 7. 影響範囲

- **コード**: `design/forge_design/geometry/{transonic,moc_inverse,moc_kernel,wall_axismach}.py`、
  `design/forge_design/evaluate/runner_axismach.py`。ソルバ本体は不変。
- **既存チェーン**: `start_line` 既定を `vertical` に保つ限りモード F / B8 は挙動不変。
- **docs**: `methods/design/overview.md`、親計画 §10 への追記。

## 8. 完了条件

- [x] A8-4 配線 + 単体テスト ALL PASS (既存テストの回帰なし)
- [x] A8-5 設計 A/B の数値を §9 に記録
- [x] A8-6 CFD 検証 run で §6 の判定基準を確認、README run 表同期
- [x] A8-7 methods 更新・commit・`accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-08-15` — 起票。前セッションで実装済みだった 3 点 (特性線・start_line 配線・
  充填ベクトル化) を計画として遡って明文化し、検証項目を定義。
- `2026-08-15` — **A8 完了**。配線 (`runner_axismach` の `geometry.start_line`)・
  単体テスト追加・A/B 検証まで。

  **1. 単体テスト** (`design/.venv-opt/bin/python`, 全て ALL PASS):
  `run_hall_tests` に特性線 6 項目 ($R$=2/5 で 軸→壁の単調性・壁足が幾何スロート
  $(0,1)$ かつ $\theta$=0 [残差 $10^{-16}$]・壁足が超音速・軸着地の $M$ が Hall 軸値と
  一致・縦線より下流・線上で $dr/dx=\tan(\theta-\mu)$ [相対誤差 1.6–2.6e-4])、
  `run_inverse_tests` に `fill_arrays` ≡ スカラー `fill` (max|Δ| 1.4e-14) と
  特性線始点の逆設計 (壁が $(0,1)$ から始まる・単調・mdot 0.9946) を追加。
  既存テスト (`run_axislaw` / `run_axismach_wall` / モード F 回帰) は不変。

  **2. 設計 A/B** ($R$=2, $M_d$=4, $n_{\rm axis}$=2000/$n_{\rm start}$=121/
  `dx_wall` 0.005, 終端特性線出口):

  | | 縦線 | スロート特性線 |
  | --- | --- | --- |
  | $x_A$ (軸着地) | 0.2477 | 0.5369 |
  | 壁流線始点 | (0.2477, 1.0153) + 骨接放物線区間 | (0, 1) 厳密 |
  | $x_E$ / $x_F$ | 9.630 / 22.268 | 10.252 / 22.895 |
  | $r_F$ 誤差 | −0.317% | −0.258% |
  | MOC mdot 比 | 0.99858 | 0.99872 |

  $M''_A$ が特性線着地点では 0.0321 (縦線 0.1131) と小さく、$a_2=\frac12L_c^2M''_A$ の
  単調ゲート上限が緩むため**同じ単一 quintic で $x_E$ が 6.5% 長くなる** (A6 knot を
  必要とする閾値が上がる)。

  **3. CFD 検証** (`case/41.wind_tunnel_design/run_0049_axismach_a8v` [縦線・対照] /
  `run_0050_axismach_a8c` [特性線・本命]、いずれも Hall アンカー pass0 = アンカー更新なし。
  両 run とも メッシュ品質 **PASS** (AR 9.8–9.9 / skew 0.413)、`residual_history.csv` の
  全 rms 列に NaN/Inf **0**、`check_convergence.py --drop 2` **PASS (converged)** (全列 2.7–2.8 桁)、
  `check_quasisteady.py` **ALL STEADY**):

  | 指標 | 縦線 (run_0049) | 特性線 (run_0050) |
  | --- | --- | --- |
  | $\|\Delta M\|_\infty$ [% $M_d$] | 0.433 | **0.353** |
  | $\Delta M$ rms | 0.0096 | **0.0082** |
  | 軸うねり amp / p2p (deg3, 窓 0.55–2.5) | 0.00147 / 0.00251 | **0.00043 / 0.00109** |
  | 出口 $\varepsilon_M$ rms | 0.041% | **0.020%** |
  | 出口 $\|\varepsilon_\theta\|_{\max}$ | 0.0112° | **0.0062°** |
  | 出口コア面/全面 | 65/65 | 64/65 |
  | overshoot | −0.015% | +0.068% |
  | **到達不能帯** $x_{\rm reach,CFD}-x_A$ | **1.533** $r_t$ | **0.0019** $r_t$ |

  判定基準 (§6) はすべて充足: (i) 到達不能帯が消滅 (1.533 → 0.0019 $r_t$ = 設計の $x_A$ と
  CFD 実測の C⁻ 着地がほぼ厳密一致)、(ii) うねりは親計画の 0.0019 (3 パス) より小さい
  0.00043 を**アンカー更新なしの 1 パス**で達成、(iii) 出口 $\varepsilon_M$/$\varepsilon_\theta$ とも改善。
  **0.5% $M_d$ ゲートを 1 パスで通過** (親計画は 3 パス・0.451%)。解析 =
  `analyze_axismach_a8.py` → `axismach_a8_ab.json`。

  **4. 高速化** (数値は不変): 充填ベクトル化 (369 → 2.4 s) の後は **Delaunay 分割が支配的**に
  なっていた (n_axis=1000 で 1 個 8.4 s を計 6 個構築)。$(\theta,\nu)$ をまとめた補間器を
  1 個だけ作って壁流線・質量流束診断・終端特性線で共有する `field_interpolator` を導入し、
  **設計 1 回 155 s → 33 s** (合計で当初比 ~12 倍)。前後で設計壁がビット一致
  (max|Δ| は CSV 往復の $10^{-15}$ のみ) することを run_0050 の壁で確認済み。

  **5. 残件**: (a) $x_A$ 以降が全域到達可能になったので **CFD アンカー更新の位置づけが
  変わった** — 必須ではなくなり、「Hall 場と CFD 実測の差を吸収する任意の追い込み」になる
  (適用可否と効果は未測定)、(b) B8 生産構成の置き換え判断は未着手、(c) A7 (粘性 $\delta^*$)、
  (d) `_CPlusMarch` (質量流束閉包) は依然 future work — 残る MOC リーク 0.13% の恒久対処。

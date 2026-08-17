# 引き継ぎ: 加熱空気 (vitiated air) 風洞 M4.19 ノズルのパラメータスタディ

- **起票**: 2026-08-18 (別セッションへの引き継ぎ用。ユーザ指示「別セッションで設計・パラスタ」)
- **case**: `case/44.vitiated_air_wt/` (問題定義 YAML 11 点を用意済み。run はまだ 0 個)
- **チェーン**: axis-Mach (Hall 遷音速 → 軸 M 則 → 逆 MOC → node Euler CFD)。
  仕様は [`methods/design/overview.md`](../../methods/design/overview.md)、
  直近の同型案件は [`case/42.isobutane_wt/README.md`](../../case/42.isobutane_wt/README.md) (M4.2 イソブタン)。

---

## 1. 案件仕様 (ユーザ指定)

| 項目 | 値 |
| --- | --- |
| 全圧 $P_t$ | 1.14 MPa (1,140,000 Pa) |
| 全温 $T_t$ | 1060 K |
| 出口径 | **800 mm** (半径 0.400 m) |
| 目標出口マッハ数 $M_d$ | **4.19** |
| 組成 (モル分率) | H₂O 0.0462711 / N₂ 0.6989119 / O₂ 0.2091779 / Ar 0.008340143 / CO₂ 0.03729886 |
| 軸 M 則 | **単一 quintic** (`axis_law: quintic`)。knot 不要 — ユーザ判断、§3 で妥当性確認済み |

**未確定 (要ユーザ確認)**: **入口配管径**。指定が無かったので `geometry.r_inlet: 3.8`
(= 3.8 $r_t$ = φ800 mm、出口径と同じ) を仮置きしている。case/42 は「入口配管径 = 出口径」だったので
それに倣った。違う場合は全 YAML の `r_inlet` を差し替える (収縮の単調判定 $\mu=L_U^2/(R\Delta r)\le20$
は現状 $\mu$=6.4 と余裕があるので、φ600–1200 mm 程度なら $L_U$ 6 のままで通る)。

## 2. 事前計算済みの導出量 (再計算不要、YAML に反映済み)

質量分率 (`gas.species` に記入済み。$MW_{\rm mix}$ = 29.0806 g/mol、$R$ = 285.911 J/(kg·K)):

```
H2O 0.0286647 / N2 0.6732631 / O2 0.2301686 / AR 0.0114568 / CO2 0.0564467   (Σ=1.0000000)
```

semi-perfect (NASA-9/CEA, frozen 組成):

| 量 | 値 |
| --- | --- |
| $\gamma^*$ (スロート, $T^*$=912.6 K) | **1.32752** ← `gas.gamma` (Hall に渡る参照値) |
| $\gamma(T_t)$ / $\gamma$(出口 253.0 K) | 1.31616 / 1.39279 |
| $c_p(T_t)$ | 1190.2 J/(kg·K) ← `gas.cp` |
| $A/A^*(4.19)$ | **14.4430** |
| $r_t$ = 0.400/√(A/A\*) | **0.105252 m** (スロート φ210.50 mm) ← `spec.r_throat` |
| 出口 $T$ / $P$ | 253.0 K / **5235 Pa** ← `spec.p_ambient` |
| $\nu(M_d)$ / $\mu_d$ | 72.06° / 13.81° |
| $r_F/r_t$ / 終端特性線の戻り $r_F/\tan\mu_d$ | 3.8004 / 15.5 $r_t$ (=1.628 m) |

- **CPG では設計できない**: 一定 $\gamma$ だと $A/A^*$ が 17.14 ($\gamma^*$) / 18.10 ($\gamma(T_t)$) となり
  出口径を **+8.9 % / +11.9 %** 誤る。必ず `gas.model: semiperfect` のまま使う。
- 出口 253 K は NASA-9 の下限 200 K より上なので **`freeze_low_T` / `T_FLOOR` の細工は不要**
  (M6/Tt1550 の案件とはここが違う)。
- 全長の目安: 全長 ≈ ($x_F$ + $L_U$ + 0.5) $r_t$。$L_c$=8 で **3.20 m**、$L_c$=10 で 3.41 m。

## 3. 単一 quintic の妥当性と $L_c$ の窓 (確認済み)

単調ゲート ($M'\ge0$) の上限は $L_c$ ≈ 10.9 (R2) / 11.7 (R3) / 12.7 (R4) と広く、**knot は不要**
(ユーザ判断どおり)。一方、**下限は壁 QA が縛る**。二分法で確定した設計可能な $L_c$ 下限:

| R | $L_c$ 下限 (壁 QA) | 単調窓上限 | 備考 |
| --- | --- | --- | --- |
| 1.5 | **6.44** | 11.32 | |
| 2.0 | **6.50** | 10.88 | |
| 3.0 | **6.72** | 11.65 | |
| 4.0 / 5.0 | — | 12.68 / 13.71 | **設計不能** (壁半径非単調 / リンギング。case/42 M4.2 で R=5 が不能だったのと同傾向) |

$L_c$=6 の棄却理由は「壁 spline リンギング 1.74°」で、`n_axis_inv` 500→1500・`axis_dx0` 0.02–0.03 の
どれでも 1.74–1.75° と**解像度非依存** = 数値ノイズではなく**設計の自己不整合** (内部衝撃波の兆候。
case/41 で $L_c\le5.5$ に崖があったのと同じ性質)。**探索範囲は $L_c \in [7, 11]$、$R \in [1.5, 3]$ とすること。**

## 4. 用意済みの問題定義 (`case/44.vitiated_air_wt/`)

11 点。case/42 のスタディ構成 (R×L_c を $L_U$=6 で振り、$L_U$ は基準点で別途振る) を踏襲:

```
problem_va_R{1.5,2,3}_LU6_Lc{7,8,10}.yaml   … R × L_c の 9 点
problem_va_R2_LU{4,9}_Lc8.yaml              … L_U 感度 2 点 (基準 R2/L_c8)
```

設計のみ (`--prepare-only`) で全点が壁 QA 合格・メッシュ品質 PASS することは確認済み
(R2/L_U6/L_c8 の実測: θ_max 20.55°、κ₀R 1.000、$x_F$ 23.91 $r_t$、AR 10.5 / skew 0.43 PASS)。

## 5. やること (この引き継ぎのスコープ)

1. **設計スイープ**: 上記 11 点 (必要なら $L_c$ を 7→11 で細かく) を `design_chain` で回し、
   θ_max・κ₀R・$x_F$・全長を表にする。case/42 の `study_design*.json` 相当。
2. **CFD スタディ**: 各点を node Euler + TP (semi-perfect) で 12000 step。投入は
   ```
   cd design && ./.venv-opt/bin/python -m forge_design.evaluate.runner_axismach \
       ../case/44.vitiated_air_wt/problem_va_R2_LU6_Lc8.yaml ../case/44.vitiated_air_wt/run_0001_va_R2_LU6_Lc8
   ```
   ドライバは `case/42.isobutane_wt/run_study_m42_tp.py` を写して使う (run 番号は 0001 から)。
   評価指標: ‖ΔM‖∞/$M_d$ (0.5 % ゲート)、出口 $\varepsilon_M$ / $\varepsilon_\theta$、overshoot。
3. **推奨点の決定**と `case/44.vitiated_air_wt/README.md` の「## 計算 run 一覧」整備
   (**AGENTS.md 必須**: run を作ったら必ず同期)。
4. 必要なら形状エクスポート: `case/42.isobutane_wt/export_wall_m42.py` を写して
   点列 CSV + gmsh .geo (Euler 版 / NS 版) を出す。

## 6. 実行時の注意 (既知の罠、そのまま踏むと時間を失う)

- **段階起動必須**: `evaluate.cfl_main: 2.0` + `mid_stage: true` (soft 1 次 cfl0.5 → mid 2 次 cfl1 → 本段 cfl2)。
  TP は cfl4 で本段が爆発する ([[cutler-cpg-vs-tp-dplur-sst]]、case/42 run_0001)。
- **`thermoHrefTemp: 298.15` が必須** (runner が自動で入れる)。絶対基準 $h$ のままだと陰解法 Jacobian の
  $\chi_{\rm eos}=c^2-\kappa h$ が桁違いになり軸近傍で発散する
  ([tooling-nozzle-semiperfect-gas.md](../../plans/accepted/tooling-nozzle-semiperfect-gas.md))。
- **残差は落ちる**: TP 入口 BC の静圧参照ブレンドは 2026-08-16 に撤去済み (rf 0.5→0) なので、
  M4.2 級なら rms_ro が 1e-6 台まで落ちる。3e-2 でプラトーしたら**古い binary** を疑う
  (`solver_density_cuda/build/forge` を再ビルド)。
- **`check_convergence.py` / `check_mesh_quality.py` の VERDICT を報告に貼る** (AGENTS.md 必須)。
  残差プラトーでも軸 M が凍っていれば相対比較には使えるが、「収束した」とは書かない。
- **run ディレクトリは毎回新規連番**。既存を使い回さない。`res_*.h5` 等は commit しない
  (入力 config と `wall_design.csv` 等はリファレンスとして明示的に `git add -f` する運用)。
- `mesh.bndFirstOrder` は**使用禁止**。

## 7. 参照

- 設計チェーン仕様: [`methods/design/overview.md`](../../methods/design/overview.md)
- semi-perfect gas: [`plans/accepted/tooling-nozzle-semiperfect-gas.md`](../../plans/accepted/tooling-nozzle-semiperfect-gas.md)
- 単一 quintic / knot の使い分け: [`plans/accepted/tooling-nozzle-axismach-knot-law.md`](../../plans/accepted/tooling-nozzle-axismach-knot-law.md)
  (**本案件は M4.19 で knot 不要**。knot が要るのは M6 級で終端特性線の戻りが 55 $r_t$ に達する場合)
- 軸則の代替案 (B/C/D) は**すべて不採用**で決着済み:
  [axislaw-smoothness](../../plans/accepted/tooling-nozzle-axislaw-smoothness.md) /
  [axislaw-onepoint](../../plans/accepted/tooling-nozzle-axislaw-onepoint.md) → **quintic/knot を使う**
- 同型の先行案件 (手順・図・run 表の書き方): [`case/42.isobutane_wt/README.md`](../../case/42.isobutane_wt/README.md)
- 凝縮を後で足す場合: `evaluate.tp_species: split_h2o` + `evaluate.condensation`
  ([tp-split-h2o-condensation](../../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md))。
  本案件は H₂O 2.87 wt% (モル 4.6 %)、出口 253 K なので凝縮の可能性はある — 必要なら dry 完了後に。

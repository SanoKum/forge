# イソブタン風洞 M5 — R×L_U×L_c スイープ (ε_M vs 全長)

## メタ

- **area**: `その他 (design campaign)`
- **status**: `done`
- **related_docs**:
  - `design/CAPABILITIES.md` (対応状況の正本)
  - `case/42.isobutane_wt/README.md` (M4.2/M6 系列の既往知見)
- **related_plans**:
  - `plans/accepted/tooling-nozzle-axismach-knot-law.md` (M6 knot 則 — M5 でも必須と判明)
  - `plans/active/tooling-nozzle-semiperfect-gas.md`
- **created**: `2026-08-30`
- **owner**: `design-intake`

## 1. 目的

case/42 イソブタン燃焼ガス風洞の M_d=5 版で、軸 M 一様性 ε_M と全長 (質量代理) の
トレードオフカーブを R×L_U×L_c パラメータスイープで取得する。MOO は使わない
(axismach への opt/driver 配管は未検証のため、既往検証済みのスイープ+直接設計経路を採用)。

## 2. スコープ

- **やる**: 27 点マトリクス (R {1.5,2,3} × L_U {4,6,9} × L_c {14,17,20}) の dry TP CFD
  評価とトレードオフ集計。パレート勝者への凝縮後段確認 (Tsat_post 診断)。
- **やらない**: axismach への MOO 配管 / ループ内凝縮 CFD (必要になったら別 plan)。

## 3. インタビュー決定事項 (Q&A 要約, 2026-08-30)

| 論点 | 決定 | 根拠 |
| --- | --- | --- |
| 問題タイプ | `wind_tunnel_axisym_axismach` (✅) | 風洞ノズル軸 M 直接設計 |
| 目的 | ε_M 最小 × 全長最小 (2 軸トレードオフ) | ユーザ指定 |
| 探索方式 | パラメータスイープ | axismach MOO 配管は未検証、既往実績経路を優先 (ユーザ承認) |
| 固定仕様 | Pt 5 MPa / 出口径 1 m / 入口配管径 1 m / ガス = C4H10/空気 φ0.9 生成物 (M4.2 系流用) | ユーザ指定 |
| M_d | 5.0 | ユーザ指定 |
| Tt | **1550 K** (M6 系列と同値) | M5 出口静温 vs H2O 飽和線 (Murphy-Koop, x_H2O=0.1356): Tt1000K→余裕 −94 K (大幅凝縮)、1400K→−5 K、1550K→**+30 K**。dry 前提が成り立つ最小格の既往 Tt |
| 凝縮 | a. 後段確認 (dry 本体 + 勝者のみ Tsat_post →必要なら凝縮 ON 再評価) | ユーザ選択。`tp_species: split_h2o` で回すので凝縮 ON への切替は physProp/condensation 差し替えのみ |
| ガスモデル | semi-perfect TP (`thermo_href_temp: 298.15` = runner 既定) | Tt>600 K・多成分。CPG は出口径を数 % 誤る (M4.2 既往) |
| dv 範囲 | L_c ∈ [14, 20] (水準 14/17/20)、R {1.5,2,3}、L_U {4,6,9} | 下限 14 = **幾何成立限界の実測** (下記) 。上限 20 = ユーザ選定範囲の上端 |
| ゲート | check_convergence (--drop 2, paste IC の warm 系)、check_quasisteady、dM/ε_M 指標、メッシュ VERDICT | 既定。plateau 注意 (下記) |

## 4. 幾何成立プローブ (DOE 前必須チェック, 2026-08-30 実施)

`--prepare-only` (設計段のみ) で四隅+中心+境界を 20 点プローブした結果:

- **単一 quintic 軸 M 則は M5 で全点不成立** (壁 spline リンギング θ 乖離 >0.2°、
  かつ L_c 許容窓が M'≥0 単調ゲートで <14.7)。M6 と同様に **knot 則 (M_knot 2.5) +
  axis_dx0 0.03 が必須**。fold 閾値は M4.2 と M6 の間ではなく M5 より下にある。
- knot MK2.5: L_c 14/16/17/18/20 成立、L_c 12 不成立 (MK2.0/3.5 なら 12 も成立するが
  水準統一のため MK2.5 固定・L_c≥14 とする)。L_c ≤10 は MK 2.0–3.5 いずれも不成立
  (スロート曲率不整合 κ₀ ≫ 1/R、リンギング)。
- マトリクス四隅 (R1.5/L_U4/L_c14, R3/L_U9/L_c14, R1.5/L_U6/L_c20, R3/L_U12/L_c20 ほか) 全て成立。
- 全長 (超音速部, throat→exit): L_c14 → 3.60 m / 17 → 3.85 m / 20 → 4.08 m。

## 5. smoke 評価 (センター点)

`case/42.isobutane_wt/run_0072_ib_m5_R2_LU6_Lc17_tp/` (problem_ib_m5_R2_LU6_Lc17.yaml):

- メッシュ: VERDICT **PASS** (AR max 10.4, skew max 0.635)。
- NaN/Inf: なし (residual_history 36000 行 + res_12000.h5 の ro/P/T/roe)。物理値健全
  (T_min 307.9 K = 等エントロピー予測出口 308.3 K、P_min 6512 Pa ≈ 背圧 6540 Pa)。
- check_quasisteady: **ALL STEADY** (shock/machmax/pmax)。
- check_convergence: **NOT CONVERGED (warm 床)** — rms_ro fin 1.5e-6。これは M4.2 TP
  スタディ (run_0031–0038) / M6 スタディ (run_0051–0061, rms_ro 1e-6〜3e-5) と同じ
  **既知の warm 床パターン** (paste IC で init 残差が既に低く drop 桁が出ない)。場は
  STEADY・軸 M は凍結しており、系列内の相対比較には使える。キャンペーン集計時は
  同床レベル (~1e-6) に揃っていることを確認する。
- 指標: dM_max 0.083 %M_d (M4.2 ベスト 0.084% と同格)、出口 ε_M_rms 2.6e-4、
  M_core 5.0013、オーバーシュート 0.072%。

## 6. キャンペーン実行 (承認後)

27 点 = problem YAML 生成 (センター YAML を基底に R/L_U/L_c 置換) → runner_axismach
逐次実行。1 評価 ≈ 65 s 実測 (設計+メッシュ+CFD 12000 step, RTX 3060) → 27 点 ≈ 30 分,
ディスク ≈ 27 × 90 MB ≈ 2.5 GB。run 番号 0073–0099 帯。終了後: ε_M vs 全長の散布
(R/L_U で色分け) + パレートフロント抽出 → 勝者に Tsat_post 診断。

## 7. キャンペーン結果 (2026-08-30 実施)

27 点全て完走・健全 (NaN 0、全メッシュ PASS、残差床 8.5e-7〜3.2e-5 の既知 warm 床帯、
場の T_min 303.8–308.1 K)。集計 = `case/42.isobutane_wt/study_cfd_m5_tp.json`
(パレートタグ一覧 `study_cfd_m5_tp_pareto.json`)。run = `run_0072`〜`run_0098` (組合せは
ディレクトリ名)。全長は亜音速部込み (len_total = (L_U+L_pipe) r_t + x_F r_t)。

- **パレートフロント (全長 vs dM_max) 9 点**。主要 3 点:
  - 短尺側: `run_0076_ib_m5_R1.5_LU6_Lc14_tp` 4.133 m / dM 0.123 % / ε_M 0.014 %
  - **エルボ (推奨)**: `run_0096_ib_m5_R3_LU9_Lc14_tp` 4.366 m / dM 0.043 % / ε_M 0.017 %
  - 品質側: `run_0088_ib_m5_R2_LU9_Lc17_tp` 4.612 m / dM 0.035 % / ε_M 0.006 %
- **L_U=4 は M5 でも明確に悪い** (dM 0.6–0.9 %M_d、R3_LU4_Lc14 は 0.5 % ゲート超え)。
  L_U 6→9 で dM が 3 倍級改善 — M4.2/M6 と同傾向で、全長ペナルティは +0.24 m。
- L_c 14→20 の効果は L_U9 では小さい (dM 0.032–0.063 % 帯内)。短尺化は L_c でなく
  L_U とのバランスで決まる。
- 準定常: 勝者 3 点 check_quasisteady **ALL STEADY**。
- **凝縮後段確認 (方式 a)**: 勝者 3 点とも全場 min(T−Tsat_post) = **+29 K**、S_max 0.16
  → 凝縮リスクなし、dry 結果で確定。凝縮 ON 再評価は不要。

## 変更ログ

- 2026-08-30: インタビュー完了、幾何プローブ、smoke run_0072 PASS 相当 (plateau 注記付き)。
- 2026-08-30: 27 点スイープ完走・パレート抽出・Tsat_post 診断まで完了。status done、
  推奨 = R3/L_U9/L_c14 (エルボ)。
- 2026-08-30 (追補): 推奨点の物理壁を NS/SST y+~1 で直接検証 (run_0099–0102)。相関 δ\* (v1) は
  出口コア −0.71 % でゲート外 → CFD 抽出 δ\* v3 (δ\* plan の確立フロー、M5 は x_lo=12) で
  **−0.29 %・ε_M 0.10 %・固定点 1.001** に収束。実機壁の正本 = `points_m5_best_ns.csv` (v3、
  run_0102 と 4 μm 一致)。教訓: 生抽出 CSV の全域採用 (run_0101) は δ\* plan 記録済みの失敗の再演。

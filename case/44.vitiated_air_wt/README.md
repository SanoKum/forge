# case/44.vitiated_air_wt — 加熱空気 (vitiated air) 風洞 M4.19 (semi-perfect, axis-Mach チェーン)

案件 (2026-08-18): $P_t$=1.14 MPa, $T_t$=1060 K, **出口径 1600 mm** (半径 0.800 m), $M_d$=**4.19**。入口配管も φ1600 mm。
組成 (モル分率) H₂O 0.0462711 / N₂ 0.6989119 / O₂ 0.2091779 / Ar 0.008340143 / CO₂ 0.03729886
→ 質量分率 H2O 0.0286647 / N2 0.6732631 / O2 0.2301686 / AR 0.0114568 / CO2 0.0564467
($MW_{\rm mix}$ 29.0806 g/mol, $R$ 285.911 J/kgK)。$r_t$=**0.210505 m** (スロート φ421.01 mm、$A/A^*$=14.4430)。
設計 (MOC) は NASA-9/CEA semi-perfect (frozen、$\gamma^*$ 1.328 / 出口 1.393)、CFD は forge の TP 単一擬似種
(`thermoHrefTemp: 298.15`)。**CPG では出口径を 9–12 % 誤るので使わない**。

**軸 M 則は単一 quintic** (knot 不要 — 単調窓上限 $L_c$≈10.9–11.7 に対し必要長は 7–10)。
設計可能 $L_c$ は **下限 6.4–6.7 (壁 QA)** 〜 **上限 10.9–11.7 (単調ゲート)**、R は **1.5–3** (R≥4 は設計不能)。
入口配管径はユーザ指定 (2026-08-18 訂正) で**出口径と同じ φ1600 mm** → `geometry.r_inlet: 3.8` (=3.8004 $r_t$)。
引き継ぎ文書 (導出量の根拠・既知の罠): [`notes/sessions/vitiated-air-m419-nozzle-session-prompt.md`](../../notes/sessions/vitiated-air-m419-nozzle-session-prompt.md)。
チェーン仕様: [`methods/design/overview.md`](../../methods/design/overview.md)、同型の先行案件: [`case/42.isobutane_wt/`](../42.isobutane_wt/README.md)。

## 結論 (2026-08-18): 推奨 **R=2 / $L_U$=6 / $L_c$=8** (`run_0005_va_R2_LU6_Lc8`)、同等次点 $L_c$=9 (`run_0013`)

- **全 15 点が 0.5 % ゲート内** (‖ΔM‖∞/$M_d$ 0.047–0.315 %)。**R2 が全 $L_c$ で最良系列**、最良は R2/$L_c$9 の
  **0.047 %** (ε_M 0.008 %)、次いで R2/$L_c$8 **0.061 %** (ε_M 0.009 %)、R3/$L_c$8 0.070 %。差 0.01–0.02 % は
  設計点間の有意差とは言い難く (軸 M の 8k→12k 変動 ~2e-6 相対よりは大きいが、入口 BC 床・出口一様性指標のばらつきと同程度)、
  **R2 の $L_c$ 8–9 を同等ベスト**と読む。全長で 0.21 m 短い **$L_c$=8 を生産基準**とし、余裕が欲しければ $L_c$=9。
- **$L_c$=7 は崖の縁**: R2 0.165 %、R3 **0.315 %** (x≈4 $r_t$ に −0.3 % の大きな谷 = 内部圧縮の兆候。壁 QA 下限 6.5–6.7 の直上)。
  **$L_c$≥8 を使うこと**。$L_c$ 10–11 は ΔM は 0.08–0.13 % に微増するが出口流れ角 ε_θ が 0.004–0.006° と最良 (長さ 6.8–7.0 m)。
- **$L_U$=4 は不可** (0.213 %、ε_M 0.058 %、overshoot +0.36 %: スロート直後 x<1 $r_t$ に +0.2 % の山 = 収縮が短く Hall アンカーが崩れる)。
  $L_U$=9 は $L_U$=6 と同等 (0.072 %) で長いだけ。**$L_U$=6 でよい**。
- R1.5 は 0.09–0.13 % で $L_c$ 依存が単調でなく (κ₀R 0.92–0.94 とスロート曲率が MOC 要求より 6–8 % 大きい)、R3 は $L_c$≥9 で
  R2 より一貫して悪い。case/42 (M4.2 イソブタン、R2/L_c8 が最良) と同じ傾向。

**推奨点の実寸 (R2/$L_U$6/$L_c$8, Euler 版)**: スロート φ421.0 mm @x=0、収縮開始 x=−1.263 m (=−6 $r_t$)、入口直管端 x=−1.368 m
(φ1600)、**出口 x=+5.032 m (φ1596.9 mm; 1D 理論比 −0.19 %)**、**全長 6.40 m** (膨張部 5.03 m)、最大壁角 20.54° @x=0.255 m。
NS 版 (δ\* オフセット、実機壁): スロート φ425.4 mm @x=−0.97 mm、出口 φ1646.9 mm (δ\* 出口 25 mm)。
CFD 実測 (Euler): 出口軸 M 4.1898、コア M 4.1893、ε_M 0.009 %、ε_θ 0.018°、overshoot +0.14 %。

## スタディ結果 (node Euler TP semi-perfect, 12000 step, 全 run メッシュ PASS・NaN 0)

`study_cfd_va.json` / `study_design_va.json` (設計のみ 21 点: $L_c$∈{7..11}×R∈{1.5,2,3} + $L_U$∈{4,9}; **R2/$L_c$11 のみ単調窓外で REJECT**、他 20 点全合格)。
図: `study_cfd_va.png` (指標 vs $L_c$)、`study_axisM_va.png` (軸 M 誤差分布)。生成: `plot_study_va.py`。
**軸 M (目標 vs CFD) と Mach コンタ**: `figs/` (`plot_mach_va.py` で生成: `axis_mach_all.png` 15 パネル、`axis_mach_overlay.png`、`mach_contour_all.png` 全 run 一覧、`mach_contour_run_0005_*.png` / `_run_0013_*.png` 推奨点)、[閲覧ページ](https://claude.ai/code/artifact/b4158f43-0a74-4115-a3fb-282e22ec5d10)。

| run | R | L_U | L_c | θ_max° | κ₀R | x_F [r_t] | 全長 [m] | ‖ΔM‖∞/M_d | ε_M | ε_θ° | over | M_exit |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `run_0001_va_R1.5_LU6_Lc7` | 1.5 | 6 | 7 | 21.01 | 0.936 | 22.94 | 6.196 | 0.098 % | 0.025 % | 0.0162 | +0.161 % | 4.1880 |
| `run_0002_va_R1.5_LU6_Lc8` | 1.5 | 6 | 8 | 20.60 | 0.925 | 23.97 | 6.414 | 0.132 % | 0.018 % | 0.0217 | +0.149 % | 4.1861 |
| `run_0012_va_R1.5_LU6_Lc9` | 1.5 | 6 | 9 | 20.26 | 0.920 | 24.96 | 6.622 | 0.115 % | 0.013 % | 0.0218 | +0.162 % | 4.1876 |
| `run_0003_va_R1.5_LU6_Lc10` | 1.5 | 6 | 10 | 20.02 | 0.919 | 25.98 | 6.837 | 0.091 % | 0.016 % | 0.0143 | +0.155 % | 4.1909 |
| `run_0010_va_R2_LU4_Lc8` | 2 | 4 | 8 | 20.55 | 1.000 | 23.91 | 5.980 | 0.213 % | 0.058 % | 0.0238 | +0.361 % | 4.1970 |
| `run_0004_va_R2_LU6_Lc7` | 2 | 6 | 7 | 20.97 | 1.012 | 22.87 | 6.183 | 0.165 % | 0.016 % | 0.0164 | +0.151 % | 4.1905 |
| **`run_0005_va_R2_LU6_Lc8`** | 2 | 6 | 8 | 20.55 | 1.000 | 23.91 | **6.401** | **0.061 %** | 0.009 % | 0.0182 | +0.140 % | 4.1898 |
| **`run_0013_va_R2_LU6_Lc9`** | 2 | 6 | 9 | 20.21 | 0.996 | 24.89 | 6.608 | **0.047 %** | 0.008 % | 0.0156 | +0.153 % | 4.1903 |
| `run_0006_va_R2_LU6_Lc10` | 2 | 6 | 10 | 20.00 | 0.995 | 25.91 | 6.823 | 0.076 % | 0.011 % | 0.0059 | +0.136 % | 4.1947 |
| `run_0011_va_R2_LU9_Lc8` | 2 | 9 | 8 | 20.55 | 1.000 | 23.91 | 7.032 | 0.072 % | 0.013 % | 0.0171 | +0.098 % | 4.1882 |
| `run_0007_va_R3_LU6_Lc7` | 3 | 6 | 7 | 20.76 | 1.044 | 22.79 | 6.165 | 0.315 % | 0.023 % | 0.0164 | +0.182 % | 4.1941 |
| `run_0008_va_R3_LU6_Lc8` | 3 | 6 | 8 | 20.23 | 1.032 | 23.82 | 6.383 | 0.070 % | 0.020 % | 0.0164 | +0.163 % | 4.1927 |
| `run_0014_va_R3_LU6_Lc9` | 3 | 6 | 9 | 19.79 | 1.027 | 24.81 | 6.591 | 0.106 % | 0.021 % | 0.0073 | +0.170 % | 4.1963 |
| `run_0009_va_R3_LU6_Lc10` | 3 | 6 | 10 | 19.46 | 1.025 | 25.83 | 6.806 | 0.128 % | 0.021 % | 0.0041 | +0.124 % | 4.1953 |
| `run_0015_va_R3_LU6_Lc11` | 3 | 6 | 11 | 19.29 | 1.026 | 26.80 | 7.009 | 0.134 % | 0.022 % | 0.0054 | +0.068 % | 4.1931 |

**収束・定常の状態 (全 run 共通)**: `check_convergence.py` は **NOT CONVERGED** (本段開始時点で soft/mid 段により rms_ro
≈5e-7 の warm 床に達しており、本段 12000 step は plateau。rms_roUx 3e-4 / rms_roUy 4–5e-4 / rms_roe 0.4–0.5 [絶対値])。
「収束した」とは書かない。ただし**設計区間の軸 M は 4k→8k→12k で |ΔM|max 5e-6〜1.4e-5 (相対 ~2e-6) に凍結**、
`check_quasisteady.py --quantity machmax,pmax` は run_0005/0013/0008 で **STEADY**。設計点間の比較には十分。
最終場は全 run で P∈[5.19 k, 1.139 M] Pa、T∈[252, 1060] K、ρ>0、NaN/Inf 0 (residual 全列・`res_12000.h5` とも)。

## 推奨形状の数式・点列・gmsh (2026-08-18)

再生成: `design/.venv-opt/bin/python case/44.vitiated_air_wt/export_wall_va.py [--problem problem_va_R2_LU6_Lc8.yaml --tag va_best]`
→ `points_va_best_{euler,ns}.csv` (1028 点 [m] + 無次元)、`nozzle_va_best_{euler,ns}.geo` (gmsh 軸対称半平面、構造化 quad 365×64、
physID inlet 1/outlet 2/wall 3/axis 4/fluid 5)。$L_c$=9 版も `points_va_R2_LU6_Lc9_*.csv` / `nozzle_va_R2_LU6_Lc9_*.geo` に出してある。
**Euler 版** = CFD が実際に解いた非粘性設計壁、**NS 版** = 非粘性壁 + 排除厚 δ\*(s) の法線オフセット (A13 `PhysicalNozzleWall`、実機の壁)。

$r_t$=0.210505 m、$R$=2 $r_t$=0.42101 m、$L_U$=6 $r_t$=1.26303 m、$L_c$=8 $r_t$、$M_d$=4.19。無次元 $\hat x=x/r_t$, $\hat r=r/r_t$:

1. **入口直管** $-6.5\le\hat x\le-6$: $\hat r=3.8$ (φ1600 mm)。
2. **収縮 U→T (5 次 Hermite)** $-6\le\hat x\le0$、$\xi=(\hat x+6)/6$:
   $$\hat r(\xi)=3.8-19\,\xi^3+24\,\xi^4-7.8\,\xi^5$$
   ($\hat r(0)=3.8,\ \hat r'=\hat r''=0$ / $\hat r(1)=1,\ \hat r'=0,\ \hat r''=1/R=0.5$ を厳密に満たす。単調収縮判定 $\mu=L_U^2/(R\Delta r)$=6.43≤20)。
3. **スロート下流 (逆 MOC 設計壁)** $0\le\hat x\le23.906$: 壁テーブル `run_0005_va_R2_LU6_Lc8/wall_design.csv` (132 点 [x, r, θ, M_wall])
   を端条件クランプ付き補間 5 次 B-spline で表現 (スロート端 $\hat r'=0,\ \hat r''=1/R$)。閉形式は無いので点列 CSV を使う。

| | Euler 版 | NS 版 (δ\* オフセット) |
| --- | --- | --- |
| スロート | x=0, r=0.210505 (φ0.42101) | x=−0.000967, r=0.212680 (φ0.42536) |
| 出口 | x=5.032269, r=0.798452 (φ1.59690) | x=5.032269, r=0.823474 (φ1.64695) |
| 入口 (直管) | x=−1.368283, r=0.799919 (φ1.6) | 同左 |
| 収縮開始 | x=−1.263030 (=−6 r_t) | 同左 |
| 全長 / 拡大部 | 6.4006 / 5.0323 | 同左 |
| 最大壁角 | 20.54° @x=0.255 | ほぼ同じ |

代表点 (Euler 版, [m]):

| x | −1.3683 | −1.0525 | −0.6315 | −0.2105 | 0 | 0.1053 | 0.2105 | 0.4210 | 0.8420 | 1.6840 | 2.5261 | 3.3681 | 4.2101 | 5.0323 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| r | 0.79992 | 0.78509 | 0.56442 | 0.26188 | 0.21050 | 0.22414 | 0.26028 | 0.33793 | 0.46975 | 0.64198 | 0.73527 | 0.78068 | 0.79669 | 0.79845 |


## NS 版 (物理壁 = inviscid + δ\*, node y+~1 SST, TP) と δ\* 反復 (2026-08-18)

推奨点 R2/$L_U$6/$L_c$8 を A12/A13 の NS チェーン ([plan](../../plans/accepted/tooling-nozzle-axismach-viscous-deltastar.md),
`runner_axismach.prepare_ns`/`run_staged_ns`) で回した。設定: node 低 Re SST (`wallTreatmentSST 0`, 断熱 no-slip)、TP semi-perfect
(`thermoHrefTemp 298.15`)、陰解法 block-DPLUR、soft cfl0.5 → mid cfl1 → 本段 **cfl 1** 48000 step、メッシュ 601×97
`wall_first_frac 3.5e-5` (AR 895 PASS; forge 実測 y+: 膨張部 x>3 で ≤1.5・x>8 で <0.4、**スロート近傍 x∈[−1,0.5] は 11–14**
[y+~1 は AR ゲートと両立しない — 相関見積 y1+=1 → スロート 1.2 μm])。**TP×SST×node は本 case が初適用**で、cfl 1 で素直に回った
(coarse 中継 `run_0016` → 本計算)。IC は Euler `run_0005` → coarse → 各 v の順に cross-mesh。ドライバ `run_ns_va.py`、
δ\* 抽出 `analyze_ns_deltastar_va.py` (case/41 規約: x≥8) / `build_dstar_v3_va.py` (x≥3 抽出 + x<3 は測定端の比で相関をスケール、
`--x-hi-trust` で末端の外挿勾配をトレンド化)。

| v | run | 壁の δ\* | ‖ΔM‖∞/M_d (設計区間) | 出口軸 M | 出口コア M (r<3.4 r_t) | overshoot | δ\* 固定点 (抽出/入力) |
|---|---|---|---|---|---|---|---|
| Euler | `run_0005` | — (非粘性壁) | 0.061 % | 4.1900 | 4.1903 ±0.0005 | +0.14 % | — |
| v1 | `run_0017_va_R2_LU6_Lc8_ns_v1` | 相関 (Eckert+乱流平板, 入口履歴) | **0.524 %** | 4.2149 | 4.168 (−0.5 %) | +0.62 % | 0.68 (x=8) → 1.0 (x=22): **相関が上流で最大 30–50 % 過大** (case/41 の 10 % 過小と逆) |
| v2 | `run_0020_va_R2_LU6_Lc8_ns_v2` | CFD (x>9) + 相関 (x<6) ブレンド (case/41 規約) | 0.525 % | 4.2149 | — | +1.9 % | 設計区間は v1 と同一 (x<9 が相関のまま); ブレンドのこぶで出口 ε_θ 0.12→0.61° 悪化 → **本 case では不適** |
| v3 | `run_0021_va_R2_LU6_Lc8_ns_v3` | CFD x≥3 全域 + x<3 相関×0.51 | 0.289 % | 4.1899 | 4.183 | +0.39 % (x≈22 のこぶ = 末端外挿勾配 0.63° の artefact) | 0.998 [0.98, 1.03] |
| **v4** | **`run_0022_va_R2_LU6_Lc8_ns_v4`** | v3 + 末端 (x>21.5) をトレンド勾配 0.006 で外挿 | **0.309 %** | 4.1876 | **4.183 ±0.002 (−0.17 %)** | **+0.004 %** | **1.001 [0.995, 1.004] = 収束** |

- **結論: v4 の物理壁が NS の到達点** (δ\* は固定点、0.5 % ゲート内、出口コア M 4.183 で ±0.05 % 一様、コア内 |θ|≤0.07°)。形状は
  `points_va_best_v4_ns.csv` / `nozzle_va_best_v4_ns.geo` (`export_wall_va.py --kind ns --dstar-csv run_0021…/dstar_v4.csv`;
  スロート φ423.3 mm @x=−0.7 mm、出口 φ1649.8 mm、δ\* 出口 26 mm)。相関 δ\* 版 (`points_va_best_ns.csv`) は **使わない**
  (上流で δ\* 過大 → 出口 M −0.5 %)。
- **残差は δ\* で表現できない**: x≈1–3 の −0.2 % と x≈7 の −0.3 % (v1/v3/v4 で共通、Euler には無い) はスロート近傍の粘性/遷音速効果で、
  δ\* を変えても動かない。消すには law 側の RANS 帰還 (A5 の RANS 版、**未実装・未確立** — case/41 plan の残件) が要るので本作業では行っていない。
  スロート近傍 y+ 11–14 の解像も一因の可能性 (未切り分け)。
- 収束: `check_convergence` は全 NS run **NOT CONVERGED** (warm 床 plateau: rms_ro 4e-7, roK 3.2–3.6 桁降下)、`check_quasisteady`
  machmax/pmax **STEADY**、軸 M 16k 以降 |ΔM|≤5e-5 で凍結、NaN 0、P∈[5.2 k, 1.139 M] Pa・T∈[253, 1060] K。
- 図: `figs/ns_dstar_iteration.png` (軸 M 誤差 / 下流軸 M / 出口断面 M — Euler・v1・v3・v4)、`figs/ns_v1_dstar_extracted.png`
  (v1 の δ\* 抽出 vs 相関)、`figs/ns_v4_dstar_fixedpoint.png` (v4 の固定点)。
- 注: `run_0018` は v2 の初回投入で **別セッションの forge 再ビルド中に binary が消え本段が rc 127** で失敗 → 削除し `run_0020` で再投入。
  新 binary の dry 経路は `run_0019` (Euler 基準の再実行) で `run_0005` と軸 M 8e-6 一致を確認済み。

## 問題定義

| ファイル | R | L_U | L_c | 備考 |
| --- | --- | --- | --- | --- |
| `problem_va_R{1.5,2,3}_LU6_Lc{7,8,9,10}.yaml` | 1.5/2/3 | 6 | 7/8/9/10 | R × L_c の 12 点 |
| `problem_va_R3_LU6_Lc11.yaml` | 3 | 6 | 11 | 単調窓上限側 (R2 は 10.88 で L_c11 不可) |
| `problem_va_R2_LU{4,9}_Lc8.yaml` | 2 | 4/9 | 8 | L_U 感度 (基準 R2/L_c8) |
| `problem_va_R2_LU6_Lc8_ns_coarse.yaml` / `problem_va_R2_LU6_Lc8_ns.yaml` | 2 | 6 | 8 | NS: coarse 中継 (365×65, frac 5e-4, cfl1) / 本計算 (601×97, frac 3.5e-5, cfl1, 48000 step) |

再生成: `cd design && ./.venv-opt/bin/python -m forge_design.evaluate.runner_axismach ../case/44.vitiated_air_wt/problem_va_R2_LU6_Lc8.yaml ../case/44.vitiated_air_wt/run_00NN_va_R2_LU6_Lc8 [--prepare-only]`
(スタディ一括は `run_study_va.py`)。CFD は 12000 step ≈ 1 分/点 (RTX 3060、段階起動込み)。

**凝縮の引き継ぎ**: H₂O 2.87 wt% (モル 4.6 %)・出口 253 K/5.2 kPa なので凝縮の可能性はある。`evaluate.tp_species: split_h2o` +
`evaluate.condensation` で回せる ([plan](../../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md)、実例 case/42 run_0063–0065)。未実施。

## 計算 run 一覧

| run | 目的・設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| **`run_0001`〜`run_0015_va_R*_LU*_Lc*`** | **R/L_U/L_c スタディ (semi-perfect 設計 + node Euler TP CFD, `thermoHrefTemp` 298.15, soft→mid→本段 cfl2 12000 step; `problem_va_*.yaml`, ドライバ `run_study_va.py`)**: R∈{1.5,2,3}×L_c∈{7,8,9,10} @L_U6 (0001–0009, 0012–0014)、R3/L_c11 (0015)、L_U∈{4,9} @R2/L_c8 (0010/0011) | 上表。**全点 0.5 % ゲート内、推奨 R2/L_U6/L_c8 (`run_0005`, 0.061 %) ≒ L_c9 (`run_0013`, 0.047 %)**。L_c7 は崖の縁 (R3 0.315 %)、L_U4 不可 (0.213 %)。全 run 品質 PASS・NaN 0・軸 M 凍結 ~1e-5、`check_convergence` は warm 床 plateau で NOT CONVERGED (収束とは書かない)。`study_cfd_va.json`、`study_cfd_va.png`、`study_axisM_va.png`。run_0005/0013 の入力・設計成果物 (config, `wall_design.csv`, `metrics.json` 等) は git にリファレンス保存 | active (**スタディの正本・推奨の根拠**) |
| **`run_0016_va_R2_LU6_Lc8_ns_coarse`** / **`run_0017_…_ns_v1`** / `run_0020_…_ns_v2` / `run_0021_…_ns_v3` / **`run_0022_…_ns_v4`** | **推奨点の NS (A12/A13: 物理壁 = inviscid + δ\*, node y+~1 低 Re SST, TP)** と δ\* 反復 (`problem_va_R2_LU6_Lc8_ns{_coarse,}.yaml`, `run_ns_va.py`)。0016 = coarse 中継 (IC 供給)、0017 = v1 相関 δ\*、0020 = v2 (case/41 規約ブレンド)、0021 = v3 (CFD δ\* x≥3 全域)、0022 = **v4 (末端勾配修正、固定点)** | 上の NS 節の表。**v4 が到達点: ‖ΔM‖∞ 0.309 %、出口コア M 4.183 ±0.002、overshoot +0.004 %、δ\* 固定点 1.001**。相関 δ\* は上流で 30–50 % 過大 (v1 出口 M −0.5 %)。残差 (x≈1–3, 7 の −0.2〜−0.3 %) は δ\* 非表現 (law 側 RANS 帰還が要る、未実装)。全 run 品質 PASS・NaN 0・STEADY、`check_convergence` NOT CONVERGED (warm 床)。0022 の入力・`wall_physical.csv`・δ\* CSV は git 保存 | active (**NS 到達点の正本**) |
| `run_0019_va_R2_LU6_Lc8_newbin` | 別セッションで再ビルドされた forge binary (2026-08-18 02:52, 凝縮改修中) の dry 経路回帰確認 (`run_0005` と同条件) | 軸 M 差 8e-6、‖ΔM‖∞ 0.061 % 同一 → 新 binary で継続可 | 破棄予定 (記録) |
| `run_0018_va_R2_LU6_Lc8_ns_v2` | v2 初回 — 本段開始時に binary 再ビルド中で rc 127 | 結果なし、削除済 (`run_0020` で再投入) | 削除済 |
| — | 設計のみスイープ 21 点 → `study_design_va.json` (`study_design_va.py`; R2/L_c11 のみ単調窓外、他全合格) | | ref |
| — | 推奨形状エクスポート `export_wall_va.py` → `points_va_best_{euler,ns}.csv` / `nozzle_va_best_{euler,ns}.geo` (+ L_c9 版) | | ref |

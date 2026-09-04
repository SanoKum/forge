# ノズル設計 SERN チェーン (⑤): 燃焼器出口 starting line + 2D 平面最大推力逆設計 (key point dv) + 多作動点 RANS 評価 + MOO

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`  <!-- 2026-09-04 起票 (branch feature/sern-design)。S0–S1 実装済み (同日)、次 = S2/S3 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「SERN チェーン」節 (現在仕様。本計画と同時に起草)
  - [`design/CAPABILITIES.md`](../../design/CAPABILITIES.md) (問題タイプ `sern_2d` を 📋 で登録)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) §4.6 ⑤ (本計画で差し替え。旧「壁圧 Bézier dv + 局所帰還 + 3D FFD」は撤回)
  - 資産元: [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md) (MOC カーネル・semi-perfect ガス・CFD-in-the-loop の運用)、[`../accepted/tooling-nozzle-moo-loop.md`](../accepted/tooling-nozzle-moo-loop.md) (`opt/` の DOE→Kriging→NSGA-II→EHVI)、[`../accepted/tooling-nozzle-axismach-viscous-deltastar.md`](../accepted/tooling-nozzle-axismach-viscous-deltastar.md) (NS 場からの δ\* 抽出)
  - 出典調査: [`notes/investigations/sern-design-method-survey.md`](../../notes/investigations/sern-design-method-survey.md) (2026-09-04。方針の根拠は全てここ)
- **created**: `2026-09-04`
- **owner**: `sano`

## 1. 目的

⑤ SERN (single expansion ramp nozzle) の設計チェーンを、調査 (上記ノート §4) の推奨どおり
**「燃焼器出口 starting line → 2D 平面最大推力理論の key point を目標にした逆 MOC → forge 2D RANS を
作動点セット (設計 NPR + オフデザイン) で評価 → 既存 `opt/` の MOO」** として実装する。
完成時には、問題定義 YAML (`type: sern_2d`) から ランプ/カウル形状・多作動点の $C_T, C_L, C_M$・
パレート (設計点推力 / オフデザイン推力 / ピッチモーメント / 長さ) が一気通貫で出て、
NASA TM X-71972 の基準形でカウル長・カウル角の傾向が再現される状態を得る。
壁圧規定は主線から外し、剥離制約の判定量と二段膨張オプション (S8) に限定する。

## 2. スコープ

- **やる**
  - 問題定義 `sern_2d` (入口 = 燃焼器出口の超音速一様状態を既定、音速スロート給気は接続点として選択可)
  - 2D 平面・非対称 2 壁 (ランプ + カウル) MOC (前進) と、平面最大推力理論 (Guderley–Hantsch/Rao の平面版) の **key point** $(M_c,\theta_c,\dot m_c/\dot m)$ を目標にしたランプ壁の逆生成
  - カウル後縁以降の自由境界 (等圧せん断層) と外部流圧の整合 (設計点は無波、オフデザインは CFD に委ねる)
  - 3 ブロック構造メッシュ (内部 / カウル下外部流 / 下流プリューム) の決定的生成 (msh4.1 直書き) と品質ゲート
  - forge 2D 平面 RANS (SST, node) を作動点セットで回す runner・力積分メトリクス ($C_T,C_L,C_M$, 剥離位置)・ゲート
  - δ\* 一発補正 (forge NS 場から抽出 → 法線オフセット → 再評価)
  - 多作動点束ねの MOO (既存 `opt/` の目的関数ベクトル化)
  - 2D パレート数点の 3D RANS 確認 (側壁・隅 R・有限スパン) と、3D 設計 (流線追跡 / FFD) の要否判定
- **やらない**
  - 壁圧 $p_w(x)$ を dv にする逆設計・局所 $p\to\theta$ 帰還 (撤回。理由は調査ノート §3)
  - NS 帰還ループ (壁を反復更新する帰還エンジン)。δ\* は一発補正のみ
  - 3D FFD in-loop (S7 で乖離が大きい場合に別 plan で再考)
  - 化学非平衡 MOC (Cain §5 の弱結合法)。凍結 semi-perfect で開始、非平衡は [`chemistry-finite-rate-h2.md`](chemistry-finite-rate-h2.md) 以降
  - 非一様入口 (接触不連続) 用の回転流 MOC (S5 オプション。初版は質量平均一様状態)
  - TBCC over/under のモード遷移・可変幾何

## 3. 前提と既存資産の照合

| 部品 | 既存資産 | ギャップ・判断 |
| --- | --- | --- |
| 平面 MOC 単位過程 | `geometry/moc_kernel.py` (`delta=0` で平面、semi-perfect ガス対応、放射源流/平面単純波則で検証済み) | **再利用**。非対称 2 壁 (上壁 = ランプ、下壁 = カウル、両角部の Prandtl–Meyer 扇、カウル後縁の等圧自由境界) の march は新規 (`moc_sern.py`) |
| 最大推力理論 | なし (③ベルは TOP 幾何 dv で理論を使っていない) | **新規** `rao_planar.py`。平面版は本文未入手のため Rao 1958 / Guderley–Hantsch 1955 / Cain 式 4.1–4.3 から自前導出し §6 の解析検算で担保する |
| 逆 MOC (目標 → 壁流線) | `geometry/moc_inverse.py` (軸目標 → C⁻ 後退 → 壁流線)。①用で軸対称・軸目標 | 構造は流用するが目標が「最終特性線上の状態」に変わる。ランプ壁は kernel 終端 C⁻ と目標 C⁻ の間を前進 MOC で埋めた場の壁流線 (質量流量法) として得る |
| starting line | `feedback/cfd_anchor.py` (CFD 場から線状態抽出)、`transonic.py` (Sauer/Hall) | 既定は一様燃焼器出口 (解析)。音速スロート給気モードでは既存 Hall + kernel MOC で $M\approx1.3$–1.5 の C⁻ まで内部ノズルを作り接続 |
| ガス | `gas/` (CPG / semi-perfect NASA-9, MOC・forge とも対応) | 再利用。凍結組成 |
| メッシュ | `meshing/mesh2d.py` (軸対称 1 ブロック、msh4.1 直書き)、case/23 の plume 6 ブロック `.geo` | **新規** `mesh_sern.py`: 3 ブロック TFI (内部 / カウル下外部流 / 下流)。mesh2d の station 配置と壁クラスタリング関数を流用 |
| 評価 runner | `evaluate/runner.py` (③ベル: 2 段起動・VERDICT・warm seed)、`runner_axismach.py` (TP ガス配管) | **新規** `runner_sern.py`: 作動点ごとに run を作る。段階起動は [`procedures/divergence-and-startup.md`](../../procedures/divergence-and-startup.md) 準拠。外部流入口 + 超音速ノズル入口の 2 入口 |
| メトリクス | `metrics/extract.py::thrust_metrics` (出口面積分)、`metrics/deltastar.py::deltastar_from_run` (NS 場から質量収支 δ\*) | **新規** `metrics/sern_forces.py`: ランプ・カウル内外面の $p,\tau_w$ 積分から $C_T,C_L,C_M$ (基準点指定)、剥離位置 ($\tau_w$ 符号)。出口面積分は検算用。δ\* 抽出は再利用 |
| MOO | `opt/` (DOE / Kriging / NSGA-II / EHVI / driver / polish) | 再利用。1 評価 = N run (作動点) の束ねと目的ベクトル化を driver に追加 |
| 3D | 3D median-dual (`discretization-median-dual-3d.md`)、case/37 Netgen ブリッジ | S7 確認 run のみ。設計側の 3D 生成は押し出し + 側壁 + 隅 R (Gmsh) |

## 4. 設計方針

### 4.1 座標・幾何 (無次元: 燃焼器出口高さ $H=1$)

平面 2D、$x$ = 流れ方向、$y$ = 上向き。入口面 $x=0$、$0\le y\le1$。**ランプ = 上壁** ($y=1$ から角部で
$\theta_{r0}$ だけ上向きに折れて膨張)、**カウル = 下壁** ($y=0$ から $\theta_{c0}$、下向き正)。カウル後縁
$x=L_{\rm cowl}$ 以降は外部流とのせん断層 (等圧自由境界 $p=p_{\rm ext}$)。ランプ後縁 $x=L_{\rm ramp}$。
幾何包絡は spec: $L_{\rm ramp}^{\max}$、後胴線 (ランプが越えてはいけない $y_{\max}(x)$)、$L_{\rm cowl}^{\max}$。
モーメント基準点 $(x_{\rm ref},y_{\rm ref})$ は spec (機体 CG 相当)。符号は頭上げ正。

### 4.2 問題定義 (`type: sern_2d`)

| 区分 | 内容 |
| --- | --- |
| spec | 入口: `inflow.mode: supersonic` ($M_{\rm in}$, $p_{\rm in}$, $T_t$, 組成 — 既定) または `sonic_throat` ($P_t,T_t$, スロート諸元 → 内部対称ノズルで $M_{\rm hand}$ まで)。作動点リスト `operating_points[]` = {名前, $M_\infty$, $p_\infty$, $T_\infty$, 入口状態の上書き, 重み $w_k$}。幾何包絡、モーメント基準点、$\delta^*_{\rm in}$ (入口境界層排除厚、既定 0) |
| derived | 設計点の $p_e/p_a$、出口高さ $H_e$ (逆設計の結果)、$C_{T,\rm ideal}$ (入口状態から $p_a$ まで等エントロピー膨張の推力 = 正規化基準) |
| dv ($d=6$) | **key point** $M_c$、$\theta_c$、質量流量比 $\dot m_c/\dot m$ / ランプ初期角 $\theta_{r0}$ / カウル初期角 $\theta_{c0}$ / カウル長 $L_{\rm cowl}$。任意固定可 |
| 目的 | $f_1=-C_T^{\rm design}$、$f_2=-\sum_k w_k C_T^{(k)}$ (オフデザイン束ね)、$f_3=C_M$ (または目標値との差)、$f_4=L_{\rm ramp}$ |
| 制約 | 幾何包絡 (逆設計結果が超えたら INFEASIBLE)、最低 NPR 作動点で剥離位置 $x_{\rm sep}/L_{\rm ramp}\ge$ 許容値、$C_L$ 符号 (任意) |

### 4.3 平面最大推力理論と key point 逆設計 (S1 で導出・検証)

制御面を「ランプ後縁 $e$ から出る最終 C⁻」に取り、質量流量一定・長さ (= $e$ の位置) 固定で推力
$\int(p+\rho u^2)\,dy$ を最大化する Lagrange 問題。Cain 式 4.1–4.3 の平面版 ($y$ 依存項が落ちる) から、
制御面上で $M,\theta$ を結ぶ乗数関係と、縁 $e$ での
$\tfrac12\rho_e w_e^2\sin2\theta_e=(p_e-p_a)\cot\mu_e$ が得られる。平面流では C⁻ 上の適合関係
($\theta+\nu=$ const) と併せて **制御面上の状態が一様** になることを S1 で確認する (成立すれば
制御面は直線 Mach 線)。kernel (入口一様流 + ランプ角部扇 + カウル角部扇 + カウル壁) の中で乗数関係を
満たす点の軌跡を求め、$c$ を「$c$–$e$ 間の質量流量 = 指定比」で決める、というのが順設計。

**逆設計 (NUAA 型)**: 順設計で反復して求める $c$ の状態 $(M_c,\theta_c,\dot m_c/\dot m)$ を **dv として
先に与え**、(i) kernel を前進 MOC で作り、(ii) $c$ を kernel 内の指定質量流量比の C⁺ 上で $M=M_c$ の点として
探し、(iii) $c$ から目標状態の C⁻ を張り、(iv) kernel 終端 C⁻ と目標 C⁻ の間を前進 MOC で埋め、
(v) 質量流量法で壁流線 = ランプ壁を抽出する。$\theta_c$ と $M_c$ から縁の関係式で設計 $p_e/p_a$ が従属に
決まる (= 設計 NPR は dv の関数)。$c$ が kernel 内に存在しない dv は INFEASIBLE として optimizer に返す。

### 4.4 カウルと外部流

カウルは $\theta_{c0}$ の直線 (S1) → 必要なら短い放物線 (S6 で dv 追加可)。カウル後縁からの波は、
設計点では後縁圧 $=p_{\rm ext}$ となるよう $p_{\rm ext}$ を derived にする (無波) か、指定 $p_{\rm ext}$ に
対する等圧自由境界 (膨張扇 / 斜め衝撃は MOC の外 → INFEASIBLE 警告) として扱う。オフデザインの
波・剥離・外部流干渉は全て forge が解く。**設計に外部流を入れない** のは NASA 1974 と同じ割り切り。

### 4.5 starting line とガス

既定は一様燃焼器出口。非一様入口 (`inflow.profile:` CSV) は S5 オプションで、初版は質量平均した
一様状態を使い、非一様の影響は forge の入口分布 BC ([`boundary-inlet-profile.md`](../accepted/boundary-inlet-profile.md))
で評価側だけに入れる。ガスは semi-perfect 凍結 (`gas.model: semiperfect`, `thermo_href_temp: 298.15` 必須)。

### 4.6 δ\* 一発補正

非粘性設計壁で forge RANS (設計点) を 1 回回し、`metrics/deltastar.py::deltastar_from_run` で
ランプ・カウル両壁の $\delta^*(x)$ を質量収支定義で抽出 → 法線オフセット → 再評価。帰還は回さない。
入口境界層は spec の $\delta^*_{\rm in}$ を入口分布 BC に反映する。効果は $\Delta C_T$ として台帳に残す。

### 4.7 評価 (forge)

- メッシュ: 3 ブロック (A 内部: ランプ–カウル間 / B カウル下の外部流: 入口 $x<0$ から / C 下流: ランプ–遠方)。
  TFI、壁クラスタリング ($y^+\approx30$–80 + `wallTreatmentSST=1`、AR ≤ 1000)、`check_mesh_quality.py` VERDICT。
  node config で変換 (RANS: no-slip 壁 bcond 必須)。
- 境界: ノズル入口 = 超音速 Dirichlet (全量指定)、外部入口 = 自由流、出口 = `outlet_statPress` + 逆流 Pt/Tt、
  遠方 = 自由流。`space.pRef` = 外部静圧。
- 段階起動: 1 次 + 低 CFL → 2 次、作動点間は warm restart (同一メッシュ = `restart_field.py`)。
- ゲート: `check_convergence.py` PASS、`check_quasisteady.py --quantity` (力係数・剥離位置)。低 NPR の
  RSS/FSS は `OSCILLATING` 前提 → 平均±振幅で報告し、振幅が閾値超なら `SUSPECT`。
- メトリクス: $C_T=F_x/(p_{\rm in}A_{\rm in}\gamma M_{\rm in}^2)$ 系の無次元 (実装時に $C_{T,\rm ideal}$ 正規化を選ぶ)、
  $C_L$、$C_M$ (基準点)、剥離位置。NASA 流の帳簿: ランプ + カウル内面 + カウル外面 (外部流側) を含め、
  制御体積を明示する。

### 4.8 MOO

`opt/driver.py` を「1 評価 = 作動点数分の run」に拡張し、目的ベクトル $(f_1..f_4)$ と制約を返す。
EHVI は 2 目的閉形式のみなので、3 目的以上は NSGA-II 側の hypervolume 近似 (既存 `moo.py` の扱いに従う)。
DOE は $10d=60$ 点、infill 30–50 点。1 評価の実測時間 (2D RANS × 3 作動点) を S6 冒頭で取る。

### 4.9 3D 確認 (S7)

2D パレートから 2–3 点を押し出し + 側壁 + 隅 R (spec) で 3D メッシュ化 (Gmsh)、3D RANS (node) で
$C_T, C_L, C_M$ を 2D と比較。推力差 < 2% なら 2D 設計を正、揚力/モーメントの差は 3D 補正表として持つ。
乖離が大きければ流線追跡 (親場に横方向膨張が要る) か FFD の plan を別途起票。

## 5. 実装ステップ

| Step | 内容 | 主要ファイル | 規模 |
| --- | --- | --- | --- |
| **S0** | `sern_2d` 問題定義 (spec/derived/dv/作動点スキーマ、過拘束検査)・幾何コンテナ (ランプ/カウル折れ線、包絡検査)・CAPABILITIES 更新 | `probdef.py`, `geometry/sern_geometry.py`, `design/CAPABILITIES.md` | 小 |
| **S1** | 平面最大推力理論の導出と key point 逆設計: 非対称 2 壁前進 MOC (角部扇・カウル後縁自由境界)、乗数関係・縁条件、$c$ 探索、目標 C⁻ 張り、壁流線抽出。解析検算 (§6) | `geometry/moc_sern.py`, `geometry/rao_planar.py`, `design/tests/run_sern_moc_tests.py` | 中 |
| **S2** | 3 ブロック TFI メッシュ生成 (msh4.1 直書き)・品質ゲート・node 変換の配管 | `meshing/mesh_sern.py` | 中 |
| **S3** | 評価 runner (作動点ごとの run 生成、段階起動、warm restart)・力積分メトリクス・ゲート | `evaluate/runner_sern.py`, `metrics/sern_forces.py`, `evaluate/health.py` | 中 |
| **S4** | 検証: (a) 設計点 Euler で MOC の $C_T$ と forge の一致、(b) NASA TM X-71972 基準形 (ランプ 20°/18.5H、カウル 6°/3.12H) の M10/M4 傾向再現、(c) 外部流ブロックの妥当性 (case/23 プリューム資産と突合) | `case/46.sern_design/` (新設、README run 一覧) | 中 |
| **S5** | δ\* 一発補正の配管 (抽出 → オフセット → 再評価) と $\Delta C_T$ 台帳。非一様入口の評価側 BC | `feedback/deltastar_sern.py` (薄いラッパ), `runner_sern.py` | 小 |
| **S6** | MOO: 多作動点束ね・目的ベクトル・制約・DOE → パレート。1 評価時間の実測と予算更新 | `opt/driver.py`, `opt/moo.py`, `design/problem_sern_*.yaml` | 中 |
| **S7** | 3D 確認 run (側壁・隅 R)・2D との差の表・3D 設計要否の判定 | `case/46.sern_design/mesh3d/`, README | 中 |
| S8 (任意) | 二段膨張オプション: 基部の壁圧プラトー規定 (④延長部の壁圧帰還と共通機構) + 延長部は最大推力 | 別 plan に切り出す | 中 |

S0→S1 は CFD 不要で先行できる。S2–S3 は S1 と並行可 (固定形状で配管)。

## 6. 検証

- **単体 (S1)**: (i) 平面単純波則 ($\theta+\nu$ 保存) 1e-10、(ii) 対称極限 (カウル = ランプの鏡像) で
  平面最小長ノズルの既知解 ($\theta_{\max}=\nu(M_e)/2$、Argrow–Emanuel) を再現、(iii) $p_e/p_a=1$ で
  制御面が一様平行流、(iv) 制御面上の推力積分 = 出口面の推力積分 (保存則、<0.1%)、(v) 逆設計で
  張った目標 C⁻ を前進 MOC で再計算した状態が目標と $M$ 0.1% / $\theta$ 0.05° 以内、(vi) 格子収束。
- **CFD (S4)**: 設計点 Euler (`visc: 0.0`, `thermCond: 0.0`) で $C_T$ が MOC 値と 1% 以内、RANS で
  全作動点 `check_convergence.py` PASS (低 NPR は `check_quasisteady.py` の統計)。NASA 基準形で
  「カウル短縮 → $C_T$ 低下・頭下げ $C_M$」「カウル角 6° 近傍で $C_M$ 最良」の傾向再現。
- **メッシュ**: 各 run で `check_mesh_quality.py` VERDICT PASS を README に記録。
- **δ\* (S5)**: 補正前後で $C_T$ の変化が符号・桁で妥当 (壁摩擦込みで −1〜−3% 程度) であること。
- **3D (S7)**: 2D vs 3D の $C_T$ 差 < 2% を「2D 設計で足りる」判定基準とする。

## 7. 影響範囲

- 新規: `design/forge_design/geometry/{sern_geometry,moc_sern,rao_planar}.py`, `meshing/mesh_sern.py`,
  `evaluate/runner_sern.py`, `metrics/sern_forces.py`, `design/tests/run_sern_moc_tests.py`, `case/46.sern_design/`
- 変更: `probdef.py` (`KNOWN_TYPES` に `sern_2d`), `opt/driver.py` / `opt/moo.py` (多作動点・多目的), `design/CAPABILITIES.md`
- forge 本体: 変更なし (2D 平面 RANS・超音速入口・外部流は既存機能)。不足が出たら別 plan
- 文書: `methods/design/overview.md` 「SERN チェーン」節を実装に同期、親 plan §4.6 ⑤ / §5 Phase 4

## 8. 未確定事項 (ユーザ確認)

1. **入口の既定**: スクラムジェット (燃焼器出口 $M_{\rm in}$ 指定) を既定とした。TBCC/ラムジェットの音速スロート給気を主対象にするなら S1 の接続 (内部対称ノズル) を先に組む。
2. **作動点セット**: 設計点 1 + オフデザイン 2 (低 NPR 過膨張・高 NPR) を想定。飛行経路重み $w_k$ の与え方 (NASA 流の推力余裕か、巡航効率か) は案件で決める。
3. **モーメント基準点と目標**: $C_M$ を最小化するのか目標値に合わせるのか。
4. **隅 R・側壁**: S7 で spec として評価するだけ。設計変数にするかは S7 の結果で判断。
5. **外部流を設計段階に入れるか**: 入れない (NASA 流) を既定。

## 9. 完了条件

- [ ] `methods/design/overview.md` の SERN 節を実装と同期
- [ ] S0–S7 完了、§6 の判定を満たす (S8 は別 plan)
- [ ] `design/CAPABILITIES.md` の `sern_2d` を ✅ に更新 (検証ケース = case/46)
- [ ] `status: done` にして §10 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 10. 変更ログ

- `2026-09-04` — 初稿。調査ノート [`sern-design-method-survey.md`](../../notes/investigations/sern-design-method-survey.md) の推奨 (§4) を計画化。親 plan §4.6 ⑤ の壁圧 Bézier dv・局所帰還・3D FFD in-loop を撤回し本計画に差し替え。branch `feature/sern-design`。
- `2026-09-04` — **S0–S1 実装** (`geometry/{sern_geometry,rao_planar,moc_sern}.py`, `probdef.KNOWN_TYPES` に `sern_2d`,
  `tests/run_sern_moc_tests.py` 全 PASS)。平面 MOC は逆 (格子) 法 + 角部扇の解析閉包 + **TE 扇と自由境界反射の
  C⁺ レイ束** (格子だけでは扇が最初の station で 1 セルに潰れ伝播しない・楔域が古い値のまま残る、の 2 件を実測して
  導入)。閉包判定の罠: 「K⁻ が入口値のまま」は TE 扇の下流でも真になるので TE 先頭レイより上流に限定した。
  §6 解析検算: 一様流機械精度 / ランプ扇 M,θ 誤差 0 (解析閉包) / カウル反射後 M(ν_in+2θ_r) 一致 /
  **対称 MLN 極限 (M_in 1.5→M_e 3)**: θ_c=0.0000°, 出口高さ = 等エントロピー面積比 (0.01%), c–e 流量 = f (0.01%),
  K⁻ 一定 1e-12, 上半分総推力 = 出口運動量流束 (0.03%) / カウル付き SERN (M_in 2.5, 15°/5°, L_cowl 1, p_ext 0.05):
  格子倍で C_T 変化 0.001%・L_ramp 0.04%。dv 掃引 33 点 (f 0.35–0.6, M_c 3.2–4.4): 等長候補間で C_T 最大点が
  縁条件残差最小側 — 粗いので最適性の確証は別途 (fine sweep か等長拘束の直接最適化)。
  **速度**: kernel march 14 s (nj 301, dx 2e-3, x 9H)、逆設計 1 本 <1 s。図: scratchpad `sern/` (artifact 化)。
  **未対応**: semi-perfect ガス (圧力比写像)、非一様入口、p_ext > p_TE の衝撃。


# node 壁関数 DOF の $\omega$ 生産閉包 (2×2 分離 → 壁法則整合 $P_\omega$)

## メタ

- **area**: `turbulence / boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5 (automatic wall treatment)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §3.7
- **related_plans**:
  - [`../accepted/turbulence-node-wall-function-coverage.md`](../accepted/turbulence-node-wall-function-coverage.md) (第一内層への `wf_pk` 適用を入れた計画 = 本件の前段)
  - [`../accepted/turbulence-sst-su2-taw-coupling.md`](../accepted/turbulence-sst-su2-taw-coupling.md) §10c (`nodeOmegaWfDirichlet` を case/40 で棄却した記録)
  - **未起票**: 一般内部場の $P_\omega$ を SST-2003 正式形へ直す件 (§6。影響範囲が大きいので別 plan)
- **created**: `2026-08-12`
- **owner**: `sano`

## 0. ★ 前提の見直し (2026-08-12(4)): node 離散化は無罪、SU2 も同じ値を出す

**本 plan の出発点だった「node 壁関数に forge 固有の欠陥がある」という前提は弱まった。**
外部相関比 $C_f/C_{f,\mathrm{KS}}$ (有効帯):

| 構成 | 壁処理 | 離散化 | $C_f/C_{f,\mathrm{KS}}$ | 収支 |
| --- | --- | --- | --- | --- |
| `run_0007` forge cell | 壁解像 | cell | 0.954 | 1.069 |
| **`run_0050` forge node** | **壁解像** | **node** | **0.943** | **1.036** |
| `run_0044` forge cell | 壁関数 | cell | 0.960 | 1.035 |
| `run_0042` forge node | 壁関数 | node | 0.762 | 0.993 |
| **`run_0048` SU2 (同一メッシュ)** | **壁関数** | **node** | **0.754** | 0.775 |

1. **壁関数を外すと node は 0.943** = cell 壁解像 (0.954) と 1.1% 差、収支はむしろ node が良い
   → **median-dual 離散化そのものは健全**。欠損は**壁関数経路に限定**される。
2. ~~SU2 が同一メッシュで 0.754 を出す~~ → **撤回**。SU2 壁関数版 (`run_0048`) は
   **運動量積分を 21% 破っており** (壁 $C_f$ / 積分 $C_f$ = 0.790)、残差も rms −4.4 止まり
   (壁解像版は −12.4)。**解として棄却**し、forge node の値を裏づける証拠には使わない。
3. **代わりに、良く収束した SU2 壁解像 (`run_0049`, rms −12.4, 収支 0.999) が
   forge node 壁解像 (`run_0050`, 0.943) と 0.5% 以内で一致**した
   (SU2 0.938–0.944)。→ **node 離散化の健全性はより強い根拠で確認された**。
4. 本ケースの壁解像基準値は **$C_f/C_{f,\mathrm{KS}}\approx0.94$** (forge node 0.943 /
   SU2 0.938–0.944 / forge cell 0.954)。**node 壁関数の 0.762 はこれに対し −19%**。

**本 plan の位置づけ (2026-08-12(5) 時点)**: node 離散化は SU2 で検証済み (無罪) なので、
調査対象は **node 壁関数経路に確定**。目標値は壁解像基準の
$C_f/C_{f,\mathrm{KS}}\approx0.94$ で、現状 0.762 = **−19%**。
SU2 壁関数との比較は、SU2 側の収束/収支が改善するまで保留 (§7.5)。

**未解明**: ①本ケースの freestream $\mu_t/\mu=65.9$ は NASA 標準 (~0.009) の 4 桁大で、
壁解像 3 者が K–S より 5–6% 低い理由の候補 (未検証)。
②`run_0007` (forge cell 壁解像、これまでの内部基準) の運動量収支が 1.076/1.039 と
3 者中で最も悪く、**内部基準自身に疑いがある**。node/SU2 壁解像 (収支 1.00) の方が素性が良い。

## 1. 現象と機構 (実測で確定している範囲)

node × SST 壁関数の $\tau_w$ が系統的に低い。**外部相関 (Kármán–Schoenherr) 比**で
有効帯 $x\in[0.60,0.80]$ では **node 0.766–0.775 / cell y+30 0.959–0.963 / cell 壁解像 0.954–0.959**
(§7)。内部基準比では node 0.843 / cell y+30 1.027。**第一内層ノードの $\omega$ が壁ノードのピン値より高く跳ねる**
(壁 1.146e5 → 第一内層 **2.645e5** → 第二 9.635e4)。

### 1.1 $\omega$ スパイクの正体 = 局所源項平衡

`run_0042` 第一内層 (x=0.6) の再構成:

| 量 | 値 |
| --- | --- |
| 解像ひずみ $S$ | 1.0750e5 s⁻¹ |
| $F_1$ | 0.621 → $\alpha$=0.5118, $\beta$=0.07796 |
| 局所平衡 $\sqrt{\alpha/\beta}\,S$ | 2.754e5 |
| **実測 $\omega$** | **2.645e5** (比 0.960) |

[`ransSource_d.cu:178`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) の
$P_\omega=\alpha\rho S^2$ と $D_\omega=\beta\rho\omega^2$ の平衡でほぼ決まっている。
**$P_\omega$ は `wf_pk` と結合していない**ので、`wf_pk` を変えても $\omega$ は動かない
(実測: `nodeKwfDirichlet:0` で実効生産 3.1 倍 → $u_\tau$ +0.55%)。

入力ひずみが壁法則に対し過大:

| $y^+$ | 解像 $S$ | $S_{\rm wf}=u_\tau^2g/\nu$ | $S/S_{\rm wf}$ | $(S/S_{\rm wf})^2$ |
| --- | --- | --- | --- | --- |
| 27.1 | 1.075e5 | 5.146e4 | 2.09 | **4.36** |
| 56.9 | 3.192e4 | 1.623e4 | 1.97 | 3.87 |
| 89.6 | 1.531e4 | 9.611e3 | 1.59 | 2.54 |

### 1.2 $\omega$ が $\mu_t$ を「直接」下げているのではない (訂正)

既定経路は $\mu_t=\rho a_1k/\max(a_1\omega,\;SF_2)$
([`turbulent_viscosity_d.cu:222`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu))。
`run_0042` の実測では**第一内層で strain 側が支配**している:

| $y^+$ | $a_1\omega$ | $SF_2$ | $SF_2/(a_1\omega)$ | 支配枝 | $\mu_t$ 予測/実測 |
| --- | --- | --- | --- | --- | --- |
| 27.1 | 8.200e4 | 1.056e5 | **1.288** | strain | 1.000 |
| 56.9 | 2.987e4 | 3.192e4 | 1.069 | strain | 1.000 |
| 89.6 | 1.570e4 | 1.531e4 | 0.975 | $a_1\omega$ | 1.000 |
| 125.5 | 1.019e4 | 9.875e3 | 0.970 | $a_1\omega$ | 1.000 |

したがって**第一内層では $\omega$ は $\mu_t$ の分母に現れない**。確実に言えるのは:

1. $\omega$ 増大 → $D_k=\beta^*\rho k\omega$ 増大 → **$k$ 低下**
2. $k$ 低下 → strain-limited $\mu_t=\rho a_1k/(SF_2)$ **低下** (分子経由)
3. $y^+\gtrsim90$ では $a_1\omega$ 枝に移るので、そこでは $\omega$ が直接効く
4. `nodeOmegaWfDirichlet` フラグは $\omega$ ピンと**同時に limiter を迂回**する

### 1.3 未分離

効いているのが (i) $\omega$ 状態 (ii) strain limiter 迂回 (iii) 両方の相互作用、のどれかは**未分離**。
$C_f$ 0.843→0.957 (`run_0046`) はこの 2 因子の合成効果であり、$\omega$ 単独の効果ではない。

## 2. 作業順序 (この順を守る)

1. **本 plan の修正** (完了)
2. **`methods/turbulence/` の現仕様更新** — 現行 $P_\omega=\alpha\rho S^2$ が SST-2003 の
   正式形と異なること (§6)、node 壁関数 DOF の構成 (どの DOF に何が課されるか) を明記
3. **収支診断の追加** (§4)
4. **2×2 分離実験** (§3)
5. **最後に E3** (§5)

E3 をいきなり本命として実装する段階にはまだ無い。

## 3. 2×2 factorial 分離実験

$\omega$ ピンと limiter 迂回は**非線形に絡む**ので「E1+E2 の和」では評価できない。
2 因子 2 水準で 4 ケースを取り、相互作用を
$\Delta_{\rm int}=Y_{11}-Y_{10}-Y_{01}+Y_{00}$ で見る。

| ケース | $\omega$ pin | limiter bypass | run |
| --- | --- | --- | --- |
| $Y_{00}$ baseline | OFF | OFF | `run_0042` (既存) |
| $Y_{10}$ | **ON** | OFF | 新規 (E1) |
| $Y_{01}$ | OFF | **ON** | 新規 (E2) |
| $Y_{11}$ combined | ON | ON | `run_0046` (既存, 現行フラグ) |

- 実装: 現行 `nodeOmegaWfDirichlet` は 2 つを束ねているので、**診断用に 2 ビットへ分解**する
  (env ゲートで可。既定経路はビット不変)。
- 応答量 $Y$: **$C_f/C_{f,\mathrm{KS}}$ (主)** / $u_\tau$ / $\theta$ / 第一内層 $\omega$・$k$・$\mu_t$ / $\mu_t$@y+100。

## 4. 収支診断 (実験の前に入れる)

### 4.1 $\omega$ 方程式の項別出力

壁近傍ノードで $P_\omega$ / $D_\omega$ / 交差拡散 $CD_\omega$ / $\omega$ 輸送 (対流+拡散) を
個別に出力する。現状は $\omega$ の**結果**しか見えず、平衡がどの項で決まっているか毎回手計算になる。

### 4.2 壁摩擦の実残差確認 (本 plan の完了条件に含める)

**`twall` 出力はモデル目標値へ再スケール済みなので使わない。** 残差へ加わる直前の flux を測る:

- $\sum$ (W–I 双対面の**接線**力) と $\sum \rho u_\tau^2 A_{\rm wall}$ の一致
- W–I で生じる**余計な法線粘性力**。
  [`viscousFlux_d.cu:174-191`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) の
  AddTauWall 再スケールは、コメントは「接線成分を再スケール」だが実装は
  `tau_x *= scale` と **traction ベクトル全体**に掛けており、法線成分 $T_n$ も同じ倍率
  (本ケースで $C_{f,\rm model}/C_{f,\rm molec}\approx2.1$) で増幅される。この量を実測する。
- inlet/outlet の運動量収支との閉合

## 5. E3: 壁関数 DOF の $P_\omega$ を壁法則整合形へ

### 5.1 定式化

壁関数が効いている**第一内層ノードのみ**で、$P_\omega$ の入力ひずみを解像 $S$ でなく壁法則整合にする:

$$S_{\rm wf}=\frac{u_\tau^2}{\nu}g,\qquad
\nu_{t,\rm wf}=\nu\left(\frac1g-1\right),\qquad
P_{\omega,\rm wf}=\alpha\rho S_{\rm wf}^2\;\left(=\alpha\frac{P_{k,\rm wf}}{\nu_{t,\rm wf}}\right)$$

$g=du^+/dy^+$ は既存の `wallLaw_reichardt_duplus_dyp`。

**実装は必ず $S_{\rm wf}^2$ を直接使う**こと。$\alpha P_{k,\rm wf}/\nu_{t,\rm wf}$ の形は
低 $y^+$ で $g\to1$ となり $P_{k,\rm wf}\propto g(1-g)\to0$、$\nu_{t,\rm wf}=\nu(1/g-1)\to0$ で
**0/0 に近づく**ため数値的に不安定。

### 5.2 適用マスク (重要)

`wf_pk>=0` は**壁ノード W と第一内層 I の両方**に立つ (実測で確認済み)。E2/E3 が狙うのは **I だけ**なので、
`isNode && wall_flag[ic]==0 && wf_pk[ic]>=0` 相当を使う。**望ましくは専用の `wf_irep_flag` を
`ransWallFunction` で立てる** — これなら cell ビット不変が構造的に保証される
(cell では `wf_irep_flag` が全 0)。

### 5.3 データ受け渡し

`ransSource` は $u_\tau$ も $g$ も持っていない (引数リストに無いことを確認済み)。また `wf_pk` から
$S_{\rm wf}$ を安定に逆算することもできない (上記 0/0)。よって
**`ransWallFunction` 側で $S_{\rm wf}^2$ を計算し `wf_sprod` 専用配列へ保存**して渡す。

### 5.4 位置づけ (SST-2003 との関係)

NASA TMR (一次情報, [tmbwg.github.io/turbmodels/sst.html](https://tmbwg.github.io/turbmodels/sst.html)) は
Menter 2003 論文の $\omega$ 生産項 $\alpha\rho S^2$ を**誤記**とし、正しくは
$\alpha\dot P_k/\nu_t$ (limited production) と明記している。また production limiter
$\min(P,\,10\beta^*\rho\omega k)$ を **$k$・$\omega$ 両方**に使うとしている。

ただし **E3 は「SST-2003 そのもの」ではない**。E3 が使うのは solver 実値の $\nu_t$ ではなく
**壁法則モデル値** $\nu_{t,\rm wf}$ なので、「SST-2003 の結合構造に合わせた**壁関数閉包**」と位置づける。
一般内部場の $P_\omega$ を正式形へ直す件は影響範囲が大きいので**別 plan** とする (§6)。

## 6. スコープ外 (別 plan へ)

- **一般内部場の $P_\omega$ を $\alpha\dot P_k/\nu_t$ へ**、および production limiter の
  $\omega$ 方程式への適用。cell/node 全ケースに効くので独立に起票し回帰を取る。

## 7. 検証と判定 (2026-08-12 改訂: 外部相関ベースへ)

### 7.1 基準の階層 (呼称も統一する)

| 役割 | 対象 | 呼び方 |
| --- | --- | --- |
| **主たる絶対物理基準** | $C_f(Re_\theta)$ Kármán–Schoenherr (NASA TMR flat-plate validation と同じ) | **外部相関** |
| 補助的な絶対基準 | Schlichting $C_f(Re_x)$、普遍速度分布 | 補助基準 |
| 収支確認 | 壁面 $C_f$ と $2\,d\theta/dx$ | 収支 |
| 内部回帰基準 | forge cell 壁解像解 (`run_0007`) | **内部基準** (「真値」と呼ばない) |
| SST 実装 verification | 条件を合わせた NASA TMR の SST/FUN3D/CFL3D 解 | verification (未実施, 別ケースが要る) |

$$C_{f,\mathrm{KS}}=\left[17.08(\log_{10}Re_\theta)^2+25.11\log_{10}Re_\theta+6.012\right]^{-1}$$

### 7.2 前提ゲート (満たさない station は判定に使わない)

実測済み (§9 の 2026-08-12(2) と [case README](../../case/26.flat_plate_sst/README.md)):

- **ZPG**: 加速パラメータ $K=(\nu/U_e^2)dU_e/dx=2.1\times10^{-9}$ ($\ll3\times10^{-6}$) → **実質 ZPG で合格**。
- **発達乱流**: $\theta(x)$ が単調増加になるのは $x\gtrsim0.05$ (node)。ただし
  **BL 内 peak $\mu_t/\mu$ が自由流値 65.9 を超えるのは $x\approx0.25$–$0.37$** なので、
  それより上流は自由流乱流に埋もれる → **$x\ge0.3$ を発達条件とする**。
- **K–S 有効域**: 本ケースの $Re_\theta$ は $x=0.90$ でも **5571–6427** 止まりで
  **上限 13000 に届かない**。$Re_\theta>4000$ は $x\gtrsim0.60$ (cell) / $x\gtrsim0.75$ (node)。
- **収支の信頼域**: $C_f/(2d\theta/dx)$ は $x=0.45$–$0.75$ で 0.99–1.05 だが
  $x\ge0.9$ で崩れる (局所 fit の station 不足 + 出口影響) → **$x\le0.8$ に限る**。

→ **判定に使う帯は $x\in[0.60,\,0.80]$** ($Re_\theta\approx3900$–$5500$)。
**K–S 主比較域の下端しかカバーできない**ので、K–S を**唯一の**合否基準にはしない。

### 7.3 数値許容差 (事前に明文化。後出ししない)

現状調査で観測したばらつき:

| 項目 | 観測されたばらつき | $C_f/C_{f,\mathrm{KS}}$ への影響 |
| --- | --- | --- |
| $\theta$ 積分上端 ($y=0.2$ / $0.1$ / $2\delta_{99}$) | 0.4–3.8% | $\pm0.8\%$ ($C_{f,\mathrm{KS}}$ の $Re_\theta$ 感度は $\partial\ln C_f/\partial\ln Re_\theta\approx-0.2$) |
| 有効帯内の station 間ばらつき | cell 0.5% / node 1.2% | 同左 |
| 収支 $C_f/(2d\theta/dx)$ | $\pm5\%$ | — |

→ **$C_f/C_{f,\mathrm{KS}}$ は $\pm0.02$ の精度で議論する**。これより小さい差は有意としない。

### 7.4 合否基準

- **主判定**: 有効帯 $x\in[0.60,0.80]$ で **$C_f/C_{f,\mathrm{KS}}\ge0.94$**。
  現状: **node 0.766–0.775 / cell y+30 0.959–0.963 / cell 壁解像 0.954–0.959**。
  cell y+30 が 0.96 帯なので、node がそこに並べば合格。
  **cell への接近そのものは合格条件にしない** (cell 自身が外部相関に対し 4–5% 低い)。
- **収支**: 有効帯で $C_f/(2\,d\theta/dx)=1.00\pm0.05$。
- **y+ 非依存**: `run_0042`/`0040`/`0043` (第一 DOF y+27/52/102) で $C_f/C_{f,\mathrm{KS}}$ の
  ばらつきが $\pm0.02$ 以内 (現状 0.766–0.788 = ±0.011 で既に満たしている)。
- **case/40 ベルノズル**: $\tau_w$/y+1 基準が現行 0.945 から**過大化しない**こと
  (`nodeOmegaWfDirichlet` は 1.237 に過大化して棄却された。ここが採否の分かれ目)。
- **OFF 回帰**: 既定でビット不変 (cell 全ケース + node 非壁関数ケース)。
- **E1/E2/E3 は §3 の 2×2 factorial で交絡を分離**してから判定する。
- **収束・準定常**: 残差は本 case で構造プラトー (全 run `NOT CONVERGED`) なので、
  `check_quasisteady.py --quantity theta,cf_retheta --cf-x 0.6` の VERDICT を根拠に貼る
  (2026-08-12 に同ツールへ `theta`/`cf_retheta` を追加済み)。

### 7.5 未実施 (別途)

**NASA SST 数値解との verification** は、本ケースの freestream $\mu_t/\mu=65.9$ に対し
NASA SST-Vm 標準ケースが約 0.009 と **4 桁違う**ため apples-to-apples にならない。
流入乱流量を合わせた別ケースを作る必要がある (本 plan のスコープ外)。

## 8. 完了条件

- [x] plan の機構記述を実測に合わせて修正 (§1.2)
- [x] 検証基準を外部相関 $C_f(Re_\theta)$ ベースへ改訂 (§7)、前提ゲートと許容差を明文化
- [ ] `methods/turbulence/` 更新 (現行 $P_\omega$ と SST-2003 の差、node 壁関数 DOF 構成)
- [ ] $\omega$ 項別収支の診断出力 (§4.1)
- [ ] 壁摩擦の実残差確認 (§4.2) — W–I 接線力 vs $\rho u_\tau^2A$、余計な法線力、inlet/outlet 収支
- [ ] 2×2 分離実験と相互作用 $\Delta_{\rm int}$ の評価 (§3)
- [ ] E3 実装 (`wf_irep_flag` + `wf_sprod`) と §7.4 の判定
- [ ] `methods/index.md` 同期、`status: done` で `plans/accepted/` へ移動、[`plans/README.md`](../README.md) 同期

## 9. 変更ログ

- `2026-08-12` — 起票。case/26 `run_0040`–`run_0046` の切り分けを反映。
- `2026-08-12 (2)` — **レビュー指摘を反映して全面改訂**。
  ①機構記述を訂正: 第一内層では **strain limiter が支配** ($SF_2/(a_1\omega)=1.288$, $\mu_t$ 予測/実測=1.000)
  なので「$\omega$ が $\mu_t$ を直接下げる」は不正確。正しくは $\omega\uparrow\to D_k\uparrow\to k\downarrow\to$
  strain-limited $\mu_t\downarrow$。②E1/E2 を「和」でなく **2×2 factorial** に変更、相互作用項を評価。
  ③適用マスクを**第一内層のみ** (`wf_irep_flag`) に限定 — `wf_pk>=0` は壁ノードにも立つ。
  ④E3 は `wf_sprod` 専用配列が必要 (`ransSource` に $u_\tau$/$g$ 無し・逆算は 0/0 で不安定)。
  ⑤**SST-2003 の erratum を NASA TMR で確認済みに更新**、ただし E3 はモデル $\nu_{t,\rm wf}$ を使うので
  「SST-2003 そのもの」ではないと明記。一般内部場の修正は別 plan へ。
  ⑥**壁摩擦の実残差確認**を完了条件へ復帰 (`twall` 出力は再スケール済みなので使わない)。
  ⑦判定目標を cell の 1.027 でなく**壁解像基準 1.000** に変更。
- `2026-08-12 (3)` — **検証基準を外部相関 $C_f(Re_\theta)$ ベースへ改訂** (§7 全面差替)。
  既存 8 run を再解析した結果:
  ①**ZPG ゲート合格** ($K=2.1\times10^{-9}$)。②**発達乱流ゲート**: BL 内 peak $\mu_t/\mu$ が
  自由流値 65.9 を超えるのは $x\approx0.25$–$0.37$ (自由流乱流が強いため上流は判定不能)。
  ③**K–S 有効域が狭い**: $Re_\theta$ は $x=0.90$ でも 5571–6427 で上限 13000 に届かない。
  ④**収支は $x\le0.8$ でのみ信頼できる**。→ 判定帯 $x\in[0.60,0.80]$。
  ⑤**cell 自身が外部相関に対し 4–5% 低い** (0.954–0.963) と判明したため、
  「cell への接近」を合格条件から外し **$C_f/C_{f,\mathrm{KS}}\ge0.94$** を主判定にした。
  ⑥許容差 $\pm0.02$ を事前明文化 ($\theta$ 積分上端感度 0.4–3.8% → 比への影響 $\pm0.8\%$)。
  ⑦`check_quasisteady.py` に `theta`/`cf_retheta` を追加し、全 run の VERDICT を取得。
  ⑧NASA SST 数値解との verification は freestream $\mu_t/\mu$ が 65.9 vs 0.009 で
  apples-to-apples でないため別ケース化 (§7.5)。
- `2026-08-12 (4)` — **SU2 クロスチェックと node 壁解像で前提を見直し** (§0 追加)。
  ①`run_0050` (node planar 細メッシュ + `wallTreatmentSST:0` + MUSCL, 現行 binary, quasisteady
  ALL STEADY) が $C_f/C_{f,\mathrm{KS}}=0.943$ = cell 壁解像 0.954 と 1.1% 差 →
  **node 離散化は無罪**。②`run_0048` (SU2 v8.5.0 + `STANDARD_WALL_FUNCTION`, forge `run_0042` と
  同一メッシュ・BC・物性・流入乱流) が **0.754** で forge node 0.762 と 1% 以内一致 →
  **forge node 壁関数に固有の欠陥があるという前提は成立しにくい**。
  ③`run_0047` (SU2 low-Re を y+27 に適用) は収支 0.181 で破綻 = low-Re の誤用は SU2 でも壊れる。
  → plan の目的を「forge のバグ修正」から「vertex-centered 壁関数の機構特定」へ再定義。

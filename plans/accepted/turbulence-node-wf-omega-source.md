# node SST 壁関数経路の低摩擦解 (実 W–I 力診断 → 2×2 分離 → 壁法則整合 $P_\omega$)

<!-- 現在の判断は §0。§1 以降は根拠・手段。過去の撤回は §9 変更ログにのみ残す。 -->

## メタ

- **area**: `turbulence / boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5 (automatic wall treatment)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §3.7
- **related_plans**:
  - [`turbulence-node-wf-representative-point.md`](turbulence-node-wf-representative-point.md) (**後継**: $\omega$ 介入 3 種でゲートを満たせず代表点側へ)
  - [`../accepted/turbulence-node-wall-function-coverage.md`](../accepted/turbulence-node-wall-function-coverage.md) (第一内層への `wf_pk` 適用を入れた計画 = 本件の前段)
  - [`../accepted/turbulence-sst-su2-taw-coupling.md`](../accepted/turbulence-sst-su2-taw-coupling.md) §10c (`nodeOmegaWfDirichlet` を case/40 で棄却した記録)
  - **未起票**: 一般内部場の $P_\omega$ を SST-2003 正式形へ直す件 (§6。影響範囲が大きいので別 plan)
- **created**: `2026-08-12`
- **owner**: `sano`

## 0. 現状認識 (2026-08-12 確定版)

> **case/26 では forge node 壁解像が SU2 壁解像と 0.5% 程度で一致するため、node 離散化一般は
> 主因ではない。forge node 壁関数では、課した低い壁応力と境界層の遅い発達が互いに整合しており、
> 欠損は壁関数経路に局在する。ただし、狙った壁応力が W–I 残差へ正しく作用しているか、また
> $\omega$ 状態・limiter・源項のどれが低摩擦解を作るかは未分離である。SU2 壁関数 run は
> 全残差未収束のため比較証拠に使用しない。**

### 0.1 確定してよいこと

| # | 内容 |
| --- | --- |
| 1 | forge node 壁解像 (`run_0050`) と **SU2 壁解像 (`run_0049`, 全残差 −5〜−12.6)** が $C_f/C_{f,\mathrm{KS}}\approx0.94$ で一致 (4 定義とも 1% 以内) |
| 2 | → **−19% 欠損は median-dual/node 離散化一般ではなく、forge の node 壁関数経路に局在** |
| 3 | forge node 壁関数は**課した壁応力 (0.762) と BL 発達 (0.759) が互いに整合**している (収支 1.003)。壁応力の出力だけの問題ではなく**解そのものが低摩擦** |
| 4 | SU2 の `Skin_Friction_Coefficient` は壁関数応力を反映 (`CFVMFlowSolverBase.inl` で `Tau_Wall/WallShearStress` 再スケール) |

### 0.2 使ってはいけないこと

- **`run_0048` (SU2 壁関数) を定量比較に使わない**: 全残差で未収束
  (`rms[RhoE]`=**+1.03**, `rms[w]`=−0.29, `rms[RhoV]`=−2.13。`rms[Rho]` −4.43 だけ見ていた)。
  未収束場では定常運動量積分式に残差項が欠けるので、(c)=1.003 は「場が正しい」根拠にならない。
- **SU2 に Reichardt 逆解きを当てない**: SU2 は **Nichols–Nelson (2004)** の圧縮性
  Spalding/White–Christoph 型 (`CNSSolver.cpp`)。整合確認には SU2 自身の `Y_Plus`
  (壁ノードに出力、x=0.6 で 27.21)・$U_\tau$・Spalding 則を使う。
- **運動量積分 $2(d\theta/dx+\mathrm{PG})$ を「定義非依存の $C_f$」と呼ばない**:
  定常性・境界層近似・$U_e,\rho_e$・積分上端・微分 fit に依存する**収支診断**である。

### 0.3 未分離のまま残っていること

1. ~~狙った壁応力が W–I 残差へ本当に作用しているか~~ ← §4.2、**実施済 (2026-08-12)。配線は健全**
2. $\omega$ 状態 / strain limiter / 源項のどれが低摩擦解を作るか ← §3 の 2×2

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
3. **実 W–I 接線力の診断** (§4.2) ← **最優先。E1/E2/E3 より先**
4. $\omega$ 項別収支の診断 (§4.1)
5. **2×2 分離実験** (§3)
6. **最後に E3** (§5)

E3 をいきなり本命として実装する段階にはまだ無い。
**§4.2 (狙った摩擦力が W–I 間で本当に作用しているか) は 2026-08-12 に実施済み — 配線は健全。**
**残る最優先は §4.1 の $\omega$ 項別収支** (E3 の設計根拠にする)。

## 3. 2×2 factorial 分離実験 — **実施済 (2026-08-12): $\omega$ pin が最大因子**

`wf_irep_flag` (第一内層マスク、cell 常時 0、node SST 壁関数時のみ受け渡し) と env
`FORGE_WF_LIMITER_BYPASS` (−1=従来/0=OFF/1=ON) で 2 因子を分離。
判定は §7.4 の 3 点セット。**全 run 共通で $Re_\theta\ge4000$ となる帯 x=[0.640, 0.784]**
(共通下限は Y00/Y01 が律速で 0.640)、同一 IC・20000 step、派生量 `STEADY`
(残差は全 run `NOT CONVERGED` = 本 case の構造プラトー)。

| run | (1) **W–I 実力** $C_f$/KS | (2) 運動量積分 | $\Sigma\lvert F_n\rvert/\Sigma\lvert F_t\rvert$ | |
| --- | --- | --- | --- | --- |
| $Y_{00}$ pin0 byp0 (`run_0053`) | **0.7615** | 0.784 | 1.97e-5 | — |
| $Y_{10}$ pin1 byp0 (`run_0054`) | **0.8594** | 0.840 | 1.96e-5 | — |
| $Y_{01}$ pin0 byp1 (`run_0055`) | **0.7670** | 0.794 | 1.92e-5 | — |
| $Y_{11}$ pin1 byp1 (`run_0056`) | **0.8864** | 0.846 | 1.87e-5 | — |

| 効果 (応答量 = (1)) | 値 |
| --- | --- |
| 単純効果 pin (byp=0 / byp=1) | +0.0979 / +0.1194 |
| 単純効果 bypass (pin=0 / pin=1) | +0.0055 / **+0.0270** |
| **平均主効果 pin** | **+0.1087** |
| **平均主効果 bypass** | **+0.0163** |
| **相互作用** | **+0.0215** |

1. **$\omega$ pin が最大因子**。ただし**「全効果の 79%」という寄与率は撤回** —
   相互作用の配分方法に依存し一意に決まらない。
2. **bypass は「ほぼ無効」ではない**: pin=0 で +0.0054 だが **pin=1 で +0.0267** と
   事前の有意差 $\pm0.02$ を超える。→ 「単独効果は小さいが $\omega$ pin 時に正の相互作用を持つ。
   それでも pin 効果より小さい」。
3. **回復率**: $Y_{11}$=0.8864 は目標帯 0.94±0.02 に対し
   **中心 0.94 まで 70% / 下限 0.92 まで 79%** (基準を必ず明記する)。
4. **「残り ~29% は $\omega$ でも limiter でもない別要因」は撤回**。言えるのは
   **「今回試した特定の $\omega$ Dirichlet 値と bypass の組合せでは目標を再現しない」**まで。
   残差分は ①pin 値 ②pin する DOF/範囲 ③$P_\omega$/拡散/交差拡散など別の $\omega$ 機構
   ④$k$–$\omega$ 結合、のいずれもありうる。**E3 自体が $\omega$ 方程式への介入**なので
   「$\omega$ でない」と断定すると自己矛盾する。
5. **次は §4.1 の $\omega$ 項別収支を出し、それを根拠に E3 を設計する**。E3 へ直行しない。

### 3b (設計) 2×2 factorial の枠組み

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

### 4.1 $\omega$ 方程式の項別出力 — **実施済 (2026-08-13): 局所 P–D 平衡が 96%**

診断 `omg_prod`/`omg_dest`/`omg_cross`/`omg_trans` (`FORGE_OMEGA_BUDGET=1` ゲート) を追加。
`run_0059_omgbudget`、x=0.694 の第一内層ノード (y+27.2):

| 項 | 値 | 比 |
| --- | --- | --- |
| 生産 $P_\omega$ | +1.793e4 | **50.0%** |
| 消滅 $-D_\omega$ | −1.652e4 | **46.1%** |
| 交差拡散 | +2.36 | 0.0% |
| 輸送 | −1.403e3 | 3.9% |
| 残差和 | 3.08e-1 | **0.00%** |

→ $D_\omega/P_\omega$=**0.922**、$|T_\omega|/P_\omega$=**0.0783**。
**局所源項が支配的で、消しきれない生産の 7.8% を輸送が釣り合わせている**
(「P–D だけで 96% 平衡」は不正確)。**収支が閉じるのは `omg_trans` の定義上ほぼ代数的**で、
符号・帳尻の確認にはなるが独立した物理検証ではない。

**E3 の事前予測**: $\omega\propto\sqrt{P_\omega}$ かつ $P_\omega=\alpha\rho S^2$ なので、
$S\to S_{\rm wf}$ で $\omega$ は $S_{\rm wf}/S\approx0.48$ 倍 = **1.27e5** (Menter ピン値 1.14e5 相当)。
**ただしこれは凍結場の初期見積り**。E3 では生産が 0.234 倍になるので、現在の輸送量が
そのままなら新生産に対する輸送比は 7.8%→**約 33%** になり、$F_1$・$\alpha/\beta$・$u_\tau$・
$S_{\rm wf}$・$k,\mu_t$・速度場もすべて再平衡する。
→ **「$Y_{10}$ に近い / 目標未達」は事前仮説に留め確定的には書かない**。
**「残り 0.86→0.94 は $\omega$ 以外」もまだ言えない** ($\omega$ の輸送・ブレンド・適用範囲・
$k$ との結合が候補に残る)。

**E3 の位置づけ**: 目標達成を狙うのではなく
**「free な $\omega$ が源項整合でどこへ再平衡するか」を検証する**実験。

**限界**: `run_0059` は全残差 `NOT CONVERGED`・snapshot 2 個で `TRANSIENT-UNSETTLED`。
着目点の $\omega$ は 500 step で 0.0048% しか動かないので**瞬時収支サンプルとしては十分**だが、
定常収支の主張には使わない。**case/40 (軸対称) で使うなら `omg_axisym` の別出力が要る**
(`omg_trans` は $1/r$ ソース加算前)。

### 4.1b (設計) $\omega$ 項別出力の仕様

壁近傍ノードで $P_\omega$ / $D_\omega$ / 交差拡散 $CD_\omega$ / $\omega$ 輸送 (対流+拡散) を
個別に出力する。現状は $\omega$ の**結果**しか見えず、平衡がどの項で決まっているか毎回手計算になる。

### 4.2 壁摩擦の実残差確認 — **実施済 (2026-08-12): 伝達は健全**

診断 `wi_ftan`/`wi_fnrm`/`wi_ftan_res` を実装 (AddTauWall 再スケール直後で壁ノードへ `atomicAdd`、
毎ステップ 0 クリア)。`run_0051_node_wf_widiag` の結果:

| 帯 | $\Sigma$`wi_ftan` [N] | $\int\rho u_\tau^2 dx$ [N] | 比 | 法線/接線 |
| --- | --- | --- | --- | --- |
| x=[0.30,0.90] | 3.8796 | 3.7980 | **1.0215** | 0.003% |
| x=[0.20,0.95] | 4.9573 | 4.8780 | **1.0163** | 0.003% |
| x=[0.45,0.75] | 1.9399 | 1.8595 | **1.0432** | 0.002% |

1. **【表現訂正】** 比 1.016–1.043 は「伝達精度」ではない。再スケールは
   **局所的に $|F_t|=\tau_w A_{\rm WI}$ をコード自身が強制**しているので、比の差は
   双対面積の和と壁側 $\int\tau_w dx$ の幾何・積分差である。正しい主張は
   **「意図した $\tau_w A_{\rm WI}$ がそのまま内部ノード残差へ加算される → 低摩擦は
   『壁応力が届かない』では説明できない」**。壁ノード側の運動量残差は
   `nodeWallDirichlet_d.cu` でゼロ化され、内部側に反力が残る構造もコード上確認済み。
2. **法線力は発達域で無害** (絶対値 `wi_fnrm_abs` で再測定):
   **x=[0.3,0.9] で $\Sigma|F_n|/\Sigma|F_t|$=2.75e-5、局所最大 6.13e-5**。
   ただし**壁全体では 6.87e-4、前縁 $x=0$ の壁∩sym コーナーで局所 0.343** なので、
   **「本ケース全体で無害」と一般化しない** (判定帯への影響は小さい)。
3. **再スケールの向きを訂正**: 係数 **0.640** = 増幅ではなく**低減**。W–I 双対面の解像 traction は
   $\mu_{\rm total}$ と双対面距離で決まるのでモデル $\tau_w$ より大きく、AddTauWall はそれを下げている。
   「$C_{f,\rm model}/C_{f,\rm molec}\approx2.1$ 倍に増幅される」という推定は誤りだった。

**残る問い**: 伝達も法線力も無罪なので、**なぜ解が低摩擦で釣り合うのか**。→ §3 の 2×2 へ。

### 4.2b (旧記述) 壁摩擦の実残差確認の設計

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
- **収支の信頼域**: $C_f/(2(d\theta/dx+\mathrm{PG}))$ は $x=0.45$–$0.75$ で概ね 1.00±0.05 だが
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

**主判定は次の 3 点セット** (単一の指標で決めない):

1. **実際に運動量残差へ入った W–I 接線力** から求めた $C_f$
2. それと**運動量積分** $2(d\theta/dx+\mathrm{PG})$ の**収支が閉じる** (比 $1.00\pm0.05$)
3. 両者と**外部相関 K–S** の比較

**目標値は検証済み壁解像帯**: forge node 壁解像 0.936–0.943 と SU2 壁解像 (全残差収束) 0.938–0.947 から

$$C_f/C_{f,\mathrm{KS}} = \mathbf{0.94\pm0.02}$$

**一方向の $\ge0.94$ ではない** (過大側も外れとする)。現状 node 壁関数は 0.762 で **−19%**。
**cell (0.960) への接近は合格条件にしない**。

- **y+ 非依存**: `run_0042`/`0040`/`0043` (第一 DOF y+27/52/102) でばらつき $\pm0.02$ 以内。
- **case/40 ベルノズル**: $\tau_w$/y+1 基準が現行 0.945 から過大化しないこと。
- **OFF 回帰**: 既定でビット不変。
- **E1/E2/E3 は §3 の 2×2 factorial で交絡を分離**してから判定。
- **収束・準定常**: `check_convergence.py` を**全残差列**で確認 (`rms_ro` だけで判断しない) +
  `check_quasisteady.py --quantity theta,cf_retheta`。
- **数値は [`tools/cf_retheta_analysis.py`](../../case/26.flat_plate_sst/tools/cf_retheta_analysis.py)
  で再生成できること** (scratchpad のスクリプトで出した数値を文書に書かない)。

### 7.5 未実施 (別途)

**NASA SST 数値解との verification** は、本ケースの freestream $\mu_t/\mu=65.9$ に対し
NASA SST-Vm 標準ケースが約 0.009 と **4 桁違う**ため apples-to-apples にならない。
流入乱流量を合わせた別ケースを作る必要がある (本 plan のスコープ外)。

## 8. 完了条件

- [x] plan の機構記述を実測に合わせて修正 (§1.2)
- [x] 検証基準を外部相関 $C_f(Re_\theta)$ ベースへ改訂 (§7)、前提ゲートと許容差を明文化
- [x] node 離散化の無罪確認 (`run_0050` vs SU2 `run_0049`, §0.1)
- [x] 後処理の正式ツール化 (`case/26.flat_plate_sst/tools/cf_retheta_analysis.py`)
- [x] **実 W–I 接線力の診断 (§4.2)**: 意図した $\tau_w A_{\rm WI}$ がそのまま内部ノード残差へ加算される
      (伝達は健全)。法線力は**発達域**で接線の 2.75e-5 (壁全体 6.87e-4、前縁コーナー局所 0.343)。
      残: inlet/outlet 収支との閉合
- [x] `methods/turbulence/theory.md` §6.5(e) に分離結果・伝達診断・診断の使い方を反映
- [x] **$\omega$ 項別収支の診断出力 (§4.1)**: 第一内層は $D/P$=0.922・$|T|/P$=0.078 で局所源項が支配。
      収支が閉じるのは `omg_trans` の定義上ほぼ代数的 (物理検証ではない)。
      E3 の事前見積り ($\omega$→1.27e5、凍結場) を立て、実測 1.345e5 で 6% 以内に的中
- [x] 2×2 分離実験と相互作用の評価 (§3): 平均主効果 pin **+0.1087** / bypass **+0.0163** /
      相互作用 **+0.0215**。$Y_{11}$=**0.8864** で目標中心まで **70%** (下限まで **79%**) の回復
- [x] **E3 実装 (`wf_irep_flag` + `wf_sprod`) と判定**: case/26 で 0.7615→**0.8445**
      ($\omega$ 2.62e5→1.35e5、事前見積り 0.484 倍に対し実測 0.513 倍)。
      **case/40 で $\tau_w$ が 1.1226 倍 = y+1 真値比 ≈1.061 に過大化** →
      **このままでは採用不可** ($\omega$ ピンの 1.237 より遥かに軽い)。
      両ケースが同方向に同程度動くので**残りは単一ノブで閉じない**
- [ ] `methods/turbulence/` の §6.5 に仕様追記、`methods/index.md` 同期
- [ ] `status: done` で `plans/accepted/` へ移動、[`plans/README.md`](../README.md) 同期

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
- `2026-08-12 (5)` — **レビュー反映で §0 を確定版へ全面整理、§7.4 を 3 点セット + $0.94\pm0.02$ に単一化**。
  ①`run_0048` (SU2 壁関数) の**全残差未収束**を確認 (`rms[RhoE]`=+1.03) → 比較証拠から除外。
  「SU2 の場は正しい」を撤回。②SU2 は **Nichols–Nelson** 壁法則なので Reichardt 逆解きを当てた
  比較は不成立 → 撤回。③SU2 の Cf 出力は `Tau_Wall` 再スケール済み = 壁関数応力を反映 (ソース確認)。
  ④運動量積分を「定義非依存の $C_f$」と呼ぶのをやめ**収支診断**に格下げ。
  ⑤**実 W–I 接線力の診断 (§4.2) を E1/E2/E3 より前に置いた** (当初からの本題で未完了)。
  ⑥後処理を `case/26.flat_plate_sst/tools/cf_retheta_analysis.py` に正式化
  (fit 窓幅/次数・積分上端・edge 基準・壁法則をパラメータ化)。
- `2026-08-12 (6)` — **§4.2 実施: 伝達は健全と判明**。`wi_ftan`/`wi_fnrm`/`wi_ftan_res` 診断を
  `viscousFlux_d.cu` に追加し `run_0051` で実測。$\Sigma$`wi_ftan`/$\int\rho u_\tau^2dx$ = **1.016–1.043**
  → 狙った壁応力は届いている。法線力は接線の 0.003% で無害。再スケール係数は **0.640** で
  「増幅」ではなく「低減」(W–I 解像 traction は $\mu_{\rm total}$ ベースでモデル $\tau_w$ より大)。
  → 伝達・法線力の 2 仮説を消去。次は §3 の 2×2 factorial。
- `2026-08-12 (7)` — **§3 の 2×2 factorial 実施**。`wf_irep_flag` (第一内層マスク) と env
  `FORGE_WF_LIMITER_BYPASS` で $\omega$ ピンと limiter 迂回を分離。
  **主効果 pin +0.100 / bypass +0.006 / 相互作用 +0.021** → **$\omega$ 状態が支配因子**。
  ただし**両方入れても 0.889** で目標帯 0.94±0.02 に届かず、**残り ~29% は別要因**。
  既定経路 (bypass=−1) は従来と同一、cell は `wf_irep_flag`≡0 でビット不変。
  次は §5 の E3 ($P_\omega$ を壁法則整合形 $S_{\rm wf}$ へ)。
- `2026-08-12 (8)` — **レビュー反映で §3/§4.2 の主張を訂正**。
  ①W–I 診断の比 1.016–1.043 は「伝達精度」ではない (コードが局所的に $|F_t|=\tau_w A_{\rm WI}$ を強制)。
  正しい主張は「意図した力がそのまま内部ノード残差へ加算される → 低摩擦は伝達では説明できない」。
  ②法線力を**絶対値** (`wi_fnrm_abs` 追加) と局所最大で再測定: 0.0027% / max 6.1e-5 → 無害で確定。
  ③再スケール係数の**分布**を確認 (min 0.6400/median 0.6405/max 0.6410 = ほぼ一様)。
  ④factorial の用語訂正 (単純効果 vs 平均主効果)、**寄与率 79% を撤回**、
  bypass は pin=1 で +0.0267 = 有意。⑤**「残り 29% は別要因」を撤回**。
  ⑥判定を帯 x=[0.60,0.80] の**3 点セット**へ (x=0.60 単点は $Re_\theta$=3902 で K–S 域外)。
  ⑦実装: `FORGE_WF_LIMITER_BYPASS` を node SST 壁関数時のみ有効化、
  W–I 診断の memset を `FORGE_WI_FORCE_DIAG` でゲート (通常経路の性能不変)。
- `2026-08-13` — **レビュー残 4 点を解消**。
  ①法線力の結論を**発達域に限定** (x=[0.3,0.9] で 2.75e-5 だが壁全体 6.87e-4・前縁コーナー局所 0.343)。
  ②**$Re_\theta$ ゲートを全 run 共通化** — 帯を x=[0.640, 0.784] に変更 (共通下限は Y00/Y01 律速)。
  ③**W–I 帯平均を正式ツールで再生成可能に** (`--wi-force/--band/--retheta-min`)。
  規約も明示: `wi_ftan` を双対幅で割って局所 $\tau_{w,\rm applied}$ にし、$C_{f,\rm KS}q$ を
  同じ壁ノード列へ内挿して**同一の台形則**で積分 (点和 vs 台形積分の端点ずれで 10% 級の偏りが出るため)。
  数値も更新: 平均主効果 pin +0.1087 / bypass +0.0163 / 相互作用 +0.0215、回復率 70%/79%。
  ④**診断 OFF 時は `wi_*` を確保も出力もしない** (`variables.cpp` で除外)、確保直後に 0 初期化。
  ⑤`methods/turbulence/theory.md` を先に更新 (開発フロー順守)。文書の残存記述も修正。
- `2026-08-13 (2)` — **§4.1 実施**。`omg_prod`/`omg_dest`/`omg_cross`/`omg_trans`
  (`FORGE_OMEGA_BUDGET=1` ゲート、OFF なら確保も出力もしない) を `ransSource_d.cu` に追加。
  第一内層 (y+27.2) は **生産 50.0% / 消滅 46.1% / 交差拡散 0.0% / 輸送 3.9%、残差和 0.00%**
  = **局所 P–D 平衡が支配**。これを根拠に **E3 の事前予測**を立てた:
  $S\to S_{\rm wf}$ で $\omega$ は 0.48 倍 = 1.27e5 (Menter ピン 1.14e5 相当) →
  **E3 単独は $Y_{10}$=0.8594 相当で目標未達の見込み**。E3 の価値は
  「Dirichlet でなく源項整合で同じ $\omega$ を得る」点にあると位置づけ直した。
- `2026-08-13 (3)` — **E3 実装・判定**。`wf_sprod` ($S_{\rm wf}^2$) を `ransWallFunction` で作り
  `ransSource` の $\omega$ 生産のみに使う (env `FORGE_WF_OMEGA_SOURCE`、cell ビット不変)。
  軸対称収支用に `omg_axisym` も追加。
  **case/26**: $C_f$/KS 0.7615→**0.8445** ($Y_{10}$=0.8594 の近く、目標 0.94±0.02 に未達)。
  $\omega$ は 2.621e5→1.345e5 で**凍結場見積り 0.484 倍に対し実測 0.513 倍** (6% 以内)。
  ただし**輸送も再平衡** (1403→491) し「輸送比 33%」の素朴外挿は起きず 7.8%→10.5%。
  **case/40**: $\tau_w$ が生産比 **1.1226 倍 = y+1 真値比 ≈1.061** に過大化 →
  **採用不可**。ただし $\omega$ ピン (1.237) より遥かに軽く、Dirichlet を使わず
  $k$/$\omega$ の非対称を解消する点は正しい方向。
  **両ケースが同方向に同程度動く** = E3 は一貫した物理変更であり、
  case/26 不足 / case/40 過剰という**符号の違い**が残る → **単一ノブでは閉じない**。
  次は $u_\tau$ を決める代表点速度側。
- `2026-08-13 (4)` — **レビュー指摘 3 点を修正し、後継 plan を起票**。
  ①**E3 env ゲートの重大バグを修正**: `FORGE_WF_OMEGA_SOURCE` が構成を確認せず有効化していたため、
  `wallTreatmentSST!=1` (壁関数初期化が走らない) の構成で誤って付けると `wf_sprod` が確保時の 0 のまま
  `>=0` 判定に通り、**$\omega$ 生産を全域で消す**危険があった。node+SST+`wallTreatmentSST==1` の
  ときだけ有効にし、それ以外は強制 0 に。壁解像 run で実地検証済み (EXPERIMENT 非表示・$\omega$ 健全)。
  `wf_sprod` も E3 OFF 時は確保・出力から除外。
  ②**`check_quasisteady.py` を修正**: 未知 quantity を黙って無視して `ALL STEADY` を返すバグを
  **エラー**に変更し、**`wall_tau` を正式対応**。case/40 の定常性を取り直し
  (`run_0066`/`run_0071` とも `wall_tau: STEADY`、tail mean 1224/1374、drift 0.0%)。
  当初の「`ALL STEADY`」は既定量 (shock/machmax/pmax) の判定で $\tau_w$ を見ていなかった。
  ③**解釈を弱めた**: 「一貫した物理変更」→「ケース固有係数の調整ではなく両ケースに
  約 11–12% の同方向応答を与える構造的介入」。case/40 の絶対誤差は −5.5%→+6.1% で**ほぼ対称**であり
  「大幅に悪化」ではなく**事前ゲート不合格**。y+1 基準との同一定義比較の正式ツール化は後継 plan へ。
  ④case/40 run 表に `run_0071_E3` を追加、`methods/turbulence/theory.md` に E3 の仕様と不採用を反映。
- `2026-08-13 (5)` — レビュー反映: ①`check_quasisteady.py` の**新規回帰を修正** —
  未知量チェックを run 依存の extractor でなく**固定集合**に対して行うようにし、
  非対称メッシュの `asym` が unknown 扱いされる問題を解消 (適用不能な量は従来どおり skip)。
  ②`wall_tau` を **`wall_model_tau`** へ改名し汎用化 — `--wall-phys-id` / `--wall-xmin` /
  `--wall-xmax` を追加 (ケース固有値のハードコード撤廃)、**同一 step の複数壁面を点数重みで集約**、
  `with h5py.File(...)` でクローズ、**utau 全ゼロなら NOT APPLICABLE** (壁関数 run 専用と明記)。
  ③`methods/turbulence/theory.md` の case/26 誤差を **−23.9%→−15.6%** に訂正
  (誤って −15.6%→−6.1% と書いていた)。④完了条件の古い数値を本文と同期。

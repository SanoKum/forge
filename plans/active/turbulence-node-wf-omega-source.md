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

## 1. 現象と機構 (実測で確定している範囲)

node × SST 壁関数の $\tau_w$ が系統的に低い (case/26 平板 x=0.6 で $C_f$/壁解像基準 =
**node 0.843 / cell 1.027**)。**第一内層ノードの $\omega$ が壁ノードのピン値より高く跳ねる**
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
- 応答量 $Y$: $u_\tau$ / $C_f$/壁解像基準 / 第一内層 $\omega$・$k$・$\mu_t$ / $\mu_t$@y+100。

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

## 7. 検証と判定

- **目標は壁解像基準の 1.000** に置く。**cell の 1.027 に合わせない** (合わせると再び較正になる)。
- **case/26 平板**: $C_f$/壁解像基準が 0.843 から 1.000 へ近づくこと。
  y+ 掃引 (`run_0042`/`0040`/`0043` = y+27/52/102) で **y+ 非依存が保たれる**こと。
- **case/40 ベルノズル**: $\tau_w$/y+1 基準が現行 0.945 から**過大化しない**こと
  (`nodeOmegaWfDirichlet` は 1.237 に過大化して棄却された。ここが採否の分かれ目)。
- **OFF 回帰**: 既定でビット不変 (cell 全ケース + node 非壁関数ケース)。
- **収束判定**: 本 case は block-DPLUR の構造プラトーなので、残差 verdict に加えて
  **派生量 $u_\tau$ の準定常**をスナップショット間ドリフトで確認する (現行 6 run は 0.001–0.075%)。

## 8. 完了条件

- [x] plan の機構記述を実測に合わせて修正 (§1.2)
- [ ] `methods/turbulence/` 更新 (現行 $P_\omega$ と SST-2003 の差、node 壁関数 DOF 構成)
- [ ] $\omega$ 項別収支の診断出力 (§4.1)
- [ ] 壁摩擦の実残差確認 (§4.2) — W–I 接線力 vs $\rho u_\tau^2A$、余計な法線力、inlet/outlet 収支
- [ ] 2×2 分離実験と相互作用 $\Delta_{\rm int}$ の評価 (§3)
- [ ] E3 実装 (`wf_irep_flag` + `wf_sprod`) と §7 の判定
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

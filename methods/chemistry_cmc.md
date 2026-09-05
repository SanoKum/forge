# 乱流‐化学相互作用: 1st-order CMC (Conditional Moment Closure)

有限速度化学 ([chemistry.md](chemistry.md)) に、セル平均でなく**混合分率に条件付けた組成分布**で反応を評価する
TCI クロージャを加える。設計判断と進捗は [`plans/active/chemistry-cmc-tci.md`](../plans/active/chemistry-cmc-tci.md)、
文献根拠は [`notes/investigations/cabra-liftoff-model-fidelity-survey.md`](../notes/investigations/cabra-liftoff-model-fidelity-survey.md)。

**状態 (2026-09-05)**: 設計確定・実装着手 (Phase A 混合分率インフラから)。

---

## 理論

### 1. なぜ平均場評価では足りないか

RANS の平均反応率は状態 $\psi=(Y_1..Y_{n_s},T)$ の確率密度 $f(\psi)$ による積分
$\bar{\dot\omega}=\int\dot\omega(\psi)f(\psi)d\psi$ であり、Arrhenius の強い非線形性のため
$\dot\omega(\tilde\psi)$ (laminar chemistry) とは一致しない。自着火で安定化する火炎 (Cabra, BK, 加熱器) では
着火は**最反応性混合分率 $\xi_{MR}$ (超希薄) かつ低スカラー散逸のポケット**で始まり、着火遅れは分布の裾で決まる。
平均場評価は可着火な平均混合気をリップ直近から即着火させ、火炎が付着する (case/48 で実証、機構交換・PaSR では直らない)。

### 2. 混合分率と条件付き平均

2 流 (燃料流 $\xi=1$ / 酸化剤流 $\xi=0$) の非予混合場では、組成・温度のばらつきの主因は「どれだけ混ざったか」であり、
同じ混合分率の流体塊同士の組成はほぼ等しい。そこで未知数を条件付き平均に取り替える:

$$Q_\alpha(\eta;\mathbf x)=\langle Y_\alpha\,|\,\xi=\eta\rangle,\qquad Q_h(\eta;\mathbf x)=\langle h\,|\,\xi=\eta\rangle$$

$h$ は forge の sensible datum の比エンタルピー (反応熱は $\dot Q=-\sum c_s\dot\omega_s$ として入る、chemistry.md §1)。
無条件平均は混合分率 PDF $\tilde P(\eta)$ で戻す:

$$\tilde Y_\alpha=\int_0^1 Q_\alpha(\eta)\tilde P(\eta)\,d\eta,\qquad
\bar{\dot\omega}_\alpha=\int_0^1 \dot\omega_\alpha\big(Q(\eta),T_Q(\eta)\big)\tilde P(\eta)\,d\eta$$

**1st-order closure**: 条件付き反応率を条件付き平均で評価する $\langle\dot\omega|\eta\rangle\approx\dot\omega(Q(\eta),T_Q(\eta))$。
無条件平均での評価と違い、同じ $\eta$ に条件付けた後の残りのばらつきは小さいので良い近似になる (2nd-order は採らない)。

### 3. 条件付きモーメント方程式

密度重み無しの標準形 (Klimenko & Bilger 1999; RANS 実装は Patwardhan et al. 2009, Kim & Huh 等):

$$\frac{\partial Q_\alpha}{\partial t}+\langle\mathbf u|\eta\rangle\cdot\nabla Q_\alpha
=\nabla\cdot\big(D_t\nabla Q_\alpha\big)+\frac{\langle\chi|\eta\rangle}{2}\frac{\partial^2 Q_\alpha}{\partial\eta^2}
+\langle\dot\omega_\alpha|\eta\rangle/\bar\rho_\eta$$

- 条件付き速度は平均速度で近似 $\langle\mathbf u|\eta\rangle\approx\tilde{\mathbf u}$ (1st-order; Gordon 2007 の PDF 計算とも整合)。
- 物理空間の乱流拡散 $D_t=\nu_t/Sc_t$ (化学種と同じ $Sc_t$)。
- $\eta$ 空間の拡散項 (右辺第 2 項) が CMC の心臓部: $\chi=2D|\nabla\xi|^2$ はスカラー散逸率で、大きいほど
  $\eta$ 空間でプロファイルが均され、生まれかけたラジカルのピークが希釈されて着火が遅れる/消える。
- $\eta=0$ で酸化剤流 (coflow) の状態、$\eta=1$ で燃料流の状態を Dirichlet 境界とする。
- 条件付き密度 $\bar\rho_\eta=p/(R_{mix}(Q)T_Q)$。エネルギーは $Q_h$ を輸送し $T_Q$ を $h(T,Y)$ の反転で得る (thermo_d と同じ NASA-9)。

$\tilde P(\eta)\to0$ の領域で方程式が退化する問題は、この形 (密度・PDF 重み無し) では起きない。
密度重み形との差は $\nabla(\bar\rho\tilde P)$ に比例する項で、RANS 1st-order の精度では無視する (文献慣行)。

### 4. 混合分率の統計量と $\chi$ のモデル

- **前提: 等拡散**。CMC は全化学種の拡散係数が同一であることを仮定する。forge の混合則拡散 (`speciesDiffusionMethod: 1`) では
  ノズルリップの層流層で H₂ が速く拡散し、平均場が「Bilger $\tilde\xi$ に対応する混合線」から外れる (T で 130 K 以上、`run_0085`)。
  この平均場に PDF 積分状態をブレンドすると毎ステップ引き戻し合って圧力が励起され発散するので、**CMC を使う計算は
  `speciesDiffusionMethod: 0` (定数 Sc) で平均場も等拡散にする** (2026-09-05 決定)。
- **平均混合分率** $\tilde\xi$ は輸送しない。全化学種の拡散係数を同一 ($Sc_t$ 共通、分子拡散は $Sc$ 一定) としているので
  元素質量分率は保存スカラーであり、Bilger の定義
  $$\xi=\frac{\beta-\beta_O}{\beta_F-\beta_O},\qquad \beta=\frac{2Z_H}{W_H}-\frac{Z_O}{W_O}$$
  ($Z$: 元素質量分率、$F$/$O$ は燃料流/酸化剤流の値) を平均組成 $\tilde Y$ から**診断**する。C 系機構では
  $\beta=2Z_C/W_C+Z_H/(2W_H)-Z_O/W_O$。元素組成は機構 YAML の `species[].composition` から読む。
- **分散** $\widetilde{\xi''^2}$ は 1 本だけ輸送する (Jones–Launder 型):
  $$\frac{\partial\bar\rho\widetilde{\xi''^2}}{\partial t}+\nabla\cdot(\bar\rho\tilde{\mathbf u}\widetilde{\xi''^2})
  =\nabla\cdot\Big(\big(\mu+\tfrac{\mu_t}{Sc_t}\big)\nabla\widetilde{\xi''^2}\Big)
  +2\frac{\mu_t}{Sc_t}|\nabla\tilde\xi|^2-\bar\rho\tilde\chi$$
- **平均散逸率**: $\tilde\chi=C_\chi\,\beta^*\tilde\omega\,\widetilde{\xi''^2}$ ($\varepsilon/k=\beta^*\omega$, SST; $C_\chi=2$)。
- **条件付き散逸率**: AMC (Amplitude Mapping Closure, O'Brien & Jiang 1991)
  $$\langle\chi|\eta\rangle=\tilde\chi\frac{G(\eta)}{\int_0^1G(\eta)\tilde P(\eta)d\eta},\qquad G(\eta)=\exp\!\big(-2[\mathrm{erf}^{-1}(2\eta-1)]^2\big)$$
- **混合分率 PDF**: $\beta$ 分布 $\tilde P(\eta)=\eta^{a-1}(1-\eta)^{b-1}/B(a,b)$、
  $a=\tilde\xi\gamma,\ b=(1-\tilde\xi)\gamma,\ \gamma=\tilde\xi(1-\tilde\xi)/\widetilde{\xi''^2}-1$ ($\gamma$ は下限でクリップ)。
  分散が微小なセルはデルタ関数扱い ($\tilde Y=Q(\tilde\xi)$)。

### 5. 流れ側との結合

**採用 (couple 2, 2026-09-05)**: 流れソルバは $\tilde Y_\alpha,\tilde h$ を輸送し続けるが、反応ソースは「PDF 積分状態への緩和」で与える:

$$\bar{\dot\omega}_\alpha=\bar\rho\,\frac{\tilde Y^{pdf}_\alpha-\tilde Y_\alpha}{\tau},\qquad
\bar{\dot Q}=\bar\rho\,\frac{\tilde h^{pdf}-\tilde h}{\tau},\qquad
\tilde Y^{pdf}_\alpha=\sum_k\Omega_k Q_\alpha(\eta_k),\ \tilde h^{pdf}=\sum_k\Omega_k Q_h(\eta_k),\ \tau=c_r\Delta\tau_{local}$$

Jacobian は $\partial\bar{\dot\omega}_\alpha/\partial(\rho Y_k)=-\delta_{\alpha k}/\tau$ で点陰解に入る。$\sum_\alpha\bar{\dot\omega}_\alpha=0$ なので質量保存を壊さず、
定常では $\tilde Y\to\tilde Y^{pdf}$ となり、文献の RANS-CMC (平均スカラーは CMC から診断) と同じ状態に落ちる。

**不採用 (couple 1, source coupling)**: $\bar{\dot\omega}_\alpha=\sum_k\Omega_k\dot\omega_\alpha(Q(\eta_k))$ で平均方程式のソースを置換する方式は、
条件付き状態が燃え切ると $\dot\omega(Q)\to0$ となり、まだ燃えていない平均場に熱が渡らない (`run_0082`: 条件付き T 1615 K に対し平均 T_max
1123 K、$\tilde Y^{pdf}-\tilde Y$ が 0.07 で定常化)。両者の時間履歴が一致する保証がないため、積分値の整合を保てない。診断用に残す。

#### (旧) source coupling の記述

流れソルバは従来どおり $\tilde Y_\alpha$ と $\tilde h$ を輸送し、**反応ソースだけ**を条件付き平均から作る:

$$\bar{\dot\omega}_\alpha=\sum_k w_k\tilde P(\eta_k)\,\dot\omega_\alpha(Q(\eta_k),T_Q(\eta_k)),\qquad
\bar{\dot Q}=\sum_k w_k\tilde P(\eta_k)\,\dot Q(Q(\eta_k),T_Q(\eta_k))$$

これは laminar chemistry の $\dot\omega(\tilde Y,\tilde T)$ を置き換えるだけなので、質量保存・既存の陰解法配管
(種ブロック LU、案C 反応熱注入) をそのまま使える。剛性は CMC 側 ($\eta$ 各点の点陰解) で吸収し、平均方程式側の
ソース Jacobian は「同じ $\eta$ 点での $\partial\dot\omega/\partial\rho Y$ の PDF 加重和」で対角近似する。
$\tilde Y$ の PDF 積分値と輸送値のずれは、定常では両者が同じ $Q$ から作られるので小さい (診断で監視する: `cmc_dY`)。
$\tilde Y$ を PDF 積分で上書きする「full coupling」は採らない (質量保存を流れ側に残すため)。

### 6. $\eta$ 格子と離散化

- $\eta_k$ は $N_\eta$ 点 (既定 41)、$\xi_{MR}$ (超希薄側) と $\xi_{st}$ に密に取る (tanh 伸長)。
- $\eta$ 拡散は 2 次中心差分、陰的 (三重対角)。物理空間の移流は 1 次風上 (条件付き量は滑らか)、拡散は汎用 `scalarTransport_d`。
- 化学は各 $\eta$ 点で `chem_source` を評価し、点陰解 (種ブロック LU, `jacobianMode 2` と同じ機構) を $\eta$ 点ごとに行う。
  `jacobianInterval` の凍結も同様に使う。
- 疑似時間: 流れと同じ局所 $\Delta\tau$。$\eta$ 方向は完全陰的なので $\chi$ の剛性は制約にならない。

### 参考文献

- Klimenko & Bilger, Prog. Energy Combust. Sci. 25 (1999) — CMC の定式。
- Patwardhan, De, Lakshmisha, Raghunandan, Proc. Combust. Inst. 32 (2009) — Cabra H₂/N₂ の RANS 1st-order CMC。
- Kim, Huh (Combust. Flame 2002) — RANS-CMC の数値実装 (AMC, β-PDF)。
- O'Brien & Jiang, Phys. Fluids A 3 (1991) — AMC。
- Mastorakos, Prog. Energy Combust. Sci. 35 (2009) — 自着火の物理と各クロージャの比較。

---

## 実装

### 1. モジュール構成 (設計)

| 処理 | ファイル | 内容 |
| --- | --- | --- |
| 元素組成 | `chemistry_mech_io.hpp` | `species[].composition` を読み、元素×種の行列 $a_{e,s}$ と Bilger 係数を `ReactionTable` に持つ |
| 混合分率診断 | `cuda_forge/cmc_d.cu` (`cmc_mixfrac_d`) | $\tilde Y$ → $\tilde\xi$ (Bilger)。`xi` を出力変数に登録 |
| 分散輸送 | `cuda_forge/cmc_d.cu` + `scalarTransport_d` | `roXiVar` (保存量), `xiVar`, `res_roXiVar`, `transport_diag_xiVar`, `src_jac_xiVar`。生成 $2\mu_t/Sc_t|\nabla\tilde\xi|^2$・散逸 $\bar\rho\tilde\chi$ を点陰解 (散逸は `src_jac`) |
| 条件付きスカラー | `cuda_forge/cmc_d.{cuh,cu}` | `Q` を `[slice][node]` (slice = var·$N_\eta$+k, `flow_float`) で保持。物理空間輸送の残差は `scalarTransport_d` を $N_\eta\times(n_s+1)$ 回、時間積分は全スライス一括 (`cmc_q_timeint_d`)。$\eta$ 拡散 + 化学は 5 カーネルに分割 (下記 §4): (A) node 毎 β-PDF 重み・AMC、(B) (node,var) 毎 Thomas、(C) node 毎 PDF 平均、(D) (node,η) 毎の化学 (活性対を圧縮リスト化; 速度定数・Jacobian は fp32 `chemistry_f32_d.cuh`、9×9 点陰解は double)、(E) node 毎の診断・ρ̄Q 再同期 |
| β-PDF・積分 | `cmc_d.cuh` (`cmc_beta_pdf`, `cmc_integrate`) | Gauss–Jacobi でなく $\eta$ 格子上の台形則 + β-PDF の解析正規化 (端の特異性は不完全 β 関数で処理) |
| ソース結合 | `chemistry_d.cu` (`chemistry_source_d` の `cmcMode` 分岐) | couple 2: $\bar\rho(\tilde Y^{pdf}-\tilde Y)/\tau$ 緩和 (既定候補), couple 1: PDF 平均 $\bar{\dot\omega}$ (不採用・診断) |
| 設定 | `physProp.chemistry.cmc: {enabled, nEta, etaPow, pdfFloor, couple (1/2), chem, fuelT, oxidizerT, dtScale, interval, relax, xiSt}` | [solver-settings.md](../procedures/solver-settings.md) |
| 出力 | `xi`, `xiVar`, `chi`, 診断 `cmc_dY` (PDF 積分 $\tilde Y$ と輸送 $\tilde Y$ の最大差), `cmc_TQmax` (η 上の最高 T), `cmc_TQst` ($\eta\approx\xi_{st}$ の T) | `res_*.h5` |

### 2. 境界条件

- 混合分率分散: 入口 0、壁 zero-gradient、出口 zero-gradient (汎用スカラの既定)。
- 条件付きスカラー: $\eta$ 端は Dirichlet (燃料流/酸化剤流の入口状態を `physProp.chemistry.cmc.fuel/oxidizer` または入口 bcond の組成・温度から自動取得)。
  物理空間の入口は「その入口の $\xi$ に対応する混合線状態」、壁・出口は zero-gradient。

### 3. 検証 (plan §6)

1. **凍結混合の整合**: 化学 OFF で $Q_\alpha(\eta)$ が混合線 (線形) に一致し、PDF 積分 $\tilde Y$ が輸送 $\tilde Y$ と一致 (`cmc_dY` < 1 %)。
2. **Cabra H₂/N₂** (case/48): $T_c$ = 1015/1030/1045/1060/1075 K の浮き上がり高さ $H(T_c)$ (Y_OH > 2e-4) を
   実験 (Cabra 2002, Wu 2003) と PDF/CMC 文献値 (Gordon 2007, Patwardhan 2009) に重ね、応答曲線の傾きと
   $T_c\pm30$ K 内での一致を判定。中心軸・半径 T の平均差 ±30 K (実験読み取り誤差) 以内を目標。
3. **Burrows–Kurkov** (case/47): 着火位置 18–25 cm。
4. **回帰**: `tci: 0/1` はビット不変。

### 4. 性能 (2026-09-05, Cabra 60k node, $N_\eta$=41, Li 9 種 21 反応, RTX 3060)

初版の node 毎 1 カーネル (`cmc_step_d`, 1 thread が 41 η 点の化学を直列) は **CMC 部分だけで ~650 ms/step** (混合のみ 61 ms) だった。
`CMC_TIMING=<steps>` (環境変数; cudaEvent 区間計測) と ncu で律速を特定し、以下の順に詰めた (数値は run_0098 と GPU 共用中の相対値。
場は各段階で double 参照と相対 $10^{-6}$–$10^{-5}$ で一致することを 40 step の A/B で確認)。

| 段階 | 変更 | `cmc_step` [ms/step] |
| --- | --- | --- |
| 初版 | node 毎 1 カーネル | 635 |
| 分割 (A–E) | (node,η) thread 化 | 608 (占有率は律速でなかった) |
| 活性対の圧縮 | フラグ → `thrust::copy_if` で (node,η) リスト化、ワープ発散を除去 (活性 28 %) | 529 |
| T 反転のウォームスタート | 前 step の $T(\eta)$ を Newton 初期値に、$\ln T$ を種ループ外へ | 384 |
| **fp32 化学** (`cmc.fp32: 1`, 既定) | ncu: FP64 パイプ 87 % が律速 (GA10x は fp64 = fp32 の 1/64)。速度定数・Jacobian・T 反転を float (`chemistry_f32_d.cuh`; 係数表は `__constant__`)、点陰解 Gauss は double | 162 |
| η 拡散 float + 時間積分の一括化 | Thomas を float (対角優位)、410 スライスの point-implicit 更新を 1 カーネル | **133** (+ 残差 41) |

平均場の化学 (`chemistry_d.cu`) は double のまま。fp32 化の精度根拠: 速度定数は対数空間で評価するので範囲は足りる、
$\omega=k_fP_f-k_bP_b$ の平衡近傍の桁落ちは $10^{-7}\max(q_f,q_b)$ で点陰解が減衰する、$g_s$ の絶対誤差 $10^{-5}$ は $K_c$ の相対 $4\times10^{-5}$。
残る大物は物理空間輸送の残差 (410 スライス × 2 カーネル, 41 ms) で、面幾何を共有する融合カーネル化が次の候補 (plan §5.1)。
旧カーネルは `CMC_STEP_LEGACY=1` で残してある (A/B 参照用)。


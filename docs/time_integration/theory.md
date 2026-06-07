# 時間積分 — 理論

forge の時間積分は次の 3 系統を提供する。

| `cfg.timeIntegration` | 概要 |
| --- | --- |
| `1` / `3` | Jameson 型多段陽解法 (擬時間ステップ) |
| `4` | 4 段 4 次 Runge–Kutta (低 storage、係数は `coef_DT_4thRunge`, `coef_Res_4thRunge`) |
| `11` | 陰解法 (block DPLUR、LU-SGS $A^\pm$)。`blockDPLUR`=1(5×5 ブロック, 既定)/0(スカラー対角)。定常・非定常 dual-time とも実装済 |

> **実装状態 (2026-06)**: `timeIntegration == 11` は **block DPLUR (`blockDPLUR == 1`)・定常 (`unsteady == 0`)** のみを正式サポートする。
> スカラー対角版 (`blockDPLUR == 0`) も有効（block より低 `cfl_pseudo`・低速）。非定常 dual-time (`unsteady=1 && dualTime=1`) も実装済（`blockDPLUR=1` のみ、物理 $\Delta t$ 固定）。

## 共通枠組み

時間 1 ステップは内部ループ (`loop = 0 .. nInner-1`) に分かれ、
1 ループあたり次の手順を踏む。

```
applyBconds → gradient → limiter → ducrosSensor
   → convectiveFlux + viscousFlux  (res_* に累積)
   → implicitCorrection (任意)
   → timeIntegration_d_wrapper (Q を更新)
   → dependentVariables (派生量)
```

`res_*` は対流・粘性両方からの符号付き面寄与の合計。これを各セルで体積で
正規化したものが $-\partial \mathbf{Q}/\partial t \cdot V$ に相当する。

## 陽解法 (Jameson 多段)

低 storage 形式: 1 ステージあたり

$$
\mathbf{Q}^{(k+1)} = c_N \mathbf{Q}^N + c_M \mathbf{Q}^M
+ c_R \frac{\Delta t_{\text{loc}}}{V}\,\mathbf{R}^{(k)},
$$

ここで $\mathbf{Q}^N$ は前回外部ステップ、$\mathbf{Q}^M$ は前回内部ステージの値。
係数 $c_N, c_M, c_R$ は `solverConfig` の `coef_N`, `coef_M`, `coef_Res` に格納する
(段数 = `cfg.nInner`)。

定常計算では局所時間刻み $\Delta t_{\text{loc}}$ を使うことで収束を加速する。

## 4 段 4 次 Runge–Kutta (low-storage)

参考: Allampalli et al. (2009)。`coef_DT_4thRunge`, `coef_Res_4thRunge` を用い、
内部 4 段で

$$
\mathbf{Q}^{(k+1)} = \mathbf{Q}^N + c_{DT}^{(k)} \frac{\Delta t_{\text{loc}}}{V}\,\mathbf{R}^{(k)},
$$

累積残差 $\mathbf{R}_m = \sum_k c_R^{(k)} \mathbf{R}^{(k)} \,\Delta t_{\text{loc}}/V$ を保持し、
最終段で $\mathbf{Q}^{N+1} = \mathbf{Q}^N + \mathbf{R}_m$ とする。LES など時間精度の必要な解析向け。

## 陰解法 block DPLUR (LU-SGS 型 $A^\pm$ 分割)

擬似時間 $\tau$ を導入し、$\Delta \mathbf{Q} = \mathbf{Q}^{n+1} - \mathbf{Q}^n$ を未知に取る。
1 次精度 upwind フラックス $\mathbf{F}_f = \tfrac12(\mathbf{F}_i+\mathbf{F}_j) - \tfrac12|\widetilde{A}_f|(\mathbf{Q}_j-\mathbf{Q}_i)$ の
線形化から、残差 $\mathbf{R}_i=\sum_f \mathbf{F}_f S_f$ の Jacobian は**通量分割行列**

$$
A^{\pm} = \tfrac12\big(\widetilde{A} \pm |\widetilde{A}|\big) = R\,\Lambda^{\pm} L,\qquad
\Lambda^{+}=\max(\Lambda,0),\ \ \Lambda^{-}=\min(\Lambda,0)
$$

で表される。これを使うと **対角ブロックは upwind スキームの厳密な自己 Jacobian** $\partial\mathbf{R}_i/\partial\mathbf{Q}_i=\sum_f A^{+}_f S_f$、
近傍ブロックは $\partial\mathbf{R}_i/\partial\mathbf{Q}_j = A^{-}_f S_f$ となり、解くべき線形系は

$$
\underbrace{\left[\frac{V}{\Delta \tau} I + \sum_f A^{+}_f S_f + \sum_f \rho^{\nu}_f S_f\, I\right]}_{D_i\ (\text{対角ブロック})} \Delta \mathbf{Q}_i
= -\mathbf{R}_i(\mathbf{Q}^n) - \sum_f A^{-}_f S_f\, \Delta \mathbf{Q}_{\text{nbr}}
$$

となる（$\rho^{\nu}_f$ は粘性スペクトル半径による対角強化）。ここで $\widetilde{A}_f=R\Lambda L$ は
セル中心状態から作る面法線フラックス Jacobian（固有値 $\Lambda=\{U{+}c,\,U,\,U,\,U,\,U{-}c\}$）で、$R\Lambda L$ は
真の流束 Jacobian を厳密に再現し、$|\widetilde{A}|=R|\Lambda|L$ の固有値は全て非負（matrix-free・近似 BC）。

> **2026-06 修正（重要）**: 旧実装は対角・近傍ともに**絶対値 Jacobian $|\widetilde{A}|$ を符号付き $A^\pm$ の代わりに誤用**しており、
> 対角が upwind 自己 Jacobian と一致せず（対流中心項を欠く）、近傍結合の符号も誤っていた。このため block DPLUR は
> 新旧コードとも **一度も収束せず発散**していた（対角のみ＝頭打ち、近傍 sweep 有効＝NaN）。
> $A^{+}$ 対角・$A^{-}$ 近傍の正しい LU-SGS 分割に置換することで反復行列のスペクトル半径が
> $\rho \approx \dfrac{V/\Delta\tau}{V/\Delta\tau+\lambda^{+}} < 1$ となり収束する。詳細は [implementation.md](implementation.md)。

### 古典 DPLUR 反復 (残差固定 + 複数 sweep)

非線形残差 $\mathbf{R}(\mathbf{Q}^n)$ と対角ブロック $D_i$ を**1 回の非線形更新あたり 1 度だけ構築**し、
$\mathbf{Q}^n$ を固定したまま $\Delta\mathbf{Q}$ を Jacobi 型に複数回緩和する（lagged neighbor）。

$$
\Delta\mathbf{Q}_i^{(s+1)} = D_i^{-1}\Big(-\mathbf{R}_i - \sum_f A^{-}_f S_f\, \Delta\mathbf{Q}_{\text{nbr}}^{(s)}\Big),
\quad s = 0 \ldots n_{\text{sweep}}-1
$$

各 sweep で隣接セルの $\Delta\mathbf{Q}$ は前 sweep の値（`dq_old`）を参照し（Gauss–Seidel ではなく Jacobi）、
sweep 間で `dq_old`/`dq_new` バッファを入れ替える。全 sweep 後に **1 度だけ** $\mathbf{Q} = \mathbf{Q}^n + \Delta\mathbf{Q}$ を commit する
（sweep 中は $\mathbf{Q}$ を固定するのが古典 DPLUR の要件）。`implicitRelax` で $\Delta\mathbf{Q}$ を緩和する。

`build_jacobian_split` がセル中心状態から $A^{+}$ と RHS 用の $-A^{-}=\tfrac12(|\widetilde{A}|-\widetilde{A})$ を同時に構築、
`solve_5x5` が部分ピボット付き Gauss 消去で $D_i^{-1}(\cdots)$ を解く（ピボット `< 1e-20` でゼロ解にフェイルセーフ）。
擬似 CFL（`cfl_pseudo`）を上げるほど $V/\Delta\tau$ が小さくなり $\rho$ が下がる（収束加速）。

### 反復回数の意味付け

| config | 意味 |
| --- | --- |
| `nStepInner` | **DPLUR 緩和 sweep 回数** $n_{\text{sweep}}$（1 非線形更新内の線形反復） |
| `nStepOuter` | 定常では擬似時間ステップ数（メインループ＝擬似時間行進） |

定常計算はメインループ（擬似時間）を回し、1 ステップあたり 1 回の非線形更新
（残差構築 → DPLUR sweep ×`nStepInner` → commit）を行う。

### 低マッハ前処理固有値 (Weiss–Smith) — 試行したが不採用

> **結論: block DPLUR では有害のため不採用**（2026-06 検証。実装後 revert）。理論メモとして残す。

着想は、$A^\pm$ 分割に使う固有値を基準速度 $U_r=\min(c,\max(|\mathbf u|,\epsilon c))$、$\beta=(U_r/c)^2$ で
前処理した

$$
\lambda_{1,2,3}' = U_n,\qquad
\lambda_{4,5}' = \tfrac12(1+\beta)U_n \pm \tfrac12\sqrt{(1-\beta)^2 U_n^2 + 4U_r^2}
$$

に差し替え、低マッハの音響 stiffness を除く、というもの。**しかし block DPLUR (LU-SGS $A^\pm$) では
これが逆効果**だった。理由: 対角 $D_i=V/\Delta\tau\,I+\sum_f A^{+}_f S_f$ の**対角優位性の源は大きい音響
固有値 $U_n\pm c$ を含む $A^{+}$** であり（2026-06 の収束化はこれに依拠）、$U_n\pm c'$ (小) に落とすと
$A^{+}$ が縮んで対角優位性が下がる。検証では、フラックス散逸前処理単独 ($\epsilon=0.15$) では安定だった
ものが、固有値前処理を加えると発散し、**安定 $\epsilon$ 範囲をむしろ狭めた**（収束加速も根治もなし）。
あわせて検討した setDT のスペクトル半径前処理 ($c\to c'$ による $\Delta t_{\rm loc}$ 増大) も、
$V/\Delta\tau$ を縮めて対角を二重に潰し、最速で発散した。

根因は、`build_jacobian_split` が Euler 固有ヤコビアンで**SLAU 低マッハ散逸項を表さない**こと。固有値だけ
前処理しても aggressive なフラックス前処理を陰的に安定化できない。真の安定化には固有ベクトルも前処理する
**完全 preconditioned-flux Jacobian**（大改修）が要る。よって現状はフラックス散逸前処理のみ採用し、LHS は
従来の $A^\pm$ のまま。詳細・データは計画
[`time_integration-lowmach-preconditioning.md`](../../.github/plans/time_integration-lowmach-preconditioning.md) §9。

## 非定常 dual-time 陰解法 (実装済み 2026-06)

非定常は物理時間 $t$ に対し BDF（後退差分）で物理時間微分項を残差に加え、各物理ステップ内で
擬似時間サブ反復により残差を収束させる二重時間刻み法を用いる。残差は

$$
\mathbf{R}^*_i = \mathbf{R}_i(\mathbf{Q}) + \frac{V}{\Delta t}\big(a\mathbf{Q} - b\mathbf{Q}^{n} + c\mathbf{Q}^{n-1}\big)
$$

の形（BDF1: $a{=}1,b{=}1,c{=}0$／BDF2: $a{=}\tfrac32,b{=}2,c{=}\tfrac12$）で、物理時間項は対角 $D_i$ にも $aV/\Delta t$ として加わる。
第 2 時間レベル $\mathbf{Q}^{n-1}$ は `roNN` 系に保持する。擬似時間サブ反復・DPLUR sweep の構造は定常と共有する
（3 階層: 物理ステップ × 擬似時間サブ反復 × DPLUR sweep。定常は前 2 つが縮退）。

| config | 意味 |
| --- | --- |
| `nStepOuter` | 物理時間ステップ数 |
| `nSubIterDualTime` | 物理ステップ内の擬似時間サブ反復数（既定 20） |
| `bdfOrder` | 物理時間 BDF 次数 (1 or 2、初回ステップは自動的に BDF1) |
| `nStepInner` | 各サブ反復内の DPLUR sweep 回数（定常と同義） |

使用条件: `unsteady=1, dualTime=1, timeIntegration=11, blockDPLUR=1, time.deltaT.control=0`（物理 $\Delta t$ 固定）。
擬似時間 $\Delta\tau$ は `cfl_pseudo` から決まり（物理 $\Delta t$ とは独立）、物理 $\Delta t$ は固定。
commit は in-place（$\mathbf Q \leftarrow \mathbf Q + \delta\mathbf Q$。`roN`=$\mathbf Q^n$ は BDF 基準で固定のため）。

> **検証 (2026-06, `case/20.naca_ml`)**: 各物理ステップ内で擬似サブ反復が $\mathbf R^*$ を ~4 桁低減
> （例: rms_roe 22.9→0.0035, サブ反復 ~10 で収束）。陽解法 CFL 上限を超える $\Delta t$ でも安定
> （$\Delta t{=}2\times10^{-6}$, 陽解法安定限界 ~$9\times10^{-7}$）。同一 $\Delta t{=}5\times10^{-7}$ では陽解法
> 時間精度解と壁面静圧が一致（平均 0.006%、最大 0.2% は RK3 と BDF2 の時間スキーム差 $O(\Delta t^2)$）。

## スカラー (k/ω) の陰解法 (segregated point-implicit)

平均流 5 式とスカラー (k/ω) は常に別カーネルで分離 (segregated) して積分する。
block 陰解法 (`timeIntegration==11`) では、平均流 block DPLUR を解いて commit した後、
**同じ擬似時間ステップで k/ω を segregated point-implicit で更新**する（旧実装はここで k/ω を凍結していた）。

SST ソースの消散項は $\mathbf{Q}$ に比例する負ソースで、時間刻みに依らず stiff:

$$
D_k = \beta^\* \rho k \omega = \beta^\*(\rho k)\,\omega,\qquad
D_\omega = \beta \rho \omega^2 = \beta \frac{(\rho\omega)^2}{\rho}
$$

$$
\frac{\partial D_k}{\partial(\rho k)} = \beta^\*\omega,\qquad
\frac{\partial D_\omega}{\partial(\rho\omega)} = 2\beta\omega
$$

これを対角に陰に取り込む point-implicit 更新を行う（移流・拡散・生産項は残差 $\text{res}_{\rho\phi}$ に
含む lagged 扱い。生産項は limiter 付きで bounded のため陰化不要）:

$$
D_\phi = \frac{V}{\Delta\tau} + V\,\frac{\partial D}{\partial(\rho\phi)},\qquad
\delta(\rho\phi) = \mathrm{relax}\cdot\frac{\text{res}_{\rho\phi}}{D_\phi},\qquad
\rho\phi \leftarrow \max(\rho\phi^{N} + \delta(\rho\phi),\ \text{floor})
$$

擬似時間 $\Delta\tau$ は平均流と同じ `dt_local` を流用する。realizability のため $\rho k \ge 0$、$\rho\omega > 0$ を課す。
近傍結合は持たない純 point-implicit（消散の stiff 性のみを陰化する最小構成）。これにより block 陰解法でも
k/ω が凍結せず収束し、乱流ケースを陰解法で回せる。

### 陽解法 RK での point-implicit 源項 (2026-06)

陽解法 RK (`timeIntegration==1/3`) で RANS (SST) を回す場合、平均流は陽的でも
**スカラー (k/ω) の消散項だけは上と同じ消散ヤコビアンで point-implicit 化**する必要がある。
純陽的に積分すると、壁近傍の大きな $\omega$ と $D_\omega=\beta\rho\omega^2$ の stiff 性により
CFL を下げても（$\mathrm{CFL}=0.05$ でも）初回サブステップで $\omega$ が発散する（陽的安定限界が
源項の時間スケール $1/(2\beta\omega)$ で律速されるため）。

陽解法 RK のスカラー更新は本来

$$
\rho\phi \leftarrow c_N\,\rho\phi^{N} + c_M\,\rho\phi^{M}
  + c_{\text{res}}\,\frac{\Delta t_{\text{loc}}}{V}\,\text{res}_{\rho\phi}
$$

だが、残差増分のみを消散ヤコビアンで割って減衰させる:

$$
\rho\phi \leftarrow c_N\,\rho\phi^{N} + c_M\,\rho\phi^{M}
  + \frac{1}{1 + c_{\text{res}}\,\Delta t_{\text{loc}}\,\partial D/\partial(\rho\phi)}
    \cdot c_{\text{res}}\,\frac{\Delta t_{\text{loc}}}{V}\,\text{res}_{\rho\phi}
$$

ここで $\partial D/\partial(\rho\phi)\ (=\beta^\*\omega,\ 2\beta\omega)$ は block 陰解法と同じ
`src_jac_k`/`src_jac_omega`。減衰係数 $1 + c_{\text{res}}\Delta t_{\text{loc}}\,\partial D/\partial(\rho\phi)\ge 1$
は陰解法対角 $D_\phi = V/\Delta\tau + V\,\partial D/\partial(\rho\phi)$ に $\Delta\tau/V$ を掛けた形と整合する
（消散項のみを陰化し、移流・拡散・生産は lagged のまま）。これにより陽解法 RK でも RANS が安定に回り、
平均流の時間積分法 (陽 RK / block 陰) を揃えた比較ができる。realizability の floor は陰解法と共通。
4 段 4 次 RK (`timeIntegration==4`) のスカラー積分は本減衰を未適用（RANS 非対応のまま）。

## 局所時間刻み

`setDT_d_wrapper` (`cuda_forge/setDT_d.cu`) が各セルの $\Delta t_{\text{loc}}$ を
CFL 条件から決定する。スペクトル半径

$$
\lambda_{\max} = \sum_f \big( |U_n| + c \big)_f A_f + \text{粘性項}
$$

から $\Delta t_{\text{loc}} = \mathrm{CFL}\,V / \lambda_{\max}$ を求める。
`cfg.dt` が指定された場合はそれを上限とする。

> **setDT の低マッハ前処理は不採用** (2026-06 検証)。スペクトル半径の $c\to c'$ 置換は低マッハ域で
> $\Delta t_{\text{loc}}$ を増大させ陰解法対角 $V/\Delta\tau$ を縮める。block DPLUR では対角優位性が崩れて
> 発散するため採用しない（LHS 固有値前処理と併用しても二重に対角が潰れさらに悪化。上記「低マッハ前処理
> 固有値」節を参照）。低マッハ前処理は対流フラックスの散逸スケールにのみ適用する。

## 参考

- 各カーネル分岐の詳細: [implementation.md](implementation.md)。
- 境界条件適用 (`applyBconds`): [`docs/boundary/`](../boundary/)。
- 残差を作る対流・粘性: [`docs/convection/`](../convection/), [`docs/diffusion/`](../diffusion/)。

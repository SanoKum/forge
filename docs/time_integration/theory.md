# 時間積分 — 理論

forge の時間積分は次の 3 系統を提供する。

| `cfg.timeIntegration` | 概要 |
| --- | --- |
| `1` / `3` | Jameson 型多段陽解法 (擬時間ステップ) |
| `4` | 4 段 4 次 Runge–Kutta (低 storage、係数は `coef_DT_4thRunge`, `coef_Res_4thRunge`) |
| `11` | 陰解法 (block DPLUR、簡易 Jacobian)。現状 `blockDPLUR == 1`（5×5 ブロック）のみ有効。定常を実装、非定常 dual-time は後続フェーズ |

> **実装状態 (2026-06)**: `timeIntegration == 11` は **block DPLUR (`blockDPLUR == 1`)・定常 (`unsteady == 0`)** のみを正式サポートする。
> スカラー対角版 (`blockDPLUR == 0`) は一時的に無効化 (config で reject)、非定常 dual-time (`unsteady == 1`) は構造のみ用意し本体は未実装 (reject)。いずれも後続フェーズで有効化する。

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

## 非定常 dual-time 陰解法 (後続フェーズ)

> 本体は未実装。定常がその縮退形になるよう構造のみ用意する。

非定常は物理時間 $t$ に対し BDF（後退差分）で物理時間微分項を残差に加え、各物理ステップ内で
擬似時間サブ反復により残差を収束させる二重時間刻み法を用いる。残差は

$$
\mathbf{R}^*_i = \mathbf{R}_i(\mathbf{Q}) + \frac{V}{\Delta t}\big(a\mathbf{Q} - b\mathbf{Q}^{n} + c\mathbf{Q}^{n-1}\big)
$$

の形（BDF1: $a{=}1,b{=}1,c{=}0$／BDF2: $a{=}\tfrac32,b{=}2,c{=}\tfrac12$）で、物理時間項は対角 $D_i$ にも $aV/\Delta t$ として加わる。
第 2 時間レベル $\mathbf{Q}^{n-1}$ は `roNN` 系に保持する。擬似時間サブ反復・DPLUR sweep の構造は定常と共有する
（3 階層: 物理ステップ × 擬似時間サブ反復 × DPLUR sweep。定常は前 2 つが縮退）。

| config (後続フェーズで追加) | 意味 |
| --- | --- |
| `nStepOuter` | 物理時間ステップ数 |
| `nSubIterDualTime` | 物理ステップ内の擬似時間サブ反復数 |
| `nStepInner` | DPLUR sweep 回数（定常と同義） |

## スカラー (k/ω) の扱い

平均流 5 式とスカラー (k/ω) は常に別カーネルで分離 (segregated) して積分する。
block 陰解法中はスカラー時間積分が休止するため **k/ω は凍結**される（frozen scalar の本格対応は後続）。
このため陰解法の検証は層流（スカラー無し）で行う。

## 局所時間刻み

`setDT_d_wrapper` (`cuda_forge/setDT_d.cu`) が各セルの $\Delta t_{\text{loc}}$ を
CFL 条件から決定する。スペクトル半径

$$
\lambda_{\max} = \sum_f \big( |U_n| + c \big)_f A_f + \text{粘性項}
$$

から $\Delta t_{\text{loc}} = \mathrm{CFL}\,V / \lambda_{\max}$ を求める。
`cfg.dt` が指定された場合はそれを上限とする。

## 参考

- 各カーネル分岐の詳細: [implementation.md](implementation.md)。
- 境界条件適用 (`applyBconds`): [`docs/boundary/`](../boundary/)。
- 残差を作る対流・粘性: [`docs/convection/`](../convection/), [`docs/diffusion/`](../diffusion/)。

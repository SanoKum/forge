# リミッタ

MUSCL 再構成 (2 / 3 次) は不連続近傍で振動を生む。forge ではセル中心勾配にスカラ係数
$\psi \in [0, 1]$ を掛けて再構成を抑制する。

$$
\phi_L = \phi_C + \psi_C \,\nabla \phi_C \cdot (\mathbf{r}_f - \mathbf{r}_C).
$$

本ドキュメントは理論(係数の定義とスキーム)と実装(ソース対応)をまとめる。

## 理論

### 評価対象

セルごとに $\rho, U_x, U_y, U_z, P$ の 5 成分について独立にリミッタ係数を持つ。
温度 $T$ と全エンタルピ $H_t$ は派生量として扱い、独立リミッタは持たない。

### $\delta^+, \delta^-, \delta_m$ の定義

セル $C$ について隣接セル値の振れ幅を

$$
\delta^+_C = \phi_C^{\max} - \phi_C ,\qquad
\delta^-_C = \phi_C^{\min} - \phi_C
$$

($\phi_C^{\max/\min}$ は近傍セル値の最大・最小)。
面 $f$ への線形再構成による予測差分を

$$
\delta_m = \tilde\phi_f - \phi_C
= \nabla \phi_C \cdot (\mathbf r_f - \mathbf r_C)
$$

とする。$\delta_m > 0$ なら $\delta^+$ を、$\delta_m < 0$ なら $\delta^-$ を比較対象に取る。

### Barth–Jespersen リミッタ

$$
\psi_f = \min\!\left(1,\, \frac{\delta^\pm}{\delta_m}\right),
\qquad \psi_C = \min_{f \in \partial C} \psi_f .
$$

単純で TVD 性が強く、滑らかな領域でも 2 次精度の頭打ちを招きやすい。

### Venkatakrishnan リミッタ

Barth–Jespersen を滑らか化したもの。$\epsilon^2 = (K |V_C|^{1/3})^3$ (`K=1.0`) として

$$
\psi(\delta^+, \delta_m) =
\frac{1}{\delta_m}\,
\frac{(\delta^{+\,2} + \epsilon^2)\,\delta_m + 2\delta_m^2 \delta^+}
     {\delta^{+\,2} + 2\delta_m^2 + \delta^+ \delta_m + \epsilon^2}.
$$

$\delta_m \to 0$ で $\psi \to 1$ (連続)。$|\delta_m|$ が体積スケール以下では実質
無リミットになり、滑らかな領域での精度低下を避ける。

### Nishikawa R1 リミッタ (未有効)

`nishikawa_r1_limiter` の実装は残されているがコメントアウト済み。
$\delta'_C$ と幾何比 $r_{ik}$ から不連続適応する形式。

### Ducros センサとの併用

[`methods/convection/`](convection/) で記述したように、
Roe カーネルでは Ducros センサ値 $\mathrm{duc}$ で

$$
\psi \leftarrow \begin{cases}
\psi & (\mathrm{duc} \le 0.8) \\
\max(0,\, (1 - \mathrm{duc})\,\psi) & (\mathrm{duc} > 0.8)
\end{cases}
$$

と二段がけする。

### 全リミッタ 1 の場合 (`cfg.limiter == 0`)

リミッタを評価せず $\psi = 1$ を一律で配るモード。低次再構成 (`convMethod = 0`)
や smooth な計算ですべての面で線形再構成を使いたい場合に選ぶ。

## 実装

forge のスカラ・リミッタ計算は GPU でのみ動作する (CPU 経路は無い)。

### ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/cuda_forge/limiter_d.cuh`](../solver_density_cuda/cuda_forge/limiter_d.cuh) | カーネル / ラッパ宣言、`deltas` 構造体 |
| [`solver_density_cuda/cuda_forge/limiter_d.cu`](../solver_density_cuda/cuda_forge/limiter_d.cu) | リミッタ device 関数とセル並列カーネル (約 325 行) |

### エントリポイント

[`limiter_d_wrapper`](../solver_density_cuda/cuda_forge/limiter_d.cu#L198)
が共通入口。次の順で動く。

1. `fill_limiter_d` で `limiter_ro, limiter_Ux, limiter_Uy, limiter_Uz, limiter_P` を
   `1.0` に初期化 (全セル `nCells_all`)。
2. `cfg.limiter == 0` なら早期 return (全リミッタ 1)。
3. 各成分について `limiter_r1_d` をセル並列で呼び、近傍探索を行って $\psi$ を確定。

呼び出しは時間積分ループの `gradient → limiter → convectiveFlux` の中間で発生する
([`main.cpp` L700](../solver_density_cuda/main.cpp#L700), [L730](../solver_density_cuda/main.cpp#L730))。

### `limiter_r1_d` 構造 ([L75](../solver_density_cuda/cuda_forge/limiter_d.cu#L75))

セル 1 並列。`ic0 < nCells` で動作する。

1. `limiter_scheme == 0` なら `limiter_Q[ic0] = 1.0` で抜ける。
2. `cell_planes_index[ic0]`, `cell_planes[*]` でセルが接する面を列挙。
3. 各内部面 (`ip < nNormalPlanes`) について反対側セル `ic1` を取得 (`plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0`)、
   `Q_max, Q_min` を更新。
4. もう一度ループして面ごとに $\delta_m$ を求め、`limiter_function` で $\psi_f$ を計算、
   最小値を取って `limiter_Q[ic0]` に書く。
5. 結果は `clamp(ψ, 0, 1)` で安全側に丸める。

#### `limiter_function` の選択

`limiter_scheme` (= `cfg.limiter`) で device 関数ポインタを差し替える。

| `limiter_scheme` | 関数 | 種別 |
| --- | --- | --- |
| `0` | — (早期 return) | 無リミッタ ($\psi = 1$) |
| `1` | `barth_Jespersen_limiter` ([L39](../solver_density_cuda/cuda_forge/limiter_d.cu#L39)) | TVD 強 |
| `2`, `-1` | `venkata_limiter` ([L17](../solver_density_cuda/cuda_forge/limiter_d.cu#L17)) | 滑らか版 |
| `3` (TODO) | `nishikawa_r1_limiter` | 実装あり・コメントアウト |

`venkata_limiter` の $\epsilon^2$ は `K = 1.0` 固定 (`K^3 * volume`)。
`|δ_m| < 1e-20` のとき $\psi = 1$ にフォールバック。

### 適用先 (面ごとの参照)

各スキームカーネル ([`convectiveFlux_d.cu`](../solver_density_cuda/cuda_forge/convectiveFlux_d.cu))
は L 側に `limiter_ro[ic0]` を、R 側に `limiter_ro[ic1]` を `interp_dispatch` に渡す
(他成分も同様)。Roe カーネルでは追加で Ducros 補正
[`apply_ducros_limiter`](../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L189)
を経由する (SLAU, HLLE 経路は Ducros 補正を直接通さず素のリミッタを使う)。

### 入出力

入力: セル中心保存量・原始量のうち $\rho, U_x, U_y, U_z, P$、対応する勾配
(`drod*, dUxd*, dUyd*, dUzd*, dPd*`)、`cell_planes_index`, `cell_planes`,
`plane_cells`, `pcx/y/z`, `ccx/y/z`, `vol`, `fx`。

出力: `limiter_ro, limiter_Ux, limiter_Uy, limiter_Uz, limiter_P` (各セル 1 値, `[0, 1]`)。

### 既知の TODO / 注意点

- `Ht` 用のリミッタ呼び出しはコメントアウト済み (現状 $H_t$ は再構成しない)。
- Nishikawa R1 は実装ありで未有効。$\delta'_C$ と $r_{ik}$ の計算を別関数 `calcDeltaIJ_dash`
  ([L19-L31](../solver_density_cuda/cuda_forge/limiter_d.cuh#L19)) で行う。
- リミッタ計算は近傍ループ 2 周回 (max/min 取得 → ψ 評価) でレジスタを多く使う。
  `dimGrid_normalcell_small`, `dimBlock_small` のブロック構成は `cudaConfig` で
  個別に調整されている。
- 全境界周期 / 滑らかな計算 (Taylor–Green 等) では `cfg.limiter = 0` (無リミッタ) が有効。

## 関連

- リミッタの参照側 (`interp_dispatch`): [`methods/convection/`](convection/)。
- 勾配の生成: [`methods/gradient/`](gradient/)。

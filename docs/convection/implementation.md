# 対流フラックス — 実装

forge の対流フラックス計算の実装と、ソース上の対応関係をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/convectiveFlux.hpp`](../../solver_density_cuda/convectiveFlux.hpp) | 旧 CPU 経路用ヘッダ (現在ほぼ宣言のみ・コメントアウト) |
| [`solver_density_cuda/convectiveFlux.cpp`](../../solver_density_cuda/convectiveFlux.cpp) | 旧 CPU 実装 (Roe 平均、ヤコビアン構築) |
| [`solver_density_cuda/cuda_forge/convectiveFlux_d.cuh`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cuh) | GPU カーネル / device 関数の宣言 |
| [`solver_density_cuda/cuda_forge/convectiveFlux_d.cu`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) | スキーム本体・ラッパ (約 2700 行) |
| [`solver_density_cuda/cuda_forge/AUSM_d.cu`](../../solver_density_cuda/cuda_forge/AUSM_d.cu) | AUSM 系の補助関数 (M±, β± など) |

実運用は GPU 経路のみ。CPU の `convectiveFlux.cpp` はリファレンス用途で、
通常のビルドでは呼ばれない。

## エントリポイント

[`convectiveFlux_d_wrapper`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L2518)
が共通入口。`cfg.solver` 文字列でディスパッチする。

| `cfg.solver` | 呼び出すカーネル | 状態 |
| --- | --- | --- |
| `"SLAU"` | `SLAU_d` (`slauVariant=1`) | 有効 |
| `"SLAU2"` | `SLAU_d` (`slauVariant=2`) | 有効 (圧力束のみ SLAU2、低マッハ改良) |
| `"HLLE"` | `HLLE_d` | 有効 |
| `"ROE"`  | `ROE_d`  | 有効 |
| `"AUSM+"`, `"AUSM+UP"`, `"KEEP_FVS"`, `"KEEP_SLAU"` | (実装あり) | ラッパでコメントアウト中、`exit(EXIT_FAILURE)` |

呼び出し直前に `res_ro`, `res_roUx`, `res_roUy`, `res_roUz`, `res_roe` を
`cudaMemset` でゼロクリアする。

## 状態再構成と dispatch

すべてのスキームが共通に [`interp_dispatch`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L142)
を呼び、`cfg.convMethod` (= `scheme`) と `cfg.limiter` (= `limit_scheme`) に応じて
次の device 関数を選ぶ。

| `scheme` | 関数 | 内容 |
| --- | --- | --- |
| `0`, `-1` | `interp_1stUp` | 1 次風上 (`phif = phiC`) |
| `1` | `interp_MUSCL_2nd` | 線形再構成: $\phi_C + \psi\,\nabla\phi_C \cdot \mathbf{cp}$ |
| `2` | `interp_MUSCL_3rd` | $k = 1/3$ で 3 次 MUSCL |
| (上記以外) | `interp_MINMOD` | 上流値 $\phi_U$ を構成して minmod $r$ を計算し制限 |

L 側は `(phiC, phiD) = (ro[ic0], ro[ic1])` ほか、距離ベクトル `dcc`, `dc0p` を渡し、
R 側は左右を反転 (`-dcc`, `dc1p`, `1-f`) して呼ぶ。
`limiter_*[ic]` は別途生成済みのスカラリミッタを各成分独立に渡す。

### Ducros センサの適用

ROE_d など一部のカーネルで [`apply_ducros_limiter`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L189) を経由して

```cpp
duc = clamp(max(ducros[ic0], ducros[ic1]), 0.0, 1.0);
limiter = (duc > 0.8) ? max(0.0, (1.0 - duc) * limiter) : limiter;
```

を適用する (Ducros 値が高い、すなわち圧縮性ショック近傍で再構成を抑制)。

## 各カーネルの構造

共通する処理パターン:

1. `ip_orig = blockDim*blockIdx + threadIdx` から `normal_halo_planes_d[ip_orig]` で
   実際の面 ID `ip` を引く (有効面のみで並列化するための間接参照)。
2. `plane_cells[2*ip+0/1]` で両側セル `ic0`, `ic1` を得る。
3. `ic1 >= nCells` (= 非周期境界のゴースト) のとき `conv_scheme = -1` に強制し
   1 次風上に落とす。
4. `interp_dispatch` で L/R 状態を作る。`ro_L = max(ro_L, small_rho)` 等で
   負値 / 0 をクランプ (Roe のみ)。
5. スキーム固有のフラックス計算 (下記)。
6. `atomicAdd(&res_*[ic0], -res_*_temp)` / `atomicAdd(&res_*[ic1], +res_*_temp)`
   で両側に符号反転加算。

### `SLAU_d` ([L218](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L218))

並列単位は `nNormal_ghst_Planes` (`normal_ghst_planes_d[]` で実面 ID へ間接参照)。

1. `interp_dispatch` で `(ro, Ux, Uy, Uz, Ps)` の L/R 状態を構築。
   `velocity2_L/R = |u|²`, `h_p/m = γP/((γ-1)ρ) + 0.5|u|² = H_t` を計算。
2. 法線速度 `Vn_p = u_L·n̂`, `Vn_m = u_R·n̂` (`/sss` で正規化)。
3. **音速平均** `c_hat = 0.5*(sonic[ic0] + sonic[ic1])` (算術平均、Roe 平均ではない)。
4. **Mach 関数** `beta_p, beta_m` を分岐式でその場計算 (`|M|≥1` と `<1` で別式)。
   AUSM ファミリの汎用ヘルパ ([`AUSM_d.cu`](../../solver_density_cuda/cuda_forge/AUSM_d.cu)) は
   呼ばず、SLAU 用に直接ベタ書き。
5. **低マッハ補正係数**:
   ```cpp
   g            = -max(min(M_p, 0), -1) * min(max(M_m, 0), 1);
   Vn_hat_abs   = (ro_L*|Vn_p| + ro_R*|Vn_m|) / (ro_L + ro_R);
   Vn_hat_p_abs = (1-g)*Vn_hat_abs + g*|Vn_p|;
   Vn_hat_m_abs = (1-g)*Vn_hat_abs + g*|Vn_m|;
   M_hat        = min(1, sqrt(0.5*(|u_L|² + |u_R|²)) / c_hat);
   chi          = (1 - M_hat)²;
   ```
6. **圧力束** `p_tilde` と **質量流束** `mdot`。圧力束の**第 3 項のみ** `slauVariant` で
   分岐する (`mdot` は両版共通):
   ```cpp
   // 第 3 項: slauVariant==2 (SLAU2) のときだけ差し替え
   if (slauVariant == 2)               // SLAU2 (低マッハ圧力散逸)
       p_third = (beta_p+beta_m-1)*sqrt(0.5*(velocity2_L+velocity2_R))
                 *0.5*(ro_L+ro_R)*c_hat;
   else                                // SLAU
       p_third = (1-chi)*(beta_p+beta_m-1)*0.5*(P_L+P_R);
   p_tilde = 0.5*(P_L+P_R) + 0.5*(beta_p-beta_m)*(P_L-P_R) + p_third;
   mdot    = sss*0.5*((ro_L*(Vn_p+Vn_hat_p_abs) + ro_R*(Vn_m-Vn_hat_m_abs))
                      - chi/c_hat * (P_R-P_L));
   ```
   SLAU2 の根拠と式は [`theory.md` の SLAU2 節](theory.md#slau2-圧力束の低マッハ改良) を参照。
   `velocity2_L/R`, `c_hat` は既存量、新規は `0.5*(ro_L+ro_R)` のみ。
7. **残差**: `0.5*(mdot ± |mdot|)` で風上選択し、運動量に `p_tilde * S_*` を加え、
   エネルギは `h_p, h_m` (全エンタルピ) を風上にとる。
8. 両側セルに `atomicAdd(±)`。

リミッタは Ducros 補正を通さず生の `limiter_*[ic0/ic1]` をそのまま渡している
(Roe との違い)。

### `ROE_d` ([L1474](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L1474))

並列単位は `nNormal_halo_Planes`。レジスタ圧が最も高いカーネル。

1. **Ducros 補正リミッタ**: `apply_ducros_limiter(limiter_*[ic0/1], duc)` で
   L/R 別の補正済みリミッタを作る (`duc = clamp(max(ducros[ic0], ducros[ic1]), 0, 1)`)。
2. `interp_dispatch` で `(ro, Ux, Uy, Uz, Ps)` の L/R 状態を構築。
   `ro_*, P_* = max(*, 1e-8)` でクランプ。`roe_*, Ht_*, ca_*` を再計算。
3. **Roe 平均状態** (`sqrt_ro_L/R` から `roa, ua, va, wa, Ha, ca, Ua`)。
4. **保存量差分** `delQ[]` の並びは `(Δρ, Δ(ρE), Δ(ρu), Δ(ρv), Δ(ρw))` で
   一般的教科書順 `(ρ, ρu, ρv, ρw, ρE)` とは異なる:
   ```cpp
   delQ[0] = ro_R - ro_L;                            // Δρ
   delQ[1] = (ro_R*Ht_R - P_R) - (ro_L*Ht_L - P_L);  // Δ(ρE) = Δ(ρe)
   delQ[2] = ro_R*Ux_R - ro_L*Ux_L;                  // Δ(ρu)
   delQ[3] = ro_R*Uy_R - ro_L*Uy_L;
   delQ[4] = ro_R*Uz_R - ro_L*Uz_L;
   ```
5. **圧力の保存変数微分**:
   ```cpp
   P_ro  = 0.5*(γ-1)*(u²+v²+w²);  P_roe = γ-1;
   P_rou = -(γ-1)*u;              P_rov = -(γ-1)*v;   P_row = -(γ-1)*w;
   z   = u²+v²+w² - P_ro/P_roe;   // = 0.5*|u|² (理想気体)
   fai = P_ro - c²;
   ```
6. **右固有ベクトル `R[5][5]`** をその場で構築 (列の順は遅い音波 / せん断 ×3 / 速い音波)。
7. **左固有ベクトル `L[5][5]`** をベタ書きし、最後に二重ループで一括 `L /= c²`
   (実装上のショートカット)。
8. **特性座標** `dW = L * delQ`、**固有値** `lam[] = (|U-c|, |U|, |U|, |U|, |U+c|)`。
9. **Harten 形エントロピーフィックス** (速度比依存型のみ `if (true)`、他はガード `if (false)`):
   ```cpp
   eta_vl = 0.1*(|Ua|/ca + 1.0);
   if (lam[i] <= eta_vl)
       lam[i] = (lam[i]² + eta_vl²) / (2*eta_vl);
   ```
10. **散逸ベクトル** `difQ1 = |lam| * dW`、`difQ2 = R * difQ1`。
11. **中央束** を質量流束 `mdot = 0.5*(ρU_L + ρU_R)*sss` 経由で書き、
    `0.5*(F_L + F_R) - 0.5*difQ2*sss` の形で残差に積算。
12. 両側セルに `atomicAdd(±)`。

### `HLLE_d` ([L1268](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L1268))

並列単位は `nNormal_halo_Planes`。Roe より軽量で堅牢。

1. `interp_dispatch` で L/R 構築 → `max(*, 1e-8)` クランプ。
2. `roe_L/R = P/(γ-1) + 0.5*ρ|u|²`, `Ht_L/R = (roe + P)/ρ`, `U_L/R = u·n̂`,
   `ca_L/R = sqrt(γP/ρ)`。
3. **Roe 平均**で `ua, va, wa, Ha, ca, Ua` を作り、**波速見積もり**:
   ```cpp
   S_L = min(U_L - ca_L, min(U_R - ca_R, Ua - ca));
   S_R = max(U_L + ca_L, max(U_R + ca_R, Ua + ca));
   ```
4. **物理フラックス** `F_L_*, F_R_*` を $(\rho U, \rho u U + P n_x, \dots, \rho H_t U)$ で構築。
5. **三分岐**:
   - `S_L ≥ 0` → `F_L * sss`
   - `S_R ≤ 0` → `F_R * sss`
   - 中間 → `(S_R*F_L - S_L*F_R + S_L*S_R*(Q_R - Q_L)) / max(S_R-S_L, 1e-8) * sss`

   中間状態のエネルギ差は `roe_R - roe_L` (= $\rho e_R - \rho e_L$) を使う。
6. 両側セルに `atomicAdd(±)`。

Ducros 補正は通さず、`limiter_*[ic0/ic1]` を直接 `interp_dispatch` に渡す。

### KEEP / AUSM 系

[`KEEP_FVS_d`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L864),
[`KEEP_SLAU_d`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L1827),
[`AUSMp_d`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L461),
[`AUSMp_UP_d`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L646)
は実装されているがラッパで現在無効化されている。`境界処理の見直し中` という
コメントが残っているので、復活時はゴーストセルとの整合 (`ic1 >= nCells` 分岐)
を再確認する。

## 並列化メモ

- 並列単位は面 1。`normal_halo_planes_d` で実際に処理する面のみ列挙し、
  `dimGrid_normal_halo = ceil(nNormal_halo_Planes / blocksize)` をグリッドサイズとする。
- 面ごとに両側セルへ `atomicAdd`。`flow_float = double` の場合 CC 6.0+ が必要。
- 面ループ内で `R[5][5], L[5][5]` などローカル配列を構築するためレジスタ圧が高い。
  Roe カーネルは特に重い。

## 入出力

入力: 保存量 / 原始量 (`ro, roUx, roUy, roUz, roe, Ux, Uy, Uz, P, Ht, sonic`)、
勾配 (`drod*`, `dUxd*`, `dUyd*`, `dUzd*`, `dPd*`)、リミッタ
(`limiter_ro, limiter_Ux, limiter_Uy, limiter_Uz, limiter_P`)、Ducros 値 (`ducros`)、
セル中心・面中心座標、面ベクトル `(sx, sy, sz, ss)`、補間重み `fx`、`plane_cells`。

出力 (累積): セル中心残差 `res_ro, res_roUx, res_roUy, res_roUz, res_roe`。

## 既知の TODO / 注意点

- KEEP / AUSM 経路が無効化中 (境界とのカップリング再整理待ち)。
- `c_hat` の中央平均は単純算術平均。Roe 平均の方が低マッハ性に望ましい場面あり。
- リミッタは `interp_dispatch` 経由で各成分 (`ro, Ux, Uy, Uz, P`) 独立。
  限界保証 (TVD) は厳密には未保証 (MUSCL 2/3 次経路は無リミット)。
- 旧 CPU 経路 (`convectiveFlux.cpp`) は使われていない。GPU と同期させる場合は
  ラッパと dispatch 構造を反映する必要がある。

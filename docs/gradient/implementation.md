# 勾配再構成 — 実装

forge の Green-Gauss 勾配計算の実装と、ソース上の対応関係をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/gradient.hpp`](../../solver_density_cuda/gradient.hpp) | CPU エントリ `gradientGauss` の宣言 |
| [`solver_density_cuda/gradient.cpp`](../../solver_density_cuda/gradient.cpp) | CPU 実装 + GPU ラッパへの委譲 |
| [`solver_density_cuda/cuda_forge/calcGradient_d.cuh`](../../solver_density_cuda/cuda_forge/calcGradient_d.cuh) | GPU カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/calcGradient_d.cu`](../../solver_density_cuda/cuda_forge/calcGradient_d.cu) | GPU カーネル実装 |

## エントリポイントとディスパッチ

`gradientGauss(cfg, cuda_cfg, msh, v)` が CPU/GPU 共通の入口。
`cfg.gpu == 1` の場合は GPU ラッパ `calcGradient_d_wrapper` に委譲し、
CPU 経路では同関数内で全処理を行う。実運用 (`main.cpp`) では GPU ラッパが
直接呼ばれている (CPU 経路は CPU デバッグ用)。

呼び出し箇所:

- 初期化直後の勾配計算 ([`main.cpp` L624](../../solver_density_cuda/main.cpp#L624))
- 時間ステップ内の各サブステップ ([`main.cpp` L697](../../solver_density_cuda/main.cpp#L697), [L727](../../solver_density_cuda/main.cpp#L727))

## CPU 実装の構造

[`gradient.cpp`](../../solver_density_cuda/gradient.cpp) の処理は次の 3 段。

1. **ゼロクリア**: 全セル (`nCells_all`、ゴースト含む) について
   $\nabla U_{x,y,z}, \nabla P, \nabla T$ を 0 に初期化。
2. **面ループ加算**: `for ip in [0, nPlanes)`
   - 両側セル `ic1 = iCells[0]`, `ic2 = iCells[1]` と面ベクトル `surfVect` を取得。
   - 重み `f = fxp[ip]` で線形補間: `Uxf = f*Ux[ic1] + (1-f)*Ux[ic2]` ほか。
   - `dXdx[ic1] += sv[0]*Xf`、`dXdx[ic2] -= sv[0]*Xf` (符号反転で対称加算)。
3. **体積除算**: `for ic in [0, nCells)` でセル体積 `msh.cells[ic].volume` で除する。

CPU 経路では密度勾配 $\nabla \rho$ は計算していない (上位で必要としていないため)。

## GPU 実装の構造

[`calcGradient_d_wrapper`](../../solver_density_cuda/cuda_forge/calcGradient_d.cu#L384) が次のカーネルを順に発行する。

1. **`cudaMemset` 群**: $\nabla U_*, \nabla \rho, \nabla P, \nabla T,
   \nabla \cdot \mathbf{u}$ をデバイスバッファ上で 0 クリア。
2. **`calcGradient_1_d`**: 面 1 並列。`atomicAdd` で
   両側セルの勾配ベクトル成分と $\nabla \cdot \mathbf{u}$ に寄与を累積。
   `plane_cells[2*ip+0/1]` でセル ID を取得し、`fx[ip]` で線形補間。
3. **`calcGradient_2_d`**: セル 1 並列。`vol[ic]` で除算して勾配を確定。

CPU 版と異なり、密度 $\rho$ の勾配 $\nabla \rho$ も同時に計算する。
$\rho U_*, \rho e, H_t$ の勾配計算用コードは保留 (コメントアウト済み)。

### 境界寄与

`calcGradient_b_d` 境界カーネルは [`calcGradient_d.cu`](../../solver_density_cuda/cuda_forge/calcGradient_d.cu) に
定義されているが、現在ラッパ内ではコメントアウトされている。
境界面の寄与はゴーストセルを内部面ループの `ic1` として処理する設計
(ゴーストセル値は境界条件モジュールが事前に書き込む)。
将来的に境界に固有の値 (壁面圧力など) を直接使う場合に
このカーネルを有効化する。

## 並列化メモ

- 面 1 並列の累積は `atomicAdd` を使用。`flow_float` が `double` の場合は
  Compute Capability 6.0 以降が必要。
- `__syncthreads()` 呼び出しは入退出の安全策 (ループ内同期はない)。
- `dimGrid_plane`, `dimGrid_cell` は `cudaConfig` で `nPlanes`, `nCells` から計算。

## 入出力データ

入力 (セル中心配列): `ro`, `Ux`, `Uy`, `Uz`, `P`, `T`, `roe`, `Ht`, セル体積、
面ベクトル (`sx, sy, sz, ss`)、補間重み `fx`、`plane_cells` 接続、
GPU 側のミラー `var.c_d[*]`, `var.p_d[*]`。

出力 (セル中心配列):

```
dUxdx, dUxdy, dUxdz
dUydx, dUydy, dUydz
dUzdx, dUzdy, dUzdz
dPdx,  dPdy,  dPdz
dTdx,  dTdy,  dTdz
drodx, drody, drodz   (GPU only)
divU
```

## 既知の TODO / 注意点

- limiter / 高次再構成は未実装。コード内では `gradientGauss` の単一経路のみ。
- $\rho U_*, \rho e, H_t$ の勾配計算用変数とカーネル分岐は `cu` 内に
  コメントとして残されている。必要時に有効化する。
- 境界カーネル `calcGradient_b_d` は無効化中。境界条件設計と合わせて再評価する。

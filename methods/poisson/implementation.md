# Poisson 求解 — 実装

forge の Poisson 求解の実装とソース対応をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/solvePoisson_amgcl.hpp`](../../solver_density_cuda/solvePoisson_amgcl.hpp) | AMGCL バックエンド宣言 |
| [`solver_density_cuda/solvePoisson_amgcl.cpp`](../../solver_density_cuda/solvePoisson_amgcl.cpp) | 行列構築・求解・補正 (約 610 行) |
| [`solver_density_cuda/amgcl_cuda/solvePoisson_amgcl_cuda.{h,cpp}`](../../solver_density_cuda/amgcl_cuda/) | AMGCL CUDA バックエンドラッパ |
| [`solver_density_cuda/solvePoisson_amgx.hpp`](../../solver_density_cuda/solvePoisson_amgx.hpp) | AMGX バックエンド宣言 |
| [`solver_density_cuda/solvePoisson_amgx.cpp`](../../solver_density_cuda/solvePoisson_amgx.cpp) | AMGX 版実装 (約 570 行) |
| [`solver_density_cuda/calcDistance_poisson.cpp`](../../solver_density_cuda/calcDistance_poisson.cpp) | 壁面距離 Poisson の独立プログラム (cuSPARSE CG) |

圧縮性メイン (`main.cpp`) では `solvePoisson` の呼び出しは現状コメントアウト
されており、SMAC 系のケースで明示的に組み込む構成。

## エントリポイント

### AMGCL 版

| 関数 | 役割 |
| --- | --- |
| [`setMatrixPoisson`](../../solver_density_cuda/solvePoisson_amgcl.cpp#L26) | 内部面・境界面を回り `mat_p.lhs[ic][*]`, `mat_p.rhs[ic]` を構築 |
| [`solvePoisson`](../../solver_density_cuda/solvePoisson_amgcl.cpp#L124) | AMGCL の `make_solver` で AMG + Krylov を解く |
| [`correctPresVel`](../../solver_density_cuda/solvePoisson_amgcl.cpp#L194) | 解 $p'$ から圧力・速度補正 |
| [`correctPresVel_BF`](../../solver_density_cuda/solvePoisson_amgcl.cpp#L382) | バックワード差分版補正 |
| [`callAmgclCuda`](../../solver_density_cuda/solvePoisson_amgcl.cpp#L570) | CUDA バックエンド用エントリ |

### AMGX 版

同じインタフェース (`setMatrixPoisson`, `solvePoisson`, `correctPresVel`,
`correctPresVel_BF`) を `solvePoisson_amgx.cpp` で提供。ビルド時にいずれかをリンクする。

### 壁面距離

[`calcDistance_poisson.cpp`](../../solver_density_cuda/calcDistance_poisson.cpp) は
独立した `main` を持つツール (`make_laplace_matrix` + `gpu_CG`)。
cuSPARSE/cuBLAS の標準サンプルに準拠した CG 実装で、外部から呼ばれるユーティリティ。

## 行列構築の概要 (`setMatrixPoisson`)

```cpp
for ip in [0, nNormalPlanes):
    ic0 = planes[ip].iCells[0];
    ic1 = planes[ip].iCells[1];
    dcc = |cells[ic1].cent - cells[ic0].cent|;
    ss  = planes[ip].surfArea;

    mat_p.lhs[ic0][0]      -= ss/dcc;
    mat_p.lhs[ic0][ip_loc0] = +ss/dcc;
    mat_p.lhs[ic1][0]      -= ss/dcc;
    mat_p.lhs[ic1][ip_loc1] = +ss/dcc;

for bc in msh.bconds:
    if bc.valueTypes["P"] == 1:  // ディリクレ境界
        for ip in bc.iPlanes:
            dn = ((pcent - c1cent) ⋅ sv) / ss;
            mat_p.lhs[ic0][0] -= ss/dn;

rhs[ic] = divU_vol[ic] * ro[ic] / cfg.dt;

if (全境界 Neumann):
    mat_p.lhs[0][0] = 1.0e+30;
    mat_p.rhs[0]    = 0.0;
```

`localPlnOfCell[ip][0/1]` でセル内ローカルインデックスを管理し、対称疎行列を構成する。

## 既知の TODO / 注意点

- 圧縮性メイン経路では Poisson は使われない (時間積分は密度ベース陽 / 陰)。
  SMAC 用は別経路として残っている。
- AMGCL CUDA バックエンド (`amgcl_cuda/`) と AMGX のどちらを使うかはビルド構成依存。
  選択肢は [`CMakeLists.txt`](../../solver_density_cuda/CMakeLists.txt) のオプションを参照。
- `calcDistance_poisson.cpp` は独立 `main` を持つため、ライブラリへの統合は今後の課題。
- 行列構築は CPU 経路。GPU で構築するルートは現状無く、ホスト側で組んでデバイスに転送する。

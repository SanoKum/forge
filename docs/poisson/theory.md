# Poisson 求解 — 理論

forge には次の 2 種類の Poisson 求解が含まれる。

1. **圧縮性ソルバの壁面距離計算** — RANS / LES の壁面距離 $d_w$ を Eikonal 方程式の
   緩和版として Poisson で求める (圧縮性メイン経路では現状未使用)。
2. **非圧縮 (SMAC) 経路の圧力 Poisson** — 旧 CPU 実装 (`solvePoisson_amgcl.cpp` 等)。
   圧縮性 (`solver_density_cuda`) のメイン時間積分では使用されないが、
   コードベースには残っており、SMAC ベースのケースに利用できる。

実装の対応は [implementation.md](implementation.md) を参照。

## 圧力 Poisson 方程式

SMAC 系では予測速度から $\nabla \cdot \mathbf{u}^*$ を計算し、補正圧力 $p'$ を

$$
\nabla \!\cdot\! \left(\frac{1}{\rho}\nabla p'\right)
= \frac{\nabla \cdot (\rho \mathbf{u}^*)}{\Delta t}
$$

として求める。これを有限体積離散化すると面ごとに

$$
\sum_{f \in \partial c} \frac{A_f}{|\mathbf{d}|_f} (p'_{c_1} - p'_{c_0})
= \frac{V_c}{\Delta t}\, (\nabla \cdot \mathbf u^*)_c \cdot \rho_c
$$

の対称 5/7 点ステンシルとなり、AMGCL/AMGX で解く。

## 壁面距離 (Eikonal 緩和)

Spalding 系の壁面距離 $d_w$ は $|\nabla d_w| = 1$ を満たすが、これは Poisson の緩和

$$
\nabla^2 \phi = -1 ,\quad d_w = -|\nabla\phi| + \sqrt{|\nabla\phi|^2 + 2\phi}
$$

から近似的に得られる。境界条件は固体壁 $\phi = 0$、その他はノイマン零。
このアプローチは大規模並列で扱いやすく、forge の壁面距離計算で採用されている。

## 境界条件

- 壁面: ディリクレ $p' = 0$ または $\phi = 0$。
- 流入/流出: 物理に応じてディリクレ (固定圧) / ノイマン (圧力勾配ゼロ)。
- 全境界が Neumann のときは行列が特異になるため、最初の行をピン留めする
  (`lhs[0][0] = 1e+30`, `rhs[0] = 0`)。

## 線形ソルバ

forge では次のバックエンドを選択できる。

| バックエンド | 役割 |
| --- | --- |
| **AMGCL** | CPU/CUDA 両対応。共役勾配 (CG / BiCGStab) + AMG 前処理 |
| **AMGX** | NVIDIA の AMG ライブラリ。GPU 専用 |
| **cuSPARSE/cuBLAS CG** | 壁面距離計算で利用 (NVIDIA サンプル準拠) |

## 参考

- 圧縮性メイン経路の時間積分は陽 / 陰の defect-correction で完結し、
  Poisson は不要。詳細は [`docs/time_integration/`](../time_integration/) 参照。
- メッシュ・境界 ID の与え方は [`docs/boundary/`](../boundary/) を参照。

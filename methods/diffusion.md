# 粘性・熱伝導フラックス

forge は Navier–Stokes の粘性応力と熱伝導束を、セル中心勾配の面平均と
セル中心値の法線差分を組み合わせた "over-relaxed" スキームで離散化する。
実装の対応は 本ドキュメントの「実装」節 を参照。

本ドキュメントは理論(係数・方程式)と実装(ソース対応)をまとめる。

## 理論

### 連続方程式系

粘性応力テンソル (Newton 流体、Stokes 仮定) と熱伝導束:

$$
\tau_{ij} = \mu \!\left( \frac{\partial u_i}{\partial x_j}
+ \frac{\partial u_j}{\partial x_i} \right)
- \frac{2}{3}\mu (\nabla \cdot \mathbf{u})\,\delta_{ij},
\qquad
q_i = -\kappa\, \frac{\partial T}{\partial x_i}.
$$

セル境界フラックスは

$$
\mathbf{F}^{\text{visc}} =
\begin{pmatrix}
0\\
\tau_{xj}\,n_j\\
\tau_{yj}\,n_j\\
\tau_{zj}\,n_j\\
\tau_{ij} u_i n_j + \kappa\, (\partial_j T)\, n_j
\end{pmatrix}.
$$

### Over-relaxed 離散化

セル $c_0, c_1$ をまたぐ面 $f$ で速度・温度の勾配を求めるとき、
セル中心間ベクトル $\mathbf{d} = \mathbf{r}_{c_1} - \mathbf{r}_{c_0}$ と面ベクトル
$\mathbf{S}$ の非直交を補正する。$|\mathbf{S}| = S$ として

$$
\boldsymbol{\Delta} = \frac{|\mathbf{S}|^2}{\mathbf{d}\cdot\mathbf{S}}\,\mathbf{d},
\qquad \mathbf{k} = \mathbf{S} - \boldsymbol{\Delta}.
$$

法線方向成分はセル中心差分から、接線成分は両側セルの勾配を面平均
($f_x$ 重み) して構成する:

$$
(\partial_j u_i)\, S_j \;\approx\;
\underbrace{\frac{u_i^{c_1} - u_i^{c_0}}{|\mathbf{d}|} \, \Delta_j}_{\text{法線差分}}
\;+\;
\underbrace{\overline{(\partial_j u_i)}_f \, k_j}_{\text{接線補正}},
$$

$\overline{(\cdot)}_f = f_x (\cdot)_{c_0} + (1-f_x)(\cdot)_{c_1}$。

この組み合わせを直接 $\tau$ に代入することで、粘性束を面ごとに評価し、
両側セルに $\pm$ で加算する。エネルギ束も同様に $\tau_{ij} u_i$ と
$\kappa (\partial_j T)$ を組み合わせる。

> **重要 — 法線項はスカラー係数で評価する。** 上式の法線差分項
> $(u_i^{c_1}-u_i^{c_0})/|\mathbf{d}|\cdot\Delta_j$ で各成分 $i$ に掛かるのは、
> over-relaxed の**スカラー** $\beta = |\mathbf{S}|^2/(\mathbf{d}\cdot\mathbf{S})$ である。
> すなわち $\mu\,\nabla u_i\cdot\mathbf{S}$ の法線寄与は
> $\mu\,\beta\,(u_i^{c_1}-u_i^{c_0})$(全成分共通のスカラー $\beta$)であり、
> 熱伝導束 $\kappa\,\beta\,(T^{c_1}-T^{c_0})$ と同型である。
> 速度成分ごとに $\boldsymbol\Delta$ の**自身の座標成分** $\Delta_i$ を掛けると
> (例: $\tau_x$ に $\Delta_x$)、軸平行な $y$ 法線面で $\Delta_x=0$ となり
> 流れ方向運動量の横方向拡散 $\mu\,\partial u_x/\partial y$ が落ちる。
> Laplacian 接線補正 $\overline{\nabla u_i}_f\cdot\mathbf{k}$ は**同一成分** $u_i$ の勾配
> $\overline{\partial u_i/\partial x_j}\,k_j$ を用いる。転置寄与 $\mu\,\partial u_j/\partial x_i\,S_j$ は
> 別項として面平均勾配にフル $\mathbf{S}$ を内積して加える。実装は
> 本ドキュメントの「実装」節 を参照。

### 粘性係数

- 層流粘性 $\mu_{\text{lam}}$: 設定値 (`cfg.visc`、定数) を使用。
- 乱流粘性 $\mu_{\text{turb}}$: `vis_turb[ic]` (SGS / RANS モデル別ファイルで生成)。
- 有効粘性 $\mu = \mu_{\text{lam}} + \mu_{\text{turb}}$。
- 熱伝導率 $\kappa$: **有効値 $\kappa_{\text{eff}} = \kappa_{\text{lam}} + c_p\,\mu_{\text{turb}}/Pr_t$**。
  $Pr_t$ は `turbulence.turbulentPrandtl` で設定 (既定 0.85)。層流 $\kappa_{\text{lam}}$ は設定値
  `cfg.thermCond`、乱流寄与は $\mu_{\text{turb}}$ から渦熱伝導として加える。$c_p$ は
  thermally-perfect の温度依存を反映するためセル配列 (`var.c_d["cp"]`、`thermalMethod==2` で $c_p(T)$)
  を面平均して使う。応力が $\mu=\mu_{\text{lam}}+\mu_{\text{turb}}$ を使う以上、熱伝導も対応する乱流寄与を
  含めないと乱流境界層でエネルギーが保存せず、断熱壁の静温が回復温度 ($\le$ 全温) を超えて発散的に上昇する。

### 壁面寄与

非滑り条件下の壁面では速度がゼロとなる。壁面用カーネル `viscousFlux_wall_d` で
壁面距離 (ghost cell までの法線距離) を使って粘性応力を片側差分で評価する。
これにより内部面ループでは未処理の壁面寄与を補う。

壁面でも粘性束は内部面と**同じ** $\tau_{ij} n_j$ の評価でなければならない。
すなわち運動量 $i$ 成分に対し

$$
F^{\text{visc}}_i = \tau_{ij}\, S_j
= \mu\!\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) S_j
- \frac{2}{3}\mu (\nabla\cdot\mathbf{u})\, S_i .
$$

平板チャネル等の軸平行壁 (法線 $\mathbf{n}=\pm\mathbf{e}_y$) では、壁摩擦の主項は
$\tau_{xy} S_y = \mu\,(\partial u_x/\partial y)\, S_y$ という**せん断 (off-diagonal) 成分**であり、
これがストリーム方向運動量 $\rho u_x$ への no-slip 抗力を与える。
したがって法線差分 $(u_i^{g}-u_i^{c})/|\mathbf{d}|$ は内部面と同じく
over-relaxed 係数 $\Delta_i$ (直交格子では $\Delta_i = S_i$、すなわち面積 $S$ に法線
単位ベクトルを掛けたもの) と組み合わせ、$\tau_{ij}S_j$ として評価する必要がある。
速度成分ごとに自身の座標成分 $S_i$ のみを掛ける ($\tau_x \propto S_x$ 等) と、
軸平行壁では主せん断項が落ちて壁摩擦がゼロになる (この不具合は 2026-06-06 に修正済み。
実装は 本ドキュメントの「実装」節 を参照)。

### 参考

- 対流フラックスとの加算: [`methods/convection/`](convection) も参照。
  両者の残差は同じ `res_*` バッファに累積される。
- 勾配の生成: [`methods/gradient/`](gradient) を参照。

## 実装

### ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/cuda_forge/viscousFlux_d.cuh`](../solver_density_cuda/cuda_forge/viscousFlux_d.cuh) | カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/viscousFlux_d.cu`](../solver_density_cuda/cuda_forge/viscousFlux_d.cu) | 内部面・壁面カーネル本体 (約 420 行) |

CPU 経路の独立実装は無く、GPU 経路のみで運用される。

### エントリポイント

[`viscousFlux_d_wrapper`](../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L315) が共通入口。
順に次を呼ぶ。

1. `viscousFlux_d` — 内部面 (`nNormalPlanes`) を並列処理。
2. `viscousFlux_wall_d` — 各壁面境界条件について壁面寄与を計算。

`res_*` バッファは対流フラックス計算後に粘性が追加加算される (ゼロクリアは対流側で実施)。

### `viscousFlux_d` 構造 ([L3](../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L3))

面 1 並列のカーネル。`ip < nNormalPlanes` のみ処理 (境界面は別カーネル)。

1. 両側セル `ic0, ic1` を取得。
2. セル中心間ベクトル `dcc` と長さ `|dcc|`。
3. 面側で速度 `Ux/y/zf` を線形補間 (`f = fx[ip]`)。
4. 各成分勾配を面平均 (`dUxdxf` ほか) し、温度勾配 `dT*f` も同様。
5. **Over-relaxed 補正係数**を構築:
   ```cpp
   delta   = dcc*sss*sss/(dcc_x*sxx + dcc_y*syy + dcc_z*szz);
   delta_x = dcc_x*sss*sss/(dcc_x*sxx + dcc_y*syy + dcc_z*szz);
   k_x     = sxx - delta_x;   // (k_y, k_z も同様)
   ```
6. 粘性係数を面平均: `mu_total = (vis_lam + vis_turb) を f で補間`。
7. 法線差分 `(U1-U0)/dcc * delta_*` と接線補正 `mu*(dU*df*k_x + ...)`、
   発散項 `-mu*(2/3)*divU*S_*` を合算して `tau_x, tau_y, tau_z` を得る。
8. 熱伝導束 `heatflux = tc_face*(T1-T0)/dcc * delta + tc_face*(dT*df * k)`。
   有効熱伝導率は **`tc_face = thermCond(層流) + cp_face*vis_turb/Pr_t`** (RANS 乱流熱伝導)。
   `Pr_t` は `turbulence.turbulentPrandtl` 設定 (既定 0.85)、`cp_face` はセル `cp` 配列の面平均
   (thermally-perfect の $c_p(T)$ 反映)。
   応力(摩擦発熱)は `mu_total = vis_lam + vis_turb` を使うので、熱伝導も同じ乱流寄与を
   含めないとエネルギーが保存せず、乱流境界層で散逸熱が逃げ場を失い静温が全温を超えて
   overshoot する (2026-06 修正。詳細 [`.github/plans/diffusion-turbulent-thermal-conductivity.md`](../plans/accepted/diffusion-turbulent-thermal-conductivity.md))。
9. 残差 `res_roUx, res_roUy, res_roUz, res_roe` を両側に符号反転で `atomicAdd`。

エネルギ残差には `tau_x*Uxf + tau_y*Uyf + tau_z*Uzf` (応力仕事) と `heatflux` を加える。

各運動量成分 $i$ の応力 `tau_i` は完全な Newton 応力 $\tau_{ij}S_j$ を3項で構成する
([L105-123](../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L105)、熱伝導束と同型):

```cpp
// 例: tau_x (i=x)
tau_x  = mu*((Ux1-Ux0)/dcc)*delta;                  // (1) Laplacian 法線 (スカラー delta)
tau_x += mu*(dUxdxf*k_x + dUxdyf*k_y + dUxdzf*k_z); // (2) Laplacian 接線補正 (同成分 ∂u_x/∂x_j)
tau_x += mu*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz); // (3) 転置 ∂u_j/∂x にフル S
tau_x += -mu*2.0/3.0*divu*sxx;                      // (4) 発散項 (成分 S_x)
```

(1)+(2) が $\mu\,\nabla u_i\!\cdot\!\mathbf{S}$ の over-relaxed 評価、(3) が転置 $(\nabla u)^T$ 寄与、
(4) が Stokes 仮定の体積粘性。`tau_y`/`tau_z` は成分を入れ替えて同型。

> **軸対称の発散項 (2026-06-14 修正)** — (4) の体積粘性で使う発散は、軸対称では
> $\nabla\!\cdot\!\mathbf u = \partial_x u_x + \partial_r u_r + u_r/r$ と**フープ項 $u_r/r$ を含む**形になる。
> 以前は planar 面の `divu` をデカルト形 `dUxdxf+dUydyf+dUzdzf` のみで評価しており、
> フープ応力 $\tau_{\theta\theta}=2\mu\,u_r/r-\tfrac23\mu(\nabla\!\cdot\!\mathbf u)$ 
> ([axisymmetricSource_d.cu](../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)、完全発散 `axisym_divU` を使用)
> と**不整合**だった ($\tau_{xx},\tau_{rr}$ の体積項だけ $u_r/r$ を落としていた)。
> 修正: `isAxisymmetric==1` のとき planar 面は `axisym_divU` を面補間
> (内部面 `f*axisym_divU[ic0]+(1-f)*axisym_divU[ic1]`、壁面はセル値 `axisym_divU[ic]`) で
> `divu` を取り、$\tau_{\theta\theta}$ と同一の完全発散に揃える。`axisym_divU` は `axisymmetricGeomTerms_d`
> が viscousFlux より前に算出済み (main ループ順)。非軸対称は従来どおりデカルト発散で**ビット不変**。
> 計画は [`.github/plans/diffusion-viscous-shear-flux.md`](../plans/accepted/diffusion-viscous-shear-flux.md) §変更ログ。

> **履歴 (2026-06-06 修正)** — 以前は法線項に**成分** `delta_x`(`=dcc_x·β`)を使い、
> 接線項の勾配添字が転置になっており、(3) の転置項もコメントアウトされていた。このため
> 軸平行な $y$ 法線面で `delta_x=0` → 流れ方向運動量の横方向拡散 $\mu\,\partial u_x/\partial y$ が
> 落ち、後述の壁面 `*sxx` 不具合と相まって境界層が形成されなかった。
> 計画は [`.github/plans/diffusion-viscous-shear-flux.md`](../plans/accepted/diffusion-viscous-shear-flux.md)。

### `viscousFlux_wall_d` ([L150](../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L150))

各壁面境界に対し、ゴーストセル `ig` と内部セル `ic` 間の片側差分で粘性応力を構築。
ラッパは `for (auto& bc : msh.bconds)` ループでこのカーネルを発行する。
`wall` / `wall_isothermal` の両 BC が同じこのカーネルを共有する。

ミラーゴースト (`Ux[ig]=-Ux[ic]`、`boundaryCond_d.cu` の `wall_d`) では $\mathbf{d}\parallel\mathbf{S}$
なので over-relaxed の `delta = sss`・`k = 0` となり、内部面と同じ $\tau_{ij}S_j$ を:

```cpp
tau_x  = mu*((Ux[ig]-Ux[ic])/dcc)*sss;              // 法線項 (面積 sss)
tau_x += mu*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz); // 転置項 (セル中心勾配にフル S)
tau_x += -mu*2.0/3.0*divu*sxx;                      // 発散項
```

として構成する。`twall_x_b = tau_x/sss` 等が物理的な壁面せん断 (ストリーム方向) を表す。

> **履歴 (2026-06-06 修正)** — 以前は法線項を成分別に `tau_x ∝ sxx`(`tau_y ∝ syy`,
> `tau_z ∝ szz`)で組み、転置項もコメントアウトされていた。軸平行な水平壁 (法線 = $y$,
> `sxx`≈0) では x 運動量の壁摩擦 `tau_x ∝ sxx = 0` となり、ストリーム方向の no-slip 抗力が
> 一切かからなかった (`twall_x ≡ 0`・`twall_y` のみ非ゼロ)。case 24 で SU2 laminar 参照の
> 放物線に対し ~18 m/s の滑り台座を持つプラグ流・流量約2倍・Mach 過大評価となっていた。
> 修正後は壁隣接セル Ux≈0.24 m/s・中心/平均比 1.53・流量 SU2 比約9%差・`twall_x` 非ゼロ。
> 計画は [`.github/plans/diffusion-viscous-shear-flux.md`](../plans/accepted/diffusion-viscous-shear-flux.md)。

#### node-centered の壁法線項 (`isNode`=`discretization=="node"`, 2026-06-20)

上記の法線項 `mu*((Ux[ig]-Ux[ic])/dcc)*sss` は **cell-centered** を前提とする。cell では
ゴースト中心が `cc_ghost = cc + 2((pc-cc)·n) n` で生成され、cell 中心が壁面から法線距離 $d_n$
だけ内側にあるため `dcc = |cc_ghost-cc| = 2 d_n`。ミラー (`Ux[ig]=-Ux[ic]`) と合わせると
法線項は $(U_{wall}-U_{cell})/d_n$ ($U_{wall}=0$) に一致し、正しい片側壁勾配になる。

**node-centered (median-dual)** では境界 CV 代表点 (=ノード) が**壁面上に乗る**。半割面重心
`pc` も壁面上にあるため $(pc-cc)\cdot n \approx 0$、すなわち `cc_ghost ≈ cc` で **`dcc ≈ 0`** に
退化する。残る微小な `dcc` は壁の曲率に由来し符号も面ごとにばらつくため、法線項
`(Ux[ig]-Ux[ic])/dcc` は $0/0$ 的に爆発し、近壁で偶奇振動・implicit 不安定を引き起こす。
勾配計算 ([`calcGradient_d.cu`](../solver_density_cuda/cuda_forge/calcGradient_d.cu)) は既に
この退化を認識し、壁半割面を cellgather から外して **bvar 弱形式**で閉じている。

そこで `cfg.discretization=="node"` のとき (`viscousFlux_wall_d` に渡す `isNode==1`)、壁法線項を
**node セル勾配ベース** $\mu\,\nabla u_i\cdot\mathbf{S}$ に置換する。$\nabla u_i$ は calcGradient で
bvar 弱形式により壁面で正しく閉じているので、退化距離 `dcc` に依存しない:

```cpp
// node: 法線項 = μ ∇u_i·S (cell 勾配は bvar 弱形式で壁閉包済み)
tnx = mu*(dUxdx[ic]*sxx + dUxdy[ic]*syy + dUxdz[ic]*szz);
// 転置項 μ ∂u_j/∂x_i S_j と体積粘性項 -2/3 μ divu S_i は cell/node 共通で不変。
hflux = (k_lam + cp*mu_turb/Prt)*(dTdx[ic]*sxx + dTdy[ic]*syy + dTdz[ic]*szz);
```

転置項・発散項はもともと node セル勾配にフル `S` を内積する形なので変更不要 (置換は法線項
のみ)。cell モードでは `isNode=0` で従来挙動をビット不変で維持する。
> 注 (2026-06-24): 当初この切替に試作した `nodeWallViscGradFlux` / `nodeWallDistFloorCoef` (dcc floor)
> は不採用となり**撤去済み** (plan §8.2.5)。現在 node の壁法線項置換は `discretization=="node"` のみで
> 決まる (専用フラグは無い)。
診断は env `FORGE_VISC_WALL_DIAG=1` で壁半割面の $d_n$, `dcc`, `dcc/(2d_n)`, 接線オフセットを
集計表示する (`viscousWallDiag_d`)。計画は
[`.github/plans/diffusion-node-wall-viscous-distance.md`](../plans/accepted/diffusion-node-wall-viscous-distance.md)。

#### node-centered の壁摩擦応力 twall = 内部双対面集約 (`wallStressForOutput_node_d`, `nodeWallStressEdgeKernel`, 2026-06-24)

> **改名 (2026-06-27)**: 旧 `viscousFlux_wall_node_d`。残差・`y+` に一切触れず `twall` 出力を後処理として
> 算出するだけのため、その役割が分かる名前 `wallStressForOutput_node_d` に変更。

node の壁ノード $W$ は壁面上に乗るため、`viscousFlux_wall_d` が出力する `twall` は退化ミラーゴースト
`dcc` / 壁ノード勾配 $\nabla U[W]\cdot\mathbf{S}$ ベースで、近壁で特異スパイク・偶奇振動する (上記 §)。

**離散化の非対称性 (重要)**: cell モードでは壁境界面は内部セル $ic$ の CV に属し、`viscousFlux_wall_d`
の残差加算が壁せん断の**唯一の供給源**なので必須。一方 node モードでは、壁せん断は**内部双対面**
$W\leftrightarrow I$ (内部ノード) を介して `viscousFlux_d` (内部カーネル, `ip<nNormalPlanes`) が既に
保存形で内部ノード $I$ の運動量に加算しており ($U[W]\approx0$ Dirichlet で $(U_I-0)/dcc$ が no-slip
せん断)、`viscousFlux_wall_d` の残差加算は**壁ノード $W$ (Dirichlet で破棄) へ捨てるだけ**。よって
node での `viscousFlux_wall_d` の実効的役割は `twall`/`y+` 出力のみで、それが上記退化で不正確。

`cfg.nodeWallStressEdgeKernel==1` (既定) のとき、node の `twall` を `wallStressForOutput_node_d` で
**置換**する (`viscousFlux_d.cu`)。1 スレッド = 1 壁半割面 $ib$ で、$W=$`bplane_cell[ib]` の入射内部
双対面 (CSR `cell_planes`) を走査し、各面の粘性運動量 flux を $W$ の CV に集約して壁半割面積で割る:

$$\boldsymbol{\tau}_w(W)=\frac{1}{A_{wall}(W)}\sum_{I:\,W\text{-}I\,\text{内部面}}\Big[\mu_f\frac{U_I-U_{wb}}{|dcc|}\,\delta\;+\;\mu_f\,\nabla u_i\!\cdot\!\mathbf{k}\;+\;\mu_f\,(\partial_i u_j)S_j\;-\tfrac{2}{3}\mu_f(\nabla\!\cdot u)S_i\Big]$$

- 壁端速度は no-slip 値 $U_{wb}=$`Ux_b/Uy_b/Uz_b`$=0$ (ghost 不使用)。
- 面勾配 $\nabla u_i$, $\mu_f$ は **$\tfrac12(\cdot[W]+\cdot[I])$ の dual 平均**。
- over-relaxed 分解 ($\delta$, $\mathbf{k}=\mathbf{S}-\boldsymbol\delta$) は `viscousFlux_d` と同形だが、距離は退化
  `dcc` でなく **$W$-$I$ 間の物理距離 $|cc_I-cc_W|$** (第一オフ壁セル厚, 退化しない)。
- **運動量残差・`y+` には触れない** (内部ノード運動量は `viscousFlux_d` が担うため二重計上回避;
  `twall` は出力専用で**場は不変**)。cell モードと `nodeWallStressEdgeKernel==0` では従来どおり。

**壁関数時の twall モデル値化 (2026-06-27)**: 上の集約は**解像** $\mu_{tot}\,\partial_n u$ なので、壁関数メッシュ
($y^+\!\sim\!30\text{–}300$, wt=1) では第一セルの $\mu_{turb}$ が大きく真の $\tau_w=\rho u_\tau^2$ を**~十数倍に過大評価**
する (`utau`/`y+` は対数則由来で物理値なのに `twall` だけ解像値で不整合だった)。そこで壁関数 active
($\texttt{Tau\_Wall}[W]=\rho u_\tau^2>0$、`ransWallFunction` が算出、AddTauWall が運動量に課す値と同一) のとき、
`wallStressForOutput_node_d` は**向きを解像 traction・大きさを $\tau_w=\rho u_\tau^2$ に再スケール**して
`twall`・`utau`・`y+` を整合させる。**壁解像 (wt≠1, `Tau_Wall<0`→nullptr) では解像値=真の $\tau_w$ なので無補正**。
→ どちらのモードでも `twall` が物理的な壁せん断になり、`twall` から正しい $C_f$ が出る。
設計判断は [`plans/accepted/output-node-wall-surface-viz.md`](../plans/accepted/output-node-wall-surface-viz.md)。

計画は [`.github/plans/diffusion-node-wall-viscous-distance.md`](../plans/accepted/diffusion-node-wall-viscous-distance.md) §11。

#### node-centered の内部面 dcc に node 座標を使う (2026-06-22)

内部面の over-relaxed 拡散/粘性係数は `delta = |dcc|·|S|²/(dcc·S) = |S|/cosθ` で、`dcc` (CV 間
ベクトル) と面ベクトル `S` の非直交角 `θ` が大きいと `1/cosθ` で発散する。**node-centered では
flux の `dcc` に `ccx`(=`CELLS/centCoords`=双対 CV 面積加重重心, `axisCentroidShift`) を使うと、
近壁の半割双対 CV で重心がノードから最大 ~0.0075 ずれ、見かけの非直交 `1/cosθ` が実測 max
1.46e6 に発散**する。median-dual の双対面は**ノードエッジに直交**するので、`dcc` に**ノード座標**を
使えば `1/cosθ`=1.00 (厳密直交) になる (検証: case/36)。

この見かけ非直交が μ_t × 近壁 omega 勾配 (restart で壁近傍 ω~4.3e8) を増幅し、**SST omega 拡散
flux を爆発**させる (case/36 node SST が step3 で roOmega→1e22)。検証: `dcc` を node 座標で測ると
`1/cosθ`=1.00、内部面拡散/粘性を node 座標 dcc にすると case/36 node SST が step3→1657 へ改善
(ncx 試作で確認、撤去済)。

**根治方針 (採用)**: node 専用配列 (ncx 等) で分岐するのでなく、**`CELLS/centCoords` 自体を「値の位置」=
ノード座標に統一**し、cell/node で同一処理にする。双対 CV 重心は別量に分離し、軸対称の r 重み/source
だけが参照する (`axisCentroidShift` 撤去)。これで内部面は自動で直交化する。ただし `centCoords=node` は
**壁ノードが壁面に乗り ghost mirror の dcc が退化** (検証: NaN 132/132 が壁) するため、**node モードの
境界を完全 ghostless 化**するのが前提。
計画: [`.github/plans/architecture-node-centroid-value-position.md`](../plans/active/architecture-node-centroid-value-position.md)
(旧 [`diffusion-node-scalar-nonortho-limit.md`](../plans/archived/diffusion-node-scalar-nonortho-limit.md) は superseded)。

**スカラ (k/ω)・化学種拡散の node 境界半割面は「加えない (skip)」**: 退化する ghost mirror も、暫定で
使っていた `∇φ·S` 弱形式 (cell 勾配の境界投影に依存・Neumann で厳密 0 にならない) も用いず、`isNode!=0` かつ
ghost を含む半割面 (`ic0>=nCells || ic1>=nCells`) では拡散流束を**スキップ**する (`scalar_diffusion_first_order_d`,
`species_diffusion_d`)。根拠: Dirichlet (固定値入口 / k=0・ω ピン) では境界ノードがピンで上書きされ半割面は無意味で、
内部ノードへの拡散は**内部双対面 W↔I が実距離で運ぶ** (主ループの非境界 plane で計算済)。Neumann (slip / zero-grad)
は物理的に $\partial_n\phi=0$ ＝半割面フラックス 0 そのもの。例外は陽に課す非ゼロ Neumann 熱流束のみ (forge は断熱/
Dirichlet/等温=Dirichlet なので不要、等温壁の熱流束は粘性流束側)。平板 `case/26.flat_plate_sst` で **Cf/u_τ/δ99 が
∇·S 版と完全一致** (差 0.00%) を確認、k は前縁上流スリップ域で skip がより正 (∇·S は Neumann 漏れで過剰生成)。
計画: [`.github/plans/diffusion-node-boundary-real-distance.md`](../plans/accepted/diffusion-node-boundary-real-distance.md)。

### 入出力

入力: セル中心保存量・原始量 (`ro, roUx, …, Ts`)、勾配 (`dUxd*, dUyd*, dUzd*, dTd*`)、
粘性係数 (`vis_lam`, `vis_turb`)、`mu`, `thermCond` (定数)。

出力 (累積): `res_ro, res_roUx, res_roUy, res_roUz, res_roe`。

### 既知の TODO / 注意点

- 法線差分に用いるのはセル中心値のみで、再構成 (MUSCL) は適用しない (拡散項として妥当)。
- 壁面カーネルは現状 `wall` / `wall_isothermal` のみ呼ばれる。`slip` には粘性束を加えない。
- `vis_lam` は配列で渡される構造になっているが、現状はラッパで `cfg.visc` を直接渡している。
  温度依存粘性 (Sutherland 法則など) を導入する場合は配列を実体化する。

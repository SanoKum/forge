# 軸対称定式化 — 実装

[theory.md](theory.md) の B 流儀 (幾何量を $r$ 重み付け) を `solver_density_cuda`
にどう落とし込むかを整理する。Phase 1 のスコープ (非粘性 + 半径運動量の圧力
ソース + 粘性 source + 自動軸 BC + RK 陽解法) のみを対象にする。

## 設定フラグ

| 場所 | 項目 | 既定 | 内容 |
| --- | --- | --- | --- |
| [solverConfig.yaml](../../case/) `physProp:` | `isAxisymmetric` | `0` | `1` で軸対称モード。対称軸は x 軸固定、半径 $r = y$ |

`isCompressible` と同じ層に置く。読込みは
[`input/solverConfig.cpp`](../../solver_density_cuda/input/solverConfig.cpp) の
`physProp` パース節で `getValidatedValue<int>(physProp, "isAxisymmetric", "physProp")`
を追加し、`isCompressible` と同じ要領で本体クラスに格納する。
未指定時は 0 (= 既存 3D ケースと完全互換) を許容する。

## データフロー全体図

```mermaid
flowchart TD
    A[solverConfig.yaml<br/>isAxisymmetric] --> B[variables::setStructuralVariables_d<br/>r 重み付け実施]
    B --> C[c_d/p_d: volume*=ȳ, S*=ȳf]
    B --> D[c_d: A_planar 保存]
    C --> E[convectiveFlux/limiter/setDT/...]
    D --> F[axisymmetricSource_d_wrapper<br/>res_roUy += (p - τ_θθ) * A_planar]
    E --> G[time integration]
    F --> G
```

ポイント:

- **既存対流フラックスカーネル群は `S_f` / `volume` の値を読むだけで、$r$ が
  乗っていることを意識しない**。重み付けの責任は幾何セットアップに閉じる。
- 半径運動量の圧力・粘性ソースは別カーネル `axisymmetricSource_d_wrapper` で加算。
- `isAxisymmetric == 0` のときはすべての分岐が早期 return か無加算で、既存
  ケースの数値挙動・性能に影響しない。

## 検証ノズル形状の位置づけ

`case/23.axi_nozzle/mesh_axisym_m2/` の geometry generator は solver core の一部ではなく、
軸対称検証ケース用のノズル形状生成である。更新後の構成は次の 3 区間で考える。

- 収束部: Witoszynski 曲線
- 喉部直後の遷移部: 5 次多項式ブレンド
- 主発散部: planar または axisymmetric の sharp-throat MOC

比較用の public reference contour として、solver core を変えずに
`case/23.axi_nozzle/mesh_axisym_m4/` 配下へ、収束部も含めて公開式だけで閉じた
bell nozzle の true 2D mesh path を追加して検証することもある。既定の public path は
「直線 chamber + conical convergent + 1.5Rt throat entrant arc + Rao bell divergent +
straight test section」で構成し、既存 in-house 収束部の影響を切り離して比較する。
この path は既存の MOC ノズルを置換するものではなく、段差や局所曲率、上流側 contour の
影響を切り分けるためのケース側比較対象である。

5 次多項式ブレンドの役割は、sharp-throat MOC が持つ喉部直後の壁角不連続を
そのまま使わず、Witoszynski 側の `y, y', y''` から、少し下流の滑らかな MOC 点の
`y, y', y''` へ $C^2$ 連続に受け渡すことである。したがって、MOC 全体を作り直す
のではなく、先頭数点だけを置換する実装になる。

axisymmetric MOC mode では、planar MOC の characteristic invariants を固定せず、
`K_- = theta + nu`, `K_+ = theta - nu` に対する一次 source correction を各 characteristic
長さで積分する。interior / axis / wall 点の幾何交点自体は planar MOC と同じく直線交点で
組み立てるが、状態量 (`theta`, `nu`, `mach`) は fixed-point で数回更新して求める。

## 幾何量の $r$ 重み付け

### 触る箇所

[`variables.cpp::setStructuralVariables_d`](../../solver_density_cuda/variables.cpp)
で host 側にコピーした直後、cudaMemcpy する前に値を上書きする。

```cpp
if (cfg.isAxisymmetric == 1) {
    const geom_float r_floor = (geom_float)1.0e-20;
    for (geom_int ip=0; ip<msh.nPlanes; ip++) {
        const geom_float r_face = (pcy[ip] > r_floor) ? pcy[ip] : r_floor;
        sx[ip] *= r_face;
        sy[ip] *= r_face;
        sz[ip] *= r_face;
        ss[ip] *= r_face;
    }
    for (geom_int ic=0; ic<msh.nCells_all; ic++) {
        const geom_float r_cell = (ccy[ic] > r_floor) ? ccy[ic] : r_floor;
        A_planar_h[ic] = volume[ic];          // 元値が planar 面積 × 単位深さ
        volume[ic]     = volume[ic] * r_cell; // V = ȳ * A_planar
    }
}
```

軸 ($r=0$) 上の face は厳密に `r_face = 0` にすると、下流の SLAU/Roe カーネルが
`n_x = S_x / |S|`, setDT が `dx = vol / |S|` を計算した瞬間に `0/0 = NaN` に
なる。これを避けるため、$r$ に小さなフロア `r_floor = 1.0e-20` を入れる。

- $S_x = S_x^{\text{orig}} \cdot r_{\text{floor}}$, $|S| = |S|^{\text{orig}} \cdot r_{\text{floor}}$
  なので法線比 $n_x = S_x/|S|$ は元の値を保つ (float32 の underflow 域 ~1e-38
  にも入らず正常域)。
- 軸面の flux 寄与 $p \cdot S \approx p \cdot |S|^{\text{orig}} \cdot 10^{-20}$ は
  他の面の寄与に対して 17 桁小さく実質ゼロ。
- ノズル系では $|S|^{\text{orig}} \sim 10^{-3}$、$\bar r \sim 10^{-2}$ が典型なので、
  非軸セル/面では `r_face = pcy[ip]` がフロアを大きく上回り通常パスに乗る。

`r_floor` を $0$ に設定すると軸対称解は走らない (NaN が出る) ので、本値は実装上
の必須パラメータとして扱う。

### `A_planar` の保存

`variables.hpp` に新しいセル変数 `"A_planar"` を追加。`isAxisymmetric == 0` の
ケースでも安全に動くよう、フラグが立っていない時は `A_planar = volume`
(つまり 3D 体積) を入れておく方針も検討したが、混乱を招くので **Phase 1 では
`isAxisymmetric == 1` のときのみ意味を持つ** とし、それ以外では 0 で初期化
する。圧力ソースカーネルもフラグで早期 return するため未参照になる。

### CPU 経路 (`setStructuralVariables`)

GPU を使わないテスト経路 (`cfg.gpu == 0` ブランチ) でも同様の処理を入れる。
ただし forge の運用は GPU 既定なので、CPU 経路の検証は最小限に留める。

## 軸 BC

### 種別追加

[`boundaryCond.hpp`](../../solver_density_cuda/boundaryCond.hpp) の
`valueTypesOfBC` に `"axis"` を **slip と同じ全 `-1` のスタブ**として追加。

```cpp
{"axis", {
    {"ro"  ,-1},  {"roUx",-1}, {"roUy",-1}, {"roUz",-1}, {"roe" ,-1},
    {"Ux"  ,-1},  {"Uy"  ,-1}, {"Uz"  ,-1},
    {"Tt"  ,-1},  {"Pt"  ,-1}, {"Ts"  ,-1}, {"Ps"  ,-1},
}},
```

### dispatch

[`boundaryCond.cpp::applyBconds`](../../solver_density_cuda/boundaryCond.cpp)
で `axis` を **slip と同じカーネル**にディスパッチする。

```cpp
else if (bc.bcondKind == "axis")
    slip_d_wrapper(cfg, cuda_cfg, bc, msh, var, mat_p);
```

理由: B 流儀で軸面 $\mathbf{S}_f = \mathbf{0}$ により対流寄与は元々ゼロだが、
LSQ 勾配計算 (`gradient.cpp`) や limiter は rfweight 前の幾何を見る経路が
一部あり、ゴーストセルの安定値を期待するため。slip ゴーストの構成
(法線速度反転) は軸面に対しても物理的に妥当 (= 鏡映対称)。

### ケース設定

[`case/23.axi_nozzle/run_slau_2d/bcondConfig.yaml`](../../case/23.axi_nozzle/run_slau_2d/bcondConfig.yaml)
の `axis: kind: slip` を `kind: axis` に置換する。検証時は別ディレクトリ
(`run_slau_axisymmetric/`) を作成して既存 run を上書きしない。

## 圧力・粘性ソースカーネル

### 新規ファイル

`solver_density_cuda/cuda_forge/axisymmetricSource_d.cu` /
`axisymmetricSource_d.cuh` を新規作成。

```cpp
__global__ void axisymmetricSource_d(
    geom_int nCells,
    flow_float* Ps,
    flow_float* A_planar,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* axisym_divU,
    flow_float* axisym_uy_over_r,
    flow_float* res_roUy)
{
    geom_int ic = blockIdx.x * blockDim.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float mu_total = vis_lam[ic] + vis_turb[ic];
    const flow_float tau_theta_theta = 2.0 * mu_total * axisym_uy_over_r[ic]
                                     - (2.0 / 3.0) * mu_total * axisym_divU[ic];
    res_roUy[ic] += (Ps[ic] - tau_theta_theta) * A_planar[ic];
}

void axisymmetricSource_d_wrapper(
    solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (cfg.isAxisymmetric != 1) return;   // 早期 return
    axisymmetricSource_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.c_d["Ps"], var.c_d["A_planar"], var.c_d["vis_lam"], var.c_d["vis_turb"],
        var.c_d["axisym_divU"], var.c_d["axisym_uy_over_r"], var.c_d["res_roUy"]);
}
```

### 呼び出し位置

[`main.cpp`](../../solver_density_cuda/main.cpp) の時間進行ループで、
`convectiveFlux_d_wrapper(...)` の **直後**、`viscousFlux_d_wrapper(...)` の
**前**に呼び出す。

```cpp
convectiveFlux_d_wrapper(cfg, cuda_cfg, msh, var, mat_ns);
axisymmetricSource_d_wrapper(cfg, cuda_cfg, msh, var);   // 新規
viscousFlux_d_wrapper(cfg, cuda_cfg, msh, var, mat_ns);
```

理由: 半径運動量 residual `res_roUy` に軸対称ソースを加える操作なので、対流
flux による residual 構築の直後が論理的に自然。face viscous flux 本体とは独立な
cell-centered source なので、`viscousFlux` よりは前に置く。

### 符号規約

forge の対流 flux カーネルは `res_roUy[ic0] -= flux_temp` の形で residual を
"流出 - 流入" として蓄積する。圧力ソースは半径運動量保存式の **右辺** に
$+(p - \tau_{\theta\theta}) \cdot A_{\text{planar}}$ で現れる。時間進行では

$$
V \frac{Q^{n+1} - Q^n}{\Delta t} = -\sum_f \mathbf{F}\cdot\mathbf{S}_f + R
$$

を `Q^{n+1} = Q^n - \Delta t / V \cdot (\sum_f \mathbf{F}\cdot\mathbf{S}_f - R)`
として書くため、`res_roUy` には **$\sum_f F\cdot S_f - R$ を蓄積**する慣習であれば
`res_roUy -= (p - \tau_{\theta\theta}) * A_planar` となる。実装時に既存 flux カーネルの
符号を確認して合わせる (Phase 1B 着手時に
[`convectiveFlux_d.cu`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu)
末尾の `atomicAdd` 符号を確認、コードコメントに符号規約を記す)。

## CFL 時間刻みと軸付近のフロア

[`cuda_forge/setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu) は
$\Delta t = \mathrm{CFL} \cdot V / \sum_f |\mathbf{F}\cdot\mathbf{S}_f|$ を計算する。
$r$ 重み付けで分子・分母とも同じ $r$ が乗るため理論的には不変だが、軸近傍
セルでは離散誤差により $\Delta t$ が病的に小さくなる懸念がある。

Phase 1 の対応: **まずは無対策で走らせ、軸近傍セルの dt が極端に小さく
なって全体収束が阻害される場合のみ**フロア $\bar r_{\text{eff}} = \max(\bar r,
\varepsilon h_{\text{local}})$ を導入する。導入する場合の touch 箇所は

- `variables::setStructuralVariables_d` の `r_cell` 計算で $\max$ 適用
- またはセル `volume` を保存する直前にフロア適用

の二択。$\varepsilon \approx 10^{-3}$ を初期値とする。

### `axisTimestepBeta`: 近軸半径音響 additive 安定化 (2026-06)

理論は [theory.md](theory.md) §"近軸半径音響モードの不足"。config
`time.deltaT.axisTimestepBeta` (既定 `0` = ビット不変・従来挙動)。
[`setDT_d.cu`](../../solver_density_cuda/cuda_forge/setDT_d.cu) の `setCFL_cell_d` で
軸対称セル (`isAxisymmetric==1`) かつ `axisTimestepBeta>0` のとき face スペクトル半径に
半径音響項を加算する:

```cpp
// face 項と同じ無次元化 (dt·λ/V)。V=r·A_planar (per-radian) なので A_planar/V=1/r。
cfl[ic] += dt * axisBeta * (fabsf(Uy[ic]) + sonic[ic]) * A_planar[ic] / vol[ic];
```

設計上の要点:

- **additive (face+axis の和) であり min クランプではない**。min は最内セルしか触れず
  CFL≳4 で不足した (検証で確認)。additive は広い near-axis 域を滑らかに減衰し頑健。
- `dt_local` が短くなると DPLUR 時間項 $V/\Delta\tau\,I$ ([`timeIntegration_d.cu`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)
  `block_dplur`) を通じ**全保存式の LHS 対角が整合的に強まる** (単一式への人工対角でない)。
- `axisTimestepBeta` は**数値スペクトル半径の重みであり物理 CFL 上限ではない**。
  $\Theta_{\text{axis}}=\Delta\tau(|u_r|+c)/\bar r \approx \mathrm{CFL}/\beta$ は face 項にも依存し厳密保証でない。
  `axisAcousticCFLMax` 等の「Θ 直接指定」名称は設定値と実効 Θ が一対一でないため**採らない**。

検証 (`case/28` 多成分 TP, 発達場から継続):

- $\beta_{\text{axis}}=2$ で **CFL 4 安定** (無対策は near-axis で step 1〜50 発散)。$\beta\approx\mathrm{CFL}/2$ は
  `case/28` の経験則 (メッシュ・流れ場依存、一般式でない)。
- **CFL 4→1 down-test で残差が再低下** (rms_ro $4.5\times10^{-6}\to5\times10^{-7}$) し、近軸の
  $u_r$・軸上 $p$ が baseline 低 CFL 解に一致 → 高 CFL の残差床は擬似時間 limit cycle で、
  **定常解は不変**。質量流量・組成場も baseline と一致。
- **高 CFL 残差床 ($\sim4\times10^{-6}$) は軸でなく He/空気 contact 混合層モード**が主因。別 issue。
  なお化学種移流は既に 1 次風上 (`scalarTransport_d.cu` `scalar_advection_first_order_d`) なので
  混合層振動は species 高次再構成ではない (flow 2 次再構成 or 物理せん断層を切り分け中)。

診断用に `FORGE_AXIS_DIAG_ALPHA` (env, 既定 0) があり、roUy 対角へ $\alpha\,A_{\text{planar}}c$ を
直接足す per-equation 版。機構特定の診断専用で **production では使わない** (additive setDT 版を使う)。

## 軸対称 r 床 (`axisRFloor`, 2026-08 — ユーザ提案)

`mesh.axisRFloor` (>0 で有効, 既定 0 = 従来床 1e-20 = ビット不変)。SU2 の y<EPS ガードの
r 重み版: **r < axisRFloor の面・セルは r 重みを床値で打ち切る** (幾何が消えない = 軸半 CV の
真空化が起きない)。床帯のセルは hoop ソース・Jacobian・u_r/r も課さない (planar 扱い)。

**閉性補正 (実装の要)**: 床を面にだけ当てると、床帯と非床帯にまたがる CV で
∮r n dA が閉じず一様圧力でも偽の力が出て即発散する (実測 step 19)。このため hoop ソース/
Jacobian の面積は解析 A_planar でなく**床適用後の離散閉性 Σ_f S_f (outward)** を使う
(`A_closure_x/y`, セットアップ時に全 plane 走査で計算)。非床領域では解析値に一致、全面床の
CV では厳密に 0 (=「ソースも入れない」が自然に導出)、混在 CV では任意の床で一様圧力が
厳密に力ゼロになる。A_planar (勾配計算の分母) は不変のまま分離する。

**検証 (case/40 node, run_0027 + soak 24000 step)**: 軸行健全 (床ピン 0)・η_CF=0.9906・
rms_ro ~1e-7 プラトー。床帯の縁に有界の k 帯 (~230, +12000 step でビット不変の平衡) が残る。
床値はメッシュ依存 (軸行重心 < axisRFloor < 第一内点行重心; case/40 では 3.0e-4)。

## SU2 流定式化 (`axisymMethod: 1`, 2026-08)

`mesh.axisymMethod` で軸対称の定式を選ぶ (既定 0 = 本文書の r 重み方式・ビット不変)。
**`isAxisymmetric`/`axisymMethod` の正本は `mesh` ブロック** (幾何/離散化の設定であり物性では
ないため。2026-08-04 に physProp から移設、physProp 配下は後方互換の deprecated 読み+警告)。
`1` は **SU2 と同じ「planar 幾何 + 1/y ソース項」方式** (plan
[axisymmetric-su2-source-formulation](../../plans/active/axisymmetric-su2-source-formulation.md)):

- 幾何は planar のまま (`variables.cpp` の r 重みをスキップ、`A_planar`=planar 体積)。
- 非粘性ソース $S = -(1/y)[\rho v,\ \rho u v,\ \rho v^2,\ \rho v H]$ と解析 Jacobian
  (block-DPLUR の対角へ)、粘性ソース (SU2 `ResidualDiffusion`: AuxVar $\mu v/y$ 系の GG 勾配
  込み)、SST の 1/y 移流拡散ソース (`rans_sst_source_d` 内、対角に $\max(v/y,0)$)。
  実装は `axisymmetricSourceSU2_d` (`axisymmetricSource_d.cu`)。
- 軸ノード (`axis_flag`) と $y\le\varepsilon$ は SU2 同様ソース 0。軸 BC は従来どおり
  slip 鏡映 + `zeroAxisRadialResidual` (= SU2 `BC_Sym_Plane` 相当)。
- **軸ノードは通常 DOF として解ける** (軸 CV は有限 planar 体積)。`nodeAxisDirichlet` は不要。
- SU2 コードの 2 点は意図的に移植しない: エネルギー粘性ソースの `+ρk` 項 (次元不整合の疑い)、
  エネルギー Jacobian の `1/2` 整数除算 (数学的に正しい形で実装)。

**検証状況 (case/40 node, 2026-08-04)**: 軸行は真に健全 (床ピン 0・k シート消滅・軸線 T 滑らか)、
12000 step 完走で η_CF=0.9904 (method 0 の 0.9896 と 0.08% 差)・quasisteady ALL STEADY。
ただし**喉部近軸の小さなリミットサイクルで rms_ro ~6e-5 プラトー** (method 0+`nodeAxisDirichlet`
は 1e-7〜1e-9 に達する)。疑いは既知の node slip 市松バグ ([[node-slip-spurious-flow]]) が
planar 面積の軸 slip で顕在化したもの。**このため生産既定は当面 method 0 + `nodeAxisDirichlet`**
とし、slip 市松の修正後に method 1 への既定切替を再評価する。

## node-centered 軸ノードの対称 Dirichlet (`nodeAxisDirichlet`, 2026-08)

node (median-dual) 軸対称では軸上ノードが半 CV の中心になる。この軸半 CV を通常どおり
解くと、強い膨張を持つノズル (case/40: Pt 4 MPa, ε9) のベル部で**軸行だけが radial 圧力
平衡からデカップルして真空まで過膨張する** (ro が隣接行の 1/10〜1/50、有効圧力負 →
EOS 床 P=1 Pa ピン)。laminar でも治癒せず (基底スキームの欠陥)、SST では崩壊軸行と第一
内点行の間の偽せん断が k 生産シート (k ~ 周囲の 10³ 倍) を作り、陰解法の大 pseudo-CFL で
この平衡が弱不安定化して遅発性発散する (詳細:
[plans/active/boundary-node-nozzle-wall-outlet-stability.md](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) §2.6)。

対策は**軸ノードを解かず、対称条件 ∂q/∂r=0・u_r=0 の 1 次離散化として radial 隣接
ノードからの Dirichlet に置換する** (`mesh.nodeAxisDirichlet: 1`、既定 0 = 従来)。
実装様式は forge の `nodeWallDirichlet`/入口スカラーピンと同じパターン。

**SU2 との関係 (正確に)**: SU2 も同じ頂点中心 median-dual だが、**軸ノードは通常 DOF として
解く** (`BC_Sym_Plane` = 鏡映 flux + 法線運動量残差射影のみで、Dirichlet 置換はしない)。
SU2 の軸が壊れないのは軸対称の定式が **planar 幾何 + 1/y ソース項方式**
(`CSourceAxisymmetric_Flow`: `Coord_i[1]>EPS` で有効・軸上はソース 0) であり、軸 CV の
体積・面積が普通の planar 値を保つため。forge は r 重み幾何 (本文書「幾何量の r 重み付け」)
のため軸半 CV の体積・面積が r̄ ∝ Δr に比例して消え、軸で離散平衡が悪条件化する。
`nodeAxisDirichlet` はこの **forge の幾何方式に固有の脆弱性への対症であり SU2 より強い措置**。
軸 CV を解けるようにする根治 (r 重み幾何の軸極限の精査) は将来課題。

- **代表点**: 軸ノード A ごとに内部双対面の相手 I のうち radial cos 最大の非軸ノードを
  変換時でなく solver 起動時に選ぶ (`mesh.cpp` → `axis_rep_d`)。
- **状態ピン** (`enforceAxisMirror_d`): 毎ステージ `assembleResidual` 冒頭で
  ro/roUx/roUz/roK/roOmega[A]←[I]、roUy[A]=0、roe[A]←roe[I]−½roUy[I]²/ro[I]。
- **残差除外** (`zeroAxisAllResiduals_d`): 軸ノードの res_ro〜res_roe (+RANS スカラー) を 0 化。
- **block-DPLUR**: axis_flag 渡し時の decouple を全 5 行単位行化に拡張し
  `nodeAxisDirichlet==1` のとき `axis_flag_d` を渡す (従来は nullptr でビット不変)。
- 軸 CV 体積は O(r̄Δr) と極小のため保存への影響は無視できる (床ピンで既に非保存だった)。

cell モード・平面 2D・`nodeAxisDirichlet=0` はビット不変。species/凝縮スカラーの軸ピンは
未対応 (必要時に拡張)。

## node 軸ノードの値の位置と軸半 CV の整合 (`nodeValueAtNode`, 2026-08-16 試作)

`nodeAxisDirichlet` が対症だった軸行欠陥の**離散的な正体**を case/42 の演算子テストで確定した
([case/42.node_axis_dof](../../case/42.node_axis_dof/README.md), plan
[architecture-node-centroid-value-position](../../plans/active/architecture-node-centroid-value-position.md))。

**演算子テスト**: $\rho=$const, $u=U_0-2ax$, $v=ar$, $p=$const は $\partial_x u+\tfrac1r\partial_r(rv)=0$ を満たす
(連続式残差 ≡ 0)。直管の一様 quad (h=0.02/0.01/0.005) で陽解法 1 step の $(\rho^1-\rho^0)/\Delta t$ から
無次元残差 $e=\mathrm{res}_\rho/(V\rho a)$ を測る (`limiter: 0` は線形場を厳密再構成するため意図的)。

| 変種 | 軸 CV | 第一内点行 | 第二内点行 | 内部 |
|---|---|---|---|---|
| 現行 (`axisCentroidShift: 1`, GG, 幾何 fx, MUSCL) | **+0.833** | −0.167 | +0.031 | 0 |
| 〃 + `nodeMidpointFx: 1` | +0.500 | −0.063 | 0 | 0 |
| 〃 + `gradLSQ: 2` | +0.500 | −0.125 | +0.031 | 0 |
| 1 次風上 (`convMethod: 0`) | 0 | 0 | 0 | 0 |
| **`nodeValueAtNode: 1`** (GG or LSQ, MUSCL) | 0 (≤2e-3) | 0 | 0 | 0 |

いずれも **h に依らない定数** = 軸半 CV の O(1) 不整合。原因は「状態はノード値 (u_r=0 のピン・出力・BC) だが
再構成基点 `cpdx`・fx・勾配は双対重心 $\bar r=\Delta r/4$」という**値の位置の混在**。軸半 CV の質量式は

$$\dot\rho_A=-\Big[\tfrac{4}{\Delta r}(\rho v)_{\rm top}+\partial_x(\rho u)\Big],\qquad (\rho v)_{\rm top}=\tfrac12(\rho v)_I\ \Rightarrow\ 2\partial_r(\rho v)$$

でノード基準なら L'Hôpital 極限に厳密一致する (1 次風上が exact な理由) が、重心基準の再構成は
$(\rho v)_{\rm top}$ を過小評価する。fx=0.5 単独や LSQ 単独では +0.5 が残る。

**`mesh.nodeValueAtNode: 1`** (node 限定、既定 0 = ビット不変): solver 読込時 (`mesh.cpp readMesh`, ゴースト生成前)
に `cells[ic].centCoords ← nodes[ic].coords`、置換前の双対重心 y を `mesh::rEff` に退避し、r 重み
(`variables.cpp` の `r_cell`) だけが `rEff` を使う (軸ノードで回転体積が消えない)。ゴースト鏡映・fx・LSQ・
`cpdx` は全てノード基準になる。境界半割面ではノードと鏡映ゴーストが面上で退化するため
`calcStructualVariables_d` の fx を `0/0 → 0.5` にガード (非退化面はビット不変)。

**ノズル検証 (case/41 生産モデル Euler, case/42 run_0002〜0010)**: 軸 DOF を解く (`nodeAxisDirichlet: 0`) と
真空崩壊は起きないが軸 M が −0.5〜−0.8% 欠損 (偶関数外挿との差 max 0.033)。`nodeValueAtNode: 1` で ≤0.007、
近軸場の cell 参照との差 0.0067 (Dirichlet) → 0.0031、残差は 1.9 → **2.7 桁**。

**境界ゴーストの回転体積潰れ (要修正点・対応済)**: 値位置=ノードでは境界ノードが境界面上に乗るため鏡映ゴースト
(`cc_g = cc + 2((pc−cc)·n)n`) が**同位置に退化**する。軸∩境界 (r=0) ではゴーストの回転体積が r 床 (1e-20) まで
潰れ、`setDT` の `dx = vol/|S|` が ~1e-20 → 局所 CFL 1e13〜1e14 → `dt_local` ~1e-22 となって
**入口軸・出口軸の 2 ノードが初期値のまま完全凍結**する (case/42 run_0003 で ro がビット一致)。
対応: `variables.cpp` の r 重みで**ゴースト CV は所有 CV の r̄ を使う** (bcond の `iCells`/`iCells_ghst` から写像)。
根治は ghost 撤廃 (plan §5 Step 1)。

**粘性 (NS) は未対応 — 2 次再構成が壁スリバーで破綻**: case/41 NS 生産モデル (y+1.5, `wall_first_frac` 4.5e-5,
SST) で `nodeValueAtNode: 1` は step 62–64 で発散する。切り分け (case/42 run_0009〜0017):
`nodeAxisDirichlet: 1` 併用でも同 step 64 = **軸 DOF は無関係**、1 次 (`convMethod: 0`) × cfl 1.0 は完走、
2 次 × cfl 0.2 は step 154 (**CFL は遅らせるだけ**)、**粗い壁格子 (5e-3) では 2 次でも完走**、
発散起点は第一内点の P=1.11 MPa > 室圧 1.0 MPa のオーバーシュート。機構は、値位置が壁面上に来ることで
W→内部面の外挿距離が第一セル厚の半分→全長に倍増し、CV 平均勾配+リミッタの外挿が高 AR スリバーで過大に
なること (CV 平均量を点値として扱う不整合。SU2 は点値ベースで LSQ・リミッタも点値)。NS 対応には
(a) 再構成の点値整合 (勾配・リミッタを含む)、(b) `wall_dist` の node 基準再計算 (変換器は双対重心基準)、
(c) ghost 撤廃、が要る。**段階起動では回避できない**ため config に警告を出す。

**軸量の抽出**: 生産 (`nodeAxisDirichlet: 1`) の軸ノードは第一内点コピーで $\tfrac12 q_{rr}r_1^2$ のバイアスを持つ
(case/41 で ~0.07% $M_d$)。設計チェーンは軸直読でなく `forge_design/evaluate/axis_extract.py` の偶関数外挿
($q=a_0+a_2r^2$, r>0 の 4 点, 軸ノード除外) を既定にした (`axis_curve_node(mode="evenfit")`)。

## Roe スキームの取扱い

軸対称固有の固有分解 ([theory.md](theory.md) §"Roe フラックスヤコビアン") は
**実装としては既存 3D Roe をそのまま流用する**。理由:

- B 流儀では $\mathbf{S}_f$ が $r$ 重みされている以外、面ごとの flux 計算
  ロジックは Cartesian と同一 (面法線方向の 1D Riemann 問題に帰着)。
- 軸対称で落ちる「$u_\theta$ 方向のせん断波」は、面ベクトルが $z$ 成分を
  持たない 2D メッシュでは波速 0 として自動的に縮退する。実装上は 5 波の
  ままで残しても数値解に影響しない。
- 軸対称専用の 4 波ヤコビアンを別途実装すると保守コストが高く、Phase 1 の
  目的 (B 流儀の動作確認) に対して過剰。

将来の Phase 3 (軸対称陰解法 Jacobian) で 4 波の固有分解を直接使う必要が
出れば、その時点で別実装を立てる。

## 陰解法 (block DPLUR) での軸対称ソースヤコビアン (2026-06)

block DPLUR 陰解法 (`timeIntegration==11`) は、平均流の `A^{+}/A^{-}` 分割修正と
$r$ 重み付け幾何の整合だけで軸対称ケースが収束する（ソースを lagged 扱いでも可）。
さらに半径運動量ソース $S_{\rho u_r}=(p-\tau_{\theta\theta})A_{\text{planar}}$ の局所ヤコビアンを
対角ブロックに加えて収束を改善する。

- 実装: [`cuda_forge/timeIntegration_d.cu`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)
  の `implicit_defect_correction_block_d` に引数 `int isAxisymmetric`・`flow_float* A_planar` を追加。
  面ループ後・`solve_5x5` 前に `isAxisymmetric==1` のとき roUy 行 ($i=2$) へ
  $D_{2,\cdot}\mathrel{+}= -A_{\text{planar}}\,\partial(p-\tau_{\theta\theta})/\partial\mathbf{Q}$ を加算
  （$r_{\text{eff}}=V/A_{\text{planar}}$、$\mu_{\text{total}}=\mu_{\text{lam}}+\mu_t$）。理論式は [theory.md](theory.md) §"陰解法との連成"。
- 粘性フープ項 $2\mu/(\rho r_{\text{eff}})$ が対角 $D_{2,2}$ を正に強化し軸近傍の stiff 性を安定化。
  非粘性 ($\mu=0$) では圧力ソースヤコビアンのみが残る。
- `timeIntegration_d_wrapper` の起動で `cfg.isAxisymmetric` と `A_planar`（非軸対称時は `volume` ダミー、
  カーネル側で未使用）を渡す。
- 検証: `case/23.axi_nozzle` M4 ノズルで陽解法収束解と壁面静圧が一致（平均 0.02%）、ソースヤコビアンで
  過渡収束 ~2 倍速・回帰なし。擬似 CFL 上限は超音速始動の case 律速（planar でも同様に発散）。

## 出力時の積分量

mdot や推力など revolved 量を物理的に正しく出すには $2\pi$ を乗じる必要が
ある。本ソルバ内部 ($V$, $\mathbf{S}_f$) は $2\pi$ を付けない規約で統一して
いるため、出力後処理で必要に応じて掛ける方針とする。

該当箇所の改修は Phase 1 の必須項目ではないが、`output/` および
`post_tool/` の関連スクリプトに「軸対称ケースでは結果に $2\pi$ を乗じる」
旨をコメントとして残す。

## ソース対応表

| 機能 | ファイル | 変更内容 |
| --- | --- | --- |
| フラグ追加 | [`input/solverConfig.hpp`](../../solver_density_cuda/input/solverConfig.hpp) | `int isAxisymmetric;` メンバ追加 |
| フラグ読込 | [`input/solverConfig.cpp`](../../solver_density_cuda/input/solverConfig.cpp) | `physProp` パースに 1 行追加 |
| $r$ 重み付け | [`variables.cpp`](../../solver_density_cuda/variables.cpp) | `setStructuralVariables_d` 末尾で sx,sy,sz,ss,volume を $r$ 重み、`A_planar` 充填 |
| `A_planar` 変数 | [`variables.hpp`](../../solver_density_cuda/variables.hpp) | セル変数リストに `"A_planar"` 追加 |
| 圧力ソース | `cuda_forge/axisymmetricSource_d.{cu,cuh}` | **新規** |
| 圧力ソース呼出 | [`main.cpp`](../../solver_density_cuda/main.cpp) | `convectiveFlux_d_wrapper` 直後に挿入 |
| 陰解法ソースヤコビアン | [`cuda_forge/timeIntegration_d.cu`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) | **2026-06** `implicit_defect_correction_block_d` に `isAxisymmetric`/`A_planar` 引数と roUy 行のソース対角を追加 |
| 軸 BC 種別 | [`boundaryCond.hpp`](../../solver_density_cuda/boundaryCond.hpp) | `"axis"` を `valueTypesOfBC` に追加 (slip 互換スタブ) |
| 軸 BC dispatch | [`boundaryCond.cpp`](../../solver_density_cuda/boundaryCond.cpp) | `applyBconds` の if-else に 1 行追加 |
| ケース設定 | [`case/23.axi_nozzle/`](../../case/23.axi_nozzle/) | `run_slau_axisymmetric/` を新規作成、`solverConfig.yaml` に `isAxisymmetric: 1`、`bcondConfig.yaml` を `kind: axis` に |

## 既存ケースへの影響

`isAxisymmetric` が未指定 / `0` の場合、上記すべての分岐が無効になり、既存
3D / 平面 2D ケース ([`case/02.airfoil`](../../case/02.airfoil/) ,
[`case/21.naca_2d`](../../case/21.naca_2d/) ほか) の数値挙動・性能は変化しない
ことを Phase 1 の検証 §"リグレッション" で確認する。

## Phase 1 検証結果と既知の課題

Phase 1 の実装について、`case/23.axi_nozzle` 周りでは次を確認した。

- `run_slau_axisymmetric_off_check/` で `isAxisymmetric: 0` に戻したケースは、
  元の `run_slau_2d/` と outer_begin residual が float32 の LSB 差を除いて一致した。
  少なくとも今回の分岐追加は既存 2D/3D 経路を壊していない。
- `run_0022_slau_axisym_m4_aximoc_2d_cfl0p8_10k_blend2/` は良好に収束しており、
  軸対称モード自体は安定に動作する基準例として扱える。
- `run_slau_axisymmetric/` は 200 step までは安定に走り、`residual_history.png`
  も生成済み。残差の絶対値が元ケースより約 100 倍小さくなるのは、
  ノズル半径スケール $r \sim 10^{-2}$ を体積・面積に掛けた B 流儀の期待通り。
- 軸近傍での大きな見かけ上の `dPdy` は、weighted geometry をそのまま
  Green-Gauss 勾配に入れたことが原因だった。現在は
  [`calcGradient_d.cu`](../../solver_density_cuda/cuda_forge/calcGradient_d.cu) で
  勾配再構成には planar geometry (`A_planar`, `sx_planar`, `sy_planar`,
  `sz_planar`, `ss_planar`) を使い、軸対称補正項は別扱いにしている。
- 診断用に `axisym_r`, `axisym_p_over_r`, `axisym_uy_over_r`,
  `axisym_divU`, `axisym_roUy_source` をセル出力へ追加し、`main.cpp` から
  CSV dump でも確認できるようにした。これにより raw gradient と
  axisymmetric correction を切り分けて点検できる。

以前の長時間ランで落ちた現象は、後から面ベクトルの向きが逆だったことに起因する
ケース側の不整合として整理している。したがって、この落ちは軸対称モデルの
剛性や $p/\bar r$ の本質的な問題としては扱わない。

現在の整理は、`run_0022_slau_axisym_m4_aximoc_2d_cfl0p8_10k_blend2/` のような
良好収束ケースを基準に、軸対称モード自体は安定に動作するという前提で進める、
というものになる。具体的には次を確認済みである。

- 低背圧 + `1stUpwind` の `run_0009_slau_axislowbp_1st/` は 6000 step 安定に収束。
- その収束解を初期値に使った `run_0010_slau_axislowbp_muscl_restart/` は MUSCL に
  戻しても安定に継続可能。
- `mesh_axisym_m2/` で作り直した軸対称ターゲット Mach 2 ノズルを用いた
  `run_0011_slau_axisym_m2_1st/` も 6000 step 安定に収束。

したがって現時点の整理は、「軸対称の幾何重み付け・圧力ソース・axis BC・raw gradient
の分離」は完了しており、`run_slau_axisymmetric/compare_isentropic.py` で
1D 等エントロピー比較も再現可能になった。残課題は、比較の自動化拡充と、
ケース側の面向き・境界設定のチェック手順を明文化することになる。Phase 2 の候補は以下。

- 1D 等エントロピー解および既存 3D 軸対称ケースとの比較を拡充する
- ケース設定の向き確認 (面ベクトル、境界法線、軸面向き) を点検項目としてまとめる
- 粘性の軸対称項 (`\tau_{\theta\theta}` を含む) を追加する

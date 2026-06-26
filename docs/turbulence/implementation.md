# SST 乱流モデルの実装方針

## 1. 目的

本章では、forge の `solver_density_cuda` に低 Re SST を組み込む際の
実装責務と変更箇所を整理する。理論は
[theory.md](theory.md) を参照し、本章ではコード上の配置と責務分担に絞る。

## 2. 実装の基本方針

- 初期導入は **explicit advection-only** に限定する
- 軸対称は SST 固有の diffusion / source を子 plan に分離する
- 既存 LES (WALE) と laminar 経路を壊さない形で dispatcher を拡張する
- 既存の `wall_dist`, `vis_turb`, 勾配再構成、出力系を最大限流用する
- 5 方程式 NS 本体と scalar transport を分離し、`k`, `omega` は後者の最初の 2 本として扱う
- scalar diffusion は face ベースで実装し、`vis_lam + sigma * vis_turb` を使う
- `scalarTransport_d.*` は共通輸送コアに限定し、SST 固有の boundary / source / closure は model-specific layer に分ける

## 3. 主要変更点

### 3.1 状態変数

`variables.hpp` のセル変数に次を追加する。

- 保存変数: `roK`, `roOmega`
- stage / backup 変数: `roKN`, `roOmegaN`, `roKM`, `roOmegaM` など
- 勾配: `dKdx`, `dKdy`, `dKdz`, `dOmegadx`, `dOmegady`, `dOmegadz`
- 残差: `res_roK`, `res_roOmega`

初期フェーズでは implicit 用の `dq_*`, `rhs_block_*`, `diag_block_*` は
拡張しない。

### 3.1.1 初期乱流場 (`kInit` / `omegaInit`)

IC (`valueFileName`) に `roK`/`roOmega` が無い場合 (非 SST 場・メッシュ変換直後など)、
`variables::readValueHDF5` は `roK = ρ·kInit`, `roOmega = ρ·omegaInit` で初期化する
(`solverConfig` の `turbulence.kInit` / `turbulence.omegaInit`、既定 0)。

**`ω=0` 初期化は `μ_t = ρ a₁ k / max(a₁ ω, S F₂)` が ill-posed で SST が cold start から
発散しやすい** (特に node `wallTreatmentSST=1`: ω₀=0 だと inlet コーナーで ω が暴走)。
非 SST IC から SST を始めるときは inlet と同じ freestream 値 (例 `kInit=1.0`, `omegaInit=1000.0`)
を設定して 0 初期化を避ける。既定 0 では `ρ·0 = 0` で従来動作 (ビット不変)。
k/ω を持つ収束場から restart する場合は不要 (IC の roK/roOmega をそのまま使う)。

### 3.2 設定

`solverConfig` の turbulence 設定は、現在の `LESorRANS` と `LESmodel` に加えて、
RANS 用モデル選択と自由流乱流量入力を持てるよう拡張する。

初期実装で想定する入力項目:

- `LESorRANS`: `0=no model`, `1=LES`, `2=RANS`
- `LESmodel`: 既存 WALE 用
- `RANSmodel`: `1=SST`
- `turbIntensity` または直接 `kInf`, `omegaInf`
- `turbLengthScale` または等価な自由流尺度

### 3.3 メインループ

メインループの基本順序は
[docs/architecture/overview.md](../architecture/overview.md) を維持する。

初回マイルストーンの概念的な順序は次である。

1. 保存変数から primitive 量を更新
2. 境界条件適用
3. NS 用の勾配計算
4. SST から `vis_turb` を更新
5. NS の対流流束
6. generic scalar transport による $k$, $\omega$ の対流流束
7. model-specific scalar boundary / source layer
8. NS の粘性流束
9. explicit time integration

既存コードでは `convectiveFlux_d_wrapper` と `viscousFlux_d_wrapper` が
NS のみを想定しているため、初期実装では `scalarTransport_d.*` を追加し,
まず 1 次 upwind の advection-only を別 wrapper として差し込んだ。
その後、同じ wrapper に face ベースの diffusion を追加した。
ただし source と壁面生成ロジックは `scalarTransport_d.*` に埋め込まず、
`ransBoundary_d.*` と `ransSource_d.*` のような model-specific layer に分離する。
さらに k/ω への適用 (descriptor 構築・有効化判定・勾配) も共通コアから切り出し、
`ransTransport_d.*` にまとめた。

### 3.4 渦粘性

`turbulent_viscosity_d_wrapper` は現在 WALE 専用である。
ここを次の dispatcher に再編する。

- `LESorRANS == 0`: `vis_turb = 0`
- `LESorRANS == 1`: 既存 WALE
- `LESorRANS == 2 && RANSmodel == 1`: SST

### 3.5 境界条件

境界条件は既存の `wall`, `wall_isothermal`, `inlet_*`, `outlet_*`, `slip` などに
対して turbulence scalar を追加定義する。

ただし責務は 2 層に分ける。

- 既存 NS BC 層: 速度・圧力・温度など主流れ変数を処理する
- turbulence scalar BC 層: `kb`, `omegab` の生成と `kg`, `omegag` の ghost 反射を処理する

初期方針:

- 壁: `wall` / `wall_isothermal` では turbulence scalar BC 層がまず `kb`, `omegab` を生成し,
	ghost セルへ反映する。等距離 ghost を前提にすると、Dirichlet wall 値 $\phi_w$ に対して
	$\phi_g = 2\phi_w - \phi_c$ を使うため、$k_w = 0$ なら $k_g = -k_c$ になる。
	`\omega_w` は壁距離と粘性から与え、ghost 側は $\omega_g = 2\omega_w - \omega_c$ で閉じる。
- 流入 / farfield: 既知の `kb`, `omegab` を与え、ghost 側は同じ反射式で閉じる
- 流出: 1 次外挿ベース
- `slip` / `axis` / `periodic`: no-gradient / copy を使う

現行コードの対応は次のとおり。

- NS BC 層が既存の `wall_d` / `wall_isothermal_d` / `inlet_*` / `outlet_*` / `slip_d` / `periodic_d` を持つ
- turbulence scalar BC はこの後段に差し込む `ransBoundary_d.*` へ移していく
- `bcondConfFormat` の `wall` / `wall_isothermal` に `kb` と `omegab`、`inlet_*` に `k` と `omega` を持てる

この説明は、壁面がセル中心と ghost セルの中点にある等距離配置を前提にしたもの。
現在の実装では、`wall_d` / `wall_isothermal_d` が `kb` / `omegab` を壁値として受け取り、
ghost セル側へ反射して入れる。
将来、非対称な距離配置にするなら、理論側の距離比の式へ切り替える。

### 3.6 出力と診断

`output_cellValNames` に次を追加する。

- `roK`, `roOmega`
- 必要なら primitive として `k`, `omega`
- 既存 `vis_turb`

`gatherResidualSnapshot` は現状 5 方程式固定なので、初回段階では既存の
主流れ residual を壊さないことを優先し、必要なら scalar residual は別ログで扱う。

### 3.7 automatic (enhanced / y⁺ 非依存) wall treatment

理論は [`theory.md`](theory.md) §6.5。`wallTreatmentSST`
(`solverConfig` の `turbulence.wallTreatmentSST`, 0:wall-resolved[既定] / 1:automatic)で
opt-in する。既定 0 では §6.1 の wall-resolved 型 (`60ν/β₁y²`) を保ち、既存結果は不変。

責務分担とソース対応:

- **`ransWallFunction_d.*` (新設)**: 壁面ごとに Reichardt 則の逆解き (Newton 数回) で
	摩擦速度 $u_\tau$ を解き、`bc.bvar_d` の `utau` / `ypls` / `twall_*` を埋める。
	$u_\tau$ は速度・`wall_dist`・$\nu$ のみで決まるため勾配計算前に評価でき、後段の
	$\omega$ 壁 BC と運動量壁せん断の両方で共有する。`applyRansScalarBoundaries`
	(=`ransBoundary_d_wrapper`) の壁ループ直前で `wallTreatmentSST==1` 時のみ呼ぶ。
	`ransWallFunction_d.*` はあわせて新セル変数 `wf_pk`
	(wall-function 生産 $P_k = \rho u_\tau^4/\nu\cdot g(1-g)$, $g=du^+/dy^+(y^+_1)$) を
	wall-adjacent セルに書く (非 wall セルは `-1`)。
- **`ransBoundary_d.*` の wall scalar カーネル**: `wallTreatmentSST`, `utau`, `roOmega` を受け、
	mode 1 で $\omega_w = \sqrt{\omega_{vis}^2 + \omega_{log}^2}$
	($\omega_{vis} = 6\nu/\beta_1 y^2$, $\omega_{log} = u_\tau/(\sqrt{\beta^*}\kappa y)$) を
	**wall-adjacent セル値 `omega[ic]`/`roOmega[ic]` にピン留め** (§6.5(b)) し、$k$ は
	zero-gradient ($k_g=k_c$, §6.5(b'))。mode 0 は現行 $60\nu/\beta_1 y^2$ と $k_w=0$ Dirichlet。
- **`ransSource_d.*`**: `wallTreatmentSST` と `wf_pk` を受け、mode 1 かつ `wf_pk>=0` のセルで
	標準 P_k を `wf_pk` に置換 (§6.5(d))。さらに ω ピン留めセルの `res_roOmega`/`src_jac_omega` を
	0 化 (Dirichlet セルの残差は 0; `rms_roOmega` 汚染回避)。
- **`viscousFlux_d.*` の `viscousFlux_wall_d`**: `wallTreatmentSST` と `utau` を受け、
	mode 1 で接線せん断を modeled $\boldsymbol{\tau}_w = \rho u_\tau^2 \hat{\mathbf e}_t$ に
	置換 (法線粘性項・熱流束は不変、no-slip なので壁せん断仕事 0)。`twall_*_b` / `ypls_b` は
	この modeled 値で上書き出力。mode 0 は現行の分子勾配式。

`utau` は `wall` / `wall_isothermal` の `bvar` 初期化リスト (`boundaryCond.hpp`) と
`mesh.hpp` の `bplaneValNames` マスターリストの**両方**に追加する (片方だと device 未確保で
illegal memory access)。`wf_pk` は `variables.hpp` の `cellValNames` に登録し、
`boundaryCond.cpp` の `applyRansScalarBoundaries` が bc ループ前に
`initWallFunctionPk_d_wrapper` で全セル `-1` に初期化する。`ypls`, `twall_x/y/z` は既存。

`wallTreatmentSST==1` の自律収束: ω ピン留め (μ_t 上限) + wall-function P_k + k zero-gradient の
3 点が揃うと、$k$ は生産・消滅平衡 $P_k=\beta^*\rho k\omega_w$ から平衡値 $u_\tau^2/\sqrt{\beta^*}$
に収束し暴走しない。いずれか欠けると k 暴走 (zero-grad 単独) または過小 (P_k 解像勾配のまま) になる。

### 3.8 SST-DES (DDES) length scale 修正

理論は [`theory.md`](theory.md) §8、計画は
[`turbulence-iddes-sst.md`](../../.github/plans/turbulence-iddes-sst.md)。`DESmode`
(`solverConfig` の `turbulence.DESmode`, 0:RANS[既定] / 1:DDES / 2:IDDES[Phase 2 未実装]) で
opt-in する。**`DESmode==0` では DES カーネルを一切呼ばず、`ransSource` も従来 $D_k$ 分岐を
通るため既存 SST とビット不変** (検証で確認済み)。

設定・変数・ソース対応:

- **`solverConfig.hpp/.cpp`**: `DESmode` (0..2 検証)、`C_DES_kw=0.78`, `C_DES_ke=0.61` を追加。
	`wallTreatmentSST` のパターンに倣い `getOptionalValidatedValue` で読む。
- **`variables.hpp`**: `cellValNames` に静的量 `delta_les` (Δmax) と毎 step 診断 `l_des` /
	`fd_shield` / `rd_des` を追加 (自動確保)。`output_cellValNames` にも 4 つ追加し HDF5 出力する
	(診断必須・float32 で $f_d$ の 0/1 張り付きは NaN より危険)。
- **`variables.cpp` (`setStructuralVariables_d`)**: 幾何セットアップ時に host で
	$\Delta_\mathrm{max}=\max_{隣接面}|\mathbf{cc}_{ic}-\mathbf{cc}_{jc}|$ を面ループで 1 回計算し
	`delta_les` へ H2D 転送。実行時 `ccx`/`ccy`/`ccz` (= CV 重心) のみ使い primal mesh を直接
	参照しないので node-centered でも双対 CV の Δ が得られる。
- **`turbulent_viscosity_d.cu` (新カーネル `compute_des_length_d` + device 関数
	`compute_rd_ddes`/`compute_fd_ddes`)**: $r_d$, $f_d$, $l_\mathrm{DDES}$ を計算し
	`l_des`/`fd_shield`/`rd_des` に書く。**依存順が最重要**: $r_d$ は $\nu_t=$`vis_turb`/ρ を
	読むため、`turbulent_viscosity_d_wrapper` 内で `sst_eddy_viscosity_d` (vis_turb 計算) の
	**直後**に `DESmode>0` でのみ呼ぶ。$l_\mathrm{RANS}=\sqrt{k}/(\beta^{*}\omega)$ ($\beta^{*}=0.09$,
	§8.1 の整合)。float32 ガード (grad/分母 floor + $r_d$ clamp[0,10]) を内蔵。
	カーネルは `cfg.discretization` を参照せず CV 抽象 (volume/wall_dist/勾配/vis_turb) のみに
	依存させ、node 化を検証だけで済むようにする。
- **`ransSource_d.cu` (`rans_sst_source_d`)**: 引数に `int DESmode`, `flow_float* l_des` を追加。
	`DESmode>0` で $D_k=\rho k^{3/2}/l_\mathrm{des}$, 陰解法対角 $\partial D_k/\partial(\rho k)=1.5\sqrt{k}/l_\mathrm{des}$
	に分岐。`DESmode==0` は従来 $\beta^{*}\rho k\omega$ / $\beta^{*}\omega$ をそのまま (ビット不変)。
	$\omega$ 方程式・wall-function P_k 置換・ω ピン留めは DES と独立に不変。

呼び出しフロー (main.cpp): `calcGradient` → `turbulent_viscosity` (vis_turb→l_des) →
`ransTransport` → `ransSource` (l_des を読む)。

> **build 注意**: `solverConfig.hpp` に `DESmode` を追加したため、差分ビルドだと CUDA obj を
> 取りこぼし step0 で dt=0/NaN 凍結することがある (既知の罠)。config 変更後は full rebuild する。

## 4. Generic Scalar Transport 基盤

`k`, `omega` の輸送は、将来の passive scalar や species でも再利用できる
generic scalar transport 層で扱う。

実装上の責務分担は次とする。

- NS 本体: `convectiveFlux_d.*`, `viscousFlux_d.*`
- scalar 共通輸送コア (物理非依存の移流・拡散・時間積分): `scalarTransport_d.*`
- RANS 適用層 (k/ω への共通コア適用・k/ω 勾配): `ransTransport_d.*`
- model-specific boundary: `ransBoundary_d.*` など
- model-specific source / closure: `ransSource_d.*` など

`scalarTransport_d.*` は `ScalarTransportDesc` を受ける汎用 launch ヘルパ
(`scalarTransportResidual_d` / `scalarTimeIntegration_d`) のみを公開し、
k/ω 固有の descriptor 構築・有効化判定・勾配計算は `ransTransport_d.*` 側に置く。
将来 species など別の scalar を足す場合も同じ共通コアを再利用する。

初回は advection-only とし、descriptor は内部実装用の metadata として扱う。
ユーザ入力に `enabled` を追加するのではなく、`LESorRANS == 2 && RANSmodel == 1`
のとき `k`, `omega` を自動有効化する。

一方で source は core に埋め込まず、main loop で `ransTransport_d_wrapper()` の後に
model-specific source wrapper を並べる。

diffusion では `k` と `omega` に別々の係数を与え、`vis_lam` と `vis_turb`
から face 係数を組み立てる。

## 5. 軸対称 SST の幾何項

軸対称 SST の幾何項は子 plan
[`architecture-axisym-sst.md`](../../.github/plans/architecture-axisym-sst.md)
で扱う。理論は [theory.md §7](theory.md) を参照。実装上の要点は次のとおり。

### 5.1 対流・拡散

対流・拡散は B 流儀の $r$ 重み付き面積 (`A_planar`, `sx_planar`) に幾何効果が
畳み込まれており、追加 source は不要 (theory.md §7.1)。
`ransGradient_d_wrapper` と `scalar_diffusion_first_order_d` は
すでに `A_planar` / `sx_planar` を使うため、コード変更なしで軸対称に対応する。
拡散の $\theta\theta$ 寄与が本当に発散形を保つかは、子 plan で
半径方向 1D 解析解との単体比較で確認する。

### 5.2 生産項・渦粘性のフープひずみ

`rans_sst_source_d` ([`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu))
と `sst_eddy_viscosity_d`
([`turbulent_viscosity_d.cu`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu))
は現在 planar 速度勾配 9 成分のみから $S^2$ を組む。軸対称では
フープひずみ $S_{\theta\theta} = u_r/r$ が欠落するため (theory.md §7.2)、次を追加する。

- 両 kernel に `int isAxisymmetric` と `flow_float* axisym_uy_over_r` を引数追加。
- `isAxisymmetric == 1` のとき
  $S^2 \mathrel{+}= 2\,(\texttt{axisym\_uy\_over\_r}[ic])^2$ を加える。
- `axisym_uy_over_r` は NS の
  [`axisymmetricSource_d.cu`](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)
  が毎ステップ更新済みなので、SST 側で再計算しない。呼び出し順序は
  `axisymmetricDiagnostics_d` → SST source の順を `main.cpp` で保証する。

### 5.3 圧縮性 (dilatation) 補正

生産項の圧縮性補正は `solverConfig` の `dilatationCorrection` で段階制御する
(theory.md §7.3)。

| 値 | 内容 |
| --- | --- |
| `0` | off (非圧縮形 $P_k = \mu_t S^2$) |
| `1` | (A) deviatoric のみ: $S^2 \mathrel{-}= \tfrac23(\nabla\!\cdot\!\mathbf u)^2$ |
| `2` | (A)+(B): さらに $P_k \mathrel{-}= \tfrac23\rho k(\nabla\!\cdot\!\mathbf u)$ を加え $P_k\ge0$ クリップ (**既定値**) |

既定は `2`。SST kernel は `LESorRANS==2 && RANSmodel==1` でのみ走るため、
laminar / LES ケースには影響しない。従来の非圧縮挙動が必要な場合のみ明示的に `0` を指定する。

実装は `rans_sst_source_d`
([`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu)) 内で完結する。

- 発散は `divU = dUxdx + dUydy + dUzdz`、軸対称時は
  さらに `+= axisym_uy_over_r[ic]`（= 完全発散 `axisym_divU` と一致）。
  既存配列を流用し新規計算しない。
- (A): $P_k, P_\omega$ に使う $S^2$ から $\tfrac23(\nabla\!\cdot\!\mathbf u)^2$ を引く。
- (B): $k$ 方程式の生産にのみ $-\tfrac23\rho k(\nabla\!\cdot\!\mathbf u)$ を加え、
  生産リミッタ適用後に `Pk = max(Pk, 0)`。$\omega$ 生産には入れない。
- `dilatationCorrection == 0` では従来挙動とビット一致。

検証は run_0091(off) を基準に run_0092(A) → run_0093(A+B) と段階比較する。

### 5.3.1 Kato–Launder 生産補正 (`katoLaunder`)

strain-based 生産の stagnation/加速アノマリー (theory.md §7.5) を抑えるオプション。
`solverConfig` の `turbulence.katoLaunder` で制御する。

| 値 | 内容 |
| --- | --- |
| `0` | 標準 $P_k=\mu_t S^2$ (**既定**, 既存挙動とビット一致) |
| `1` | Kato–Launder $P_k=\mu_t S\,\Omega$ ($\Omega=\sqrt{2\Omega_{ij}\Omega_{ij}}=|\boldsymbol\omega|$) |

実装も `rans_sst_source_d`
([`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu)) 内で完結:

- 渦度の大きさ二乗を既存の速度勾配から組む (新規物理量なし):
  `Om_sq = (dUxdy-dUydx)^2 + (dUxdz-dUzdx)^2 + (dUydz-dUzdy)^2` $= 2\Omega_{ij}\Omega_{ij}=|\boldsymbol\omega|^2$。
  軸対称(無旋回)でもフープは渦度を持たないため**補正不要**(ひずみと非対称)。
- `dilatationCorrection` 適用後の `S_sq` を用い、生産係数を
  $S^2 \to \sqrt{S\_sq}\cdot\sqrt{Om\_sq}$ に置換 (`katoLaunder==1` のとき)。
  $P_k=\mu_t\cdot S\Omega$、$P_\omega=\alpha\rho\cdot S\Omega$ の双方に適用。
- リミタ ($10\beta^*\rho k\omega$) と (B) 等方項は従来どおり後段で適用。
- `katoLaunder == 0` では従来挙動とビット一致。`dilatationCorrection` と直交・併用可。

**検証結果** (case 29 conical 2D, `run_0017`(off) vs `run_0018_conical_kl`(on)):
`katoLaunder:1` で軸中心線 $k$ **17.0→1.93** (μt/μlam 37→10.7), 核 $k$ **5.0→0.47** と偽生産を除去。
一方 **壁 BL は不変** (x=50mm で μt/μlam ピーク 13.0→12.9, $k$ 5374→5369), **推力も不変**
(F=1953.0→1953.1, mdot・λ 一致)。$S\approx\Omega$ の BL を保ち $S\gg\Omega$ の核生産のみ落とす狙い通り。

### 5.4 非軸対称ケースの不変性

`isAxisymmetric == 0` のときは追加項をすべて 0 とし、既存の 2D / 3D SST と
ビット単位で同一の挙動を維持する (子 plan の回帰判定基準)。

## 6. 陰解法との連成 (segregated point-implicit, 2026-06)

平均流 5 式は block DPLUR (5×5)、SST (k/ω) は **segregated point-implicit** で別に解く
（7×7 連成はしない）。block 陰解法 (`timeIntegration==11`) では平均流を commit した後、
同じ擬似時間ステップで k/ω を point-implicit 更新する（旧実装の「凍結」を解除）。

- 消散項のみ陰化（stiff、$\partial D_k/\partial(\rho k)=\beta^\*\omega$, $\partial D_\omega/\partial(\rho\omega)=2\beta\omega$）。
  `ransSource_d.cu` が残差に加えて `src_jac_k`/`src_jac_omega` を出力。
- `update_d.cu` の `applySSTPointImplicit_d` が $D_\phi=V/\Delta\tau+V\,\partial D/\partial(\rho\phi)$、
  $\delta(\rho\phi)=\text{relax}\cdot\text{res}_{\rho\phi}/D_\phi$ を作り $\rho\phi=\max(\rho\phi^N+\delta,\text{floor})$ を適用。
- `main.cpp` `implicitNonlinearUpdate` で `scalarResidualEnabled` のとき `applyBlockImplicitCorrection` 直後に呼ぶ。
- 移流・拡散・生産は残差 lagged。詳細は [`docs/time_integration/`](../time_integration/implementation.md)。

将来、強連成が必要なら `gpu-implicit-plan` の子 plan として 7×7 化を検討する。

## 7. 検証方針

標準の検証起点は [../../.github/forge-verification-cases.md](../../.github/forge-verification-cases.md)
に従う。

推奨順序:

1. `case/23.axi_nozzle/run_0033_slau_axisym_m4_publicrao_full_10k` の複製 run で advection-only の通りを確認
2. laminar / LES case で null regression
3. その後、airfoil 系で壁面境界層と `vis_turb` の立ち上がり確認
4. 最後に 3D case と軸対称子 plan へ展開
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
`ransScalarBoundary_d.*` と `ransSource_d.*` のような model-specific layer に分離する。

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
- turbulence scalar BC はこの後段に差し込む `ransScalarBoundary_d.*` へ移していく
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

## 4. Generic Scalar Transport 基盤

`k`, `omega` の輸送は、将来の passive scalar や species でも再利用できる
generic scalar transport 層で扱う。

実装上の責務分担は次とする。

- NS 本体: `convectiveFlux_d.*`, `viscousFlux_d.*`
- scalar 共通: `scalarTransport_d.*`
- model-specific boundary: `ransScalarBoundary_d.*` など
- model-specific source / closure: `ransSource_d.*` など

初回は advection-only とし、descriptor は内部実装用の metadata として扱う。
ユーザ入力に `enabled` を追加するのではなく、`LESorRANS == 2 && RANSmodel == 1`
のとき `k`, `omega` を自動有効化する。

一方で source は core に埋め込まず、main loop で `scalarTransport_d_wrapper()` の後に
model-specific source wrapper を並べる。

diffusion では `k` と `omega` に別々の係数を与え、`vis_lam` と `vis_turb`
から face 係数を組み立てる。

## 5. 軸対称との切り分け

軸対称は本実装文書の初期スコープ外とする。理由は次のとおり。

- 既存 plan で RANS が別 plan と明示されている
- B 流儀の幾何 source を scalar 輸送へどう入れるかを整理する必要がある
- 軸近傍での $\omega$ 境界・source の剛性確認が別途必要

したがって、親 plan では 2D / 3D explicit SST を完了条件とし、
軸対称は子 plan に委ねる。

## 6. 陰解法との切り分け

`timeIntegration_d.cu` とその周辺は 5 方程式固定の correction / Jacobian を
持っているため、SST を同時に陰解法へ入れると設計面が大きく広がる。

そのため初期実装では:

- explicit storage と residual 更新だけを先に広げる
- implicit block system は不変のままにする
- 将来必要なら `gpu-implicit-plan` の子 plan として 7x7 化を扱う

## 7. 検証方針

標準の検証起点は [../../.github/forge-verification-cases.md](../../.github/forge-verification-cases.md)
に従う。

推奨順序:

1. `case/23.axi_nozzle/run_0033_slau_axisym_m4_publicrao_full_10k` の複製 run で advection-only の通りを確認
2. laminar / LES case で null regression
3. その後、airfoil 系で壁面境界層と `vis_turb` の立ち上がり確認
4. 最後に 3D case と軸対称子 plan へ展開
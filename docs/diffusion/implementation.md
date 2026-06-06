# 粘性・熱伝導フラックス — 実装

forge の粘性・熱伝導フラックス計算の実装と、ソース上の対応関係。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/cuda_forge/viscousFlux_d.cuh`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cuh) | カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) | 内部面・壁面カーネル本体 (約 420 行) |

CPU 経路の独立実装は無く、GPU 経路のみで運用される。

## エントリポイント

[`viscousFlux_d_wrapper`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L315) が共通入口。
順に次を呼ぶ。

1. `viscousFlux_d` — 内部面 (`nNormalPlanes`) を並列処理。
2. `viscousFlux_wall_d` — 各壁面境界条件について壁面寄与を計算。

`res_*` バッファは対流フラックス計算後に粘性が追加加算される (ゼロクリアは対流側で実施)。

## `viscousFlux_d` 構造 ([L3](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L3))

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
8. 熱伝導束 `heatflux = thermCond*(T1-T0)/dcc * delta + thermCond*(dT*df * k)`。
9. 残差 `res_roUx, res_roUy, res_roUz, res_roe` を両側に符号反転で `atomicAdd`。

エネルギ残差には `tau_x*Uxf + tau_y*Uyf + tau_z*Uzf` (応力仕事) と `heatflux` を加える。

## `viscousFlux_wall_d` ([L150](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L150))

各壁面境界に対し、ゴーストセル `ig` と内部セル `ic` 間の片側差分で粘性応力を構築。
ラッパは `for (auto& bc : msh.bconds)` ループでこのカーネルを発行する。

## 入出力

入力: セル中心保存量・原始量 (`ro, roUx, …, Ts`)、勾配 (`dUxd*, dUyd*, dUzd*, dTd*`)、
粘性係数 (`vis_lam`, `vis_turb`)、`mu`, `thermCond` (定数)。

出力 (累積): `res_ro, res_roUx, res_roUy, res_roUz, res_roe`。

## 既知の TODO / 注意点

- 法線差分に用いるのはセル中心値のみで、再構成 (MUSCL) は適用しない (拡散項として妥当)。
- 壁面カーネルは現状 `wall` / `wall_isothermal` のみ呼ばれる。`slip` には粘性束を加えない。
- `vis_lam` は配列で渡される構造になっているが、現状はラッパで `cfg.visc` を直接渡している。
  温度依存粘性 (Sutherland 法則など) を導入する場合は配列を実体化する。

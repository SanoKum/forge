# 境界条件 — 実装

forge の境界条件 (BC) 実装とソース対応をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/boundaryCond.hpp`](../../solver_density_cuda/boundaryCond.hpp) | `bcond`, `bcondConfFormat` 構造体、各 BC kind の値型テーブル |
| [`solver_density_cuda/boundaryCond.cpp`](../../solver_density_cuda/boundaryCond.cpp) | YAML 読み込みと `applyBconds` ディスパッチ |
| [`solver_density_cuda/cuda_forge/boundaryCond_d.cuh`](../../solver_density_cuda/cuda_forge/boundaryCond_d.cuh) | GPU カーネル / ラッパ宣言 |
| [`solver_density_cuda/cuda_forge/boundaryCond_d.cu`](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu) | 全 BC kind のカーネル本体 (約 1500 行) |

CPU 経路 (`cfg.gpu == 0`) は未対応 (`applyBconds` 冒頭で `exit`)。

## YAML 設定

ケース直下の `bcondConfig.yaml` で各物理 ID の BC を定義する。

```yaml
1:
  physID: 1
  kind: wall
  outputHDFflg: 1
  ints: {}
  floats: { Ts: 300.0 }
```

`boundaryCond.cpp` の [`readBcondConfig`](../../solver_density_cuda/boundaryCond.cpp#L11) が
YAML を読み、各 `bcond` (メッシュ生成側で物理 ID 付き) に `kind`, パラメータ,
値型テーブル (`valueTypes`) を紐付け、`bcondInitVariables` で
GPU/CPU 双方の境界変数バッファを確保する。

## ディスパッチ

[`applyBconds`](../../solver_density_cuda/boundaryCond.cpp#L116) が
`msh.bconds` をループし、各 `bc.bcondKind` で次のラッパを呼ぶ。

| `bcondKind` | ラッパ | カーネル位置 |
| --- | --- | --- |
| `slip` | `slip_d_wrapper` | [L86](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L86) |
| `wall` | `wall_d_wrapper` | [L214](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L214) |
| `wall_isothermal` | `wall_isothermal_d_wrapper` | [L363](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L363) |
| `outlet_statPress` | `outlet_statPress_d_wrapper` | [L531](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L531) |
| `inlet_uniformVelocity` | `inlet_uniformVelocity_d_wrapper` | [L651](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L651) |
| `inlet_fluctVelocity` | (`fluct_variables_d` 経由) | — |
| `inlet_Pressure` | `inlet_Pressure_d_wrapper` | [L863](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L863) |
| `inlet_Pressure_dir` | `inlet_Pressure_dir_d_wrapper` | [L1033](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1033) |
| `outflow` | `outflow_d_wrapper` | [L1184](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1184) |
| `periodic` | `periodic_d_wrapper` | [L1310](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1310) |

呼び出し後に `cudaPeekAtLastError + cudaDeviceSynchronize` で同期。

### `outlet_statPress` の特性ベース構成と逆流統一

`outlet_statPress_d` は亜音速流出で指定静圧 `Psb[ib]` のみを境界値とし、$\rho$・速度は内部セル `ic` の
エントロピーと外向き Riemann 不変量から構成する (理論は [theory.md](theory.md) 参照)。**逆流 (`Un<0`) も
同じ静圧アンカーで扱う** (`Vn_exit<0` 許容)。値の取得元に注意: 静圧のみ境界 `Psb[ib]`、エントロピー・
Riemann・接線速度はすべて内部 `ic`、面法線は plane の `sx/sy/sz`。乱流スカラーは `rans_neumann_scalar_boundary_d`
でゼロ勾配 (`k[ig]=k[ic]`)。

旧実装は逆流時に全圧 `Ptb/Ttb` の stagnation 流入 (=`inlet_Pressure` 構成) へ切替え、かつ速度を
`-Ux[ic]*nx` で構成していたため、剥離 BL が出口に達するケース (擬似衝撃ダクト, 壁∩出口コーナー) で
過加圧・forward↔backflow のばたつきにより発散していた。検証は
[`.github/plans/boundary-outlet-characteristic-backflow.md`](../../design/accepted/boundary-outlet-characteristic-backflow.md)。

## 入力データ構造

各 `bcond` (境界グループ) は次の GPU データを持つ。

- `map_bplane_plane_d[ib]` — グループ内インデックス `ib` → メッシュ面 ID `ip`
- `map_bplane_cell_d[ib]` — `ip` の内部セル ID
- `map_bplane_cell_ghst_d[ib]` — 対応するゴーストセル ID (`>= nCells`)
- `bvar_d[*]` — 境界面値 (`ro`, `roUx`, …, `Ts`, `Psb` ほか)
- `inputInts`, `inputFloats` — YAML から渡された定数

## カーネル構造 (例: `slip_d`)

[`slip_d`](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L3) は境界面 1 並列。

```cpp
ip = bplane_plane[ib];
ic = bplane_cell[ib];
ig = bplane_cell_ghst[ib];

Un = (sx[ip]*Ux[ic] + sy[ip]*Uy[ic] + sz[ip]*Uz[ic]) / ss[ip];

// ghost cell  ← 法線速度を反転
ro[ig]    = ro[ic];
Ux[ig]    = Ux[ic] - 2*Un*sx[ip]/ss[ip];
P[ig]     = P[ic];
roe[ig]   = P[ic]/(ga-1) + 0.5*ro[ic]*|U_L|^2;
sonic[ig] = sqrt(ga*P[ig]/ro[ig]);

// boundary value  ← 接線速度のみ (Un を 1 回引く)
Uxb[ib] = Ux[ic] - Un*sx[ip]/ss[ip];
…
```

他の BC kind も同じパターンで、`ig` (内部面で R 状態として参照) と `ib`
(境界面値で粘性束等から参照) の双方を書き込む。

## copyBcondsGradient (現状無効)

[`copyBcondsGradient_d_wrapper`](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu#L1439)
は境界面の勾配をゴーストセルにコピーする補助カーネルだが、現状の `boundaryCond.cpp`
からは呼び出されていない (コメントアウト)。境界勾配を厳密に扱う実装に
切り替える場合に有効化する。

## 既知の TODO / 注意点

- CPU 経路 (`cfg.gpu == 0`) は未対応。
- `inlet_fluctVelocity` は `fluct_variables_d` モジュールの変動生成と組み合わせて使う。
- 周期境界は対流再構成で 1 次風上強制が外れる (`ic1 < nCells` を維持)。
  ペアリングはメッシュ生成側で行う。
- `bcondConfig.yaml` の `kind` を新設する場合は、`boundaryCond.hpp` の
  `valueTypesOfBC` に値型テーブルを追加し、`applyBconds` にディスパッチ分岐を増やす。

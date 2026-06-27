# 乱流熱伝導 (turbulent thermal conductivity) の追加

## メタ

- **area**: `diffusion`
- **status**: `done`
- **related_docs**:
  - `methods/diffusion/theory.md`
  - `methods/diffusion/implementation.md`
- **created**: `2026-06-14`
- **owner**: `CFD Dev`

## 1. 目的

RANS (k-ω SST 等) のエネルギー方程式に乱流熱伝導 $\kappa_t = c_p\,\mu_t/Pr_t$ が
欠落しており、乱流境界層で散逸熱 (摩擦発熱) を運ぶ熱伝導が層流値のみだったため、
断熱壁の静温が回復温度 ($\le$ 全温) を大きく超えて overshoot していた。
有効熱伝導率に乱流寄与を加え、エネルギー保存を回復する。

## 2. スコープ

- **やる**: `viscousFlux_d.cu` の内部面・壁面カーネルの熱伝導束で
  $\kappa_{\text{eff}} = \kappa_{\text{lam}} + c_p\,\mu_t/Pr_t$ ($Pr_t=0.9$) を使う。
- **やらない**: 乱流 $Pr_t$ の config 化 (当面 0.9 固定)、$(2/3)\rho k$ 乱流法線応力の
  運動量/エネルギーへの厳密追加 (寄与が小さく別 plan)、可変 $c_p$ (thermalMethod==2) 対応。

## 3. 関連 docs と前提

- 応力は既に $\mu = \mu_{\text{lam}} + \mu_{\text{turb}}$ (`mu_total`) を使用 (`methods/diffusion/`)。
- `vis_turb[ic]` は SGS/RANS モデルが毎ステップ更新済み (`ransSource_d.cu`)。
- 熱伝導束は内部面 `viscousFlux_d` と壁面 `viscousFlux_wall_d` の2カーネルにある。

## 4. 設計方針

摩擦発熱は $\mu_{\text{total}}$ で評価されるのに対し、従来の熱伝導束は層流 `thermCond`
(定数) のみだった。乱流境界層では $\mu_t$ が層流の数十倍になり、対応する渦熱伝導
$\kappa_t = c_p\mu_t/Pr_t$ が層流 $\kappa$ の ~50 倍に達する。この熱を運ぶ項が無いと
散逸熱が局所に溜まり静温が全温を超える (本ケースで 449 K vs 全温 293 K)。

実装は熱伝導束で有効伝導率を使うだけ:

```cpp
const flow_float Prt = 0.9;
tc_face = f*thermCond[ic0] + (1-f)*thermCond[ic1];   // 層流
tc_face += cp*v_turb/Prt;                            // 乱流 (v_turb は面平均済み)
heatflux = tc_face*((Ts1-Ts0)/dcc)*delta + tc_face*(dT*df·k);
```

カーネルに `cp` (= `cfg.cp`) を引数追加。壁面カーネルも同形で `tc_w = thermCond + cp*vis_turb/Pr_t`
(断熱壁はミラーゴーストで $\Delta T=0$ のため寄与 0、isothermal 壁で有効)。

## 5. 実装ステップ

1. `viscousFlux_d` / `viscousFlux_wall_d` の signature に `flow_float cp` を追加。
2. 両カーネルの熱伝導束を $\kappa_{\text{eff}}$ 版に変更。
3. wrapper の2つの launch site で `cfg.cp` を渡す。
4. methods/diffusion の theory/implementation を更新。

## 変更ログ

- **2026-06-14**: 実装・検証完了。`case/16.nozzle_wys` の 2D SST (Wyslouzil ノズル,
  Pt=101325/Tt=293.15, 超音速流出) で検証。
  - 修正前: 静温 max **449 K** (全温 293.15 K を 156 K 超過)、`T>Tt` セル 7254 個 (近壁発散)。
  - 修正後: 静温 max **293.2 K** (全温と一致)、`T>Tt` セル 10 個 (丸め誤差レベル)。
  - 中心線 Mach 分布・残差収束は修正前後で大きく変わらず、過熱が近壁の乱流境界層に
    限定されていたことと整合。`run_0004_slau_2d_sst/` を根拠とする。

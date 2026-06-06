# 境界条件 — 理論

forge は密度ベース有限体積で **ゴーストセル方式** を採用し、各境界条件は
内部セル状態と境界条件パラメータからゴーストセル値を構成する。
内部面ループでは境界面もそのまま処理されるため、ゴースト値が境界条件を
そのまま表現する形になる。

## 共通モデル

非周期境界面 $f$ について内部セルを $L$、ゴーストセルを $R$ と書く。
境界条件は次の二系統を返す。

1. **ゴーストセル状態** $(\rho, \rho \mathbf u, \rho e, U, P, T, H_t, c)_R$ —
   対流フラックス計算で R 状態として使用される。
2. **境界面値 `bvar_d[*]`** — 粘性束計算・post 処理 (壁面摩擦・$y^+$) で参照される。

非周期境界では対流再構成が `scheme = -1` (1 次風上) に強制されるため、
ゴーストセル値はそのまま面値として使われる ([`docs/convection/`](../convection/) 参照)。

## 提供される境界種別

| `kind` | 用途 | ゴースト構成の概要 |
| --- | --- | --- |
| `slip` | 滑り壁 (オイラー壁) | 法線速度を反転 ($\mathbf u_R = \mathbf u_L - 2 U_n \hat{\mathbf n}$)。圧力・密度は同値 |
| `wall` | 非滑り断熱壁 | $\mathbf u_R = -\mathbf u_L$、$P_R = P_L$、$T_R = T_L$ (断熱) |
| `wall_isothermal` | 非滑り等温壁 | $\mathbf u_R = -\mathbf u_L$、$T_R = 2 T_{\text{wall}} - T_L$ |
| `inlet_uniformVelocity` | 均一速度流入 | $\mathbf u_R, \rho_R$ を指定値に固定、$P_R = P_L$ |
| `inlet_fluctVelocity` | 速度変動つき流入 | uniformVelocity に変動成分を加算 (`fluct_variables`) |
| `outlet_statPress` | 静圧固定流出 | $P_R = P_{\text{back}}$、その他は外挿 |
| `inlet_Pressure` | 全圧・全温固定流入 | 全条件 ($P_t, T_t$) から内部マッハで $P, T$ を再構成 |
| `inlet_Pressure_dir` | 方向指定全圧流入 | inlet_Pressure に流入方向ベクトルを併用 |
| `outflow` | サブソニック流出 | リーマン不変量に基づく Non-reflecting 流出 |
| `periodic` | 周期境界 | 対応するペア面のセル値をコピー (`scheme` 強制なし) |

## 例: 滑り壁

法線方向速度 $U_n = \mathbf u_L \cdot \hat{\mathbf n}$ を用いて

$$
\mathbf u_R = \mathbf u_L - 2 U_n \hat{\mathbf n},
\quad P_R = P_L, \quad \rho_R = \rho_L,
\quad \rho e_R = \frac{P_L}{\gamma - 1} + \tfrac{1}{2}\rho_L |\mathbf u_L|^2.
$$

これは法線フラックスから運動量の壁面垂直成分が消え、圧力束のみ残る効果を持つ。

## 例: 非滑り壁

$\mathbf u_R = -\mathbf u_L$ とすることで、面値 $\mathbf u_f = \tfrac{1}{2}(\mathbf u_L + \mathbf u_R) = 0$
(無滑り)。粘性束は壁面用カーネル ([`docs/diffusion/`](../diffusion/)) で別途加算する。

## 例: 全圧固定流入

総温・総圧 $T_t, P_t$ と局所マッハ $M_L$ を用い、

$$
T_R = \frac{T_t}{1 + \tfrac{\gamma - 1}{2} M_L^2},\quad
P_R = P_t \left(\frac{T_R}{T_t}\right)^{\gamma/(\gamma-1)},\quad
\rho_R = \frac{P_R}{R T_R},
$$

速度方向は外挿 (`inlet_Pressure`) または指定方向 (`inlet_Pressure_dir`)。

## 周期境界

周期境界は他境界と異なり、対応するペアセル値を直接コピーする。
対流再構成は内部面と同じ MUSCL 経路で処理されるため、`scheme` の強制 1 次降格は無い。

## 物理 ID と YAML 設定

各境界面はメッシュ生成時に物理 ID を持ち、`bcondConfig.yaml` で
ID → `kind` → 数値パラメータ (`ints`, `floats`) を紐付ける。
パラメータの種類 (定数 / 配列入力) は `bcondConfFormat::valueTypesOfBC` で定義する。

## 参考

- 適用順序: [`docs/time_integration/`](../time_integration/) の "ループ全体" を参照。
- 粘性壁面寄与: [`docs/diffusion/implementation.md`](../diffusion/implementation.md) の
  `viscousFlux_wall_d` を参照。
- 対流側の境界処理 (ゴースト → 1 次風上): [`docs/convection/`](../convection/) を参照。

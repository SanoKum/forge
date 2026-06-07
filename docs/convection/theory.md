# 対流フラックス — 理論

forge は圧縮性 Navier–Stokes 方程式の対流項を有限体積法のセル境界フラックスとして
評価する。本ファイルではセル両側の状態量再構成 (MUSCL 系) と、フラックス分割
(SLAU / Roe / HLLE / KEEP / AUSM 系) の理論を整理する。
実装の対応は [implementation.md](implementation.md) を参照。

## 共通フレームワーク

保存変数 $\mathbf{Q} = (\rho,\, \rho u,\, \rho v,\, \rho w,\, \rho e)^\top$
の対流フラックスをセル境界 $f$ で評価し、隣接 2 セル $L, R$ に符号反転して加算する。

$$
\frac{d \mathbf{Q}_c}{dt} V_c
= -\sum_{f \in \partial c} \mathbf{F}^{\text{conv}}_f \cdot \mathbf{S}_f
\;+\; (\text{粘性・湧出項}).
$$

各面で次の手順を踏む。

1. **L/R 状態の再構成** — セル中心値とセル中心勾配 $\nabla \phi$ から面側で
   $\phi_L, \phi_R$ を作る (MUSCL 第 2 / 第 3 次)。
2. **リミッタ** — Ducros センサ・スカラリミッタにより不連続近傍で次数を落とす。
3. **数値フラックス** — 上記 L/R 状態から $\mathbf{F}^{\text{conv}}_f$ を組む。

## 状態再構成

セル $C$ から面側へ向かう線形再構成は

$$
\phi_L = \phi_C + \psi_C \,\nabla \phi_C \cdot (\mathbf{r}_f - \mathbf{r}_C),
$$

$\psi_C$ がリミッタ (1 で無リミット)。第 3 次 MUSCL では隣接セル $D$ の差分を加える。

$$
\phi_L = \phi_C + \psi_C \left[\,\tfrac{1}{2} k\, (\phi_D - \phi_C) +
(1-k)\, \nabla \phi_C \cdot (\mathbf{r}_f - \mathbf{r}_C)\right],
\qquad k = \tfrac{1}{3}.
$$

$\phi_R$ は対称に、$D \to C$ を反転して構成する。
非周期境界面では右側 $R$ がゴーストセルになり、ゴーストには勾配・リミッタが
存在しない。このため境界面では再構成スキームを 1 次風上に落とす
(`scheme = -1`)。

### Ducros センサ

衝撃波近傍で次数を落とすため、Ducros センサ

$$
\mathrm{duc} = \frac{(\nabla \cdot \mathbf{u})^2}
{(\nabla \cdot \mathbf{u})^2 + \|\nabla \times \mathbf{u}\|^2 + \varepsilon}
$$

を用い、$\mathrm{duc} > 0.8$ の領域で $\psi \leftarrow (1 - \mathrm{duc}) \psi$ と
減衰させる。

## 数値フラックススキーム

forge には以下が実装されている。括弧内は現状のラッパで有効化されているか。

| スキーム | 種別 | 状態 |
| --- | --- | --- |
| **Roe (FDS)** | 線形化 Riemann 解 | 有効 |
| **SLAU** | AUSM 系 (低マッハ性) | 有効 |
| **HLLE** | HLL 系 (中央波スキップ) | 有効 |
| KEEP_FVS | 中心型 (運動エネルギ保存) + FVS 散逸 | コード有・ラッパで無効 |
| KEEP_SLAU | KEEP + SLAU 散逸 | コード有・ラッパで無効 |
| AUSM+ | AUSM 系 | コード有・ラッパで無効 |
| AUSM+UP | AUSM 系 (低マッハ補正) | コード有・ラッパで無効 |

### Roe (FDS)

#### 基本式

中央フラックスから Roe 平均ヤコビアンに基づく散逸項を引いた近似 Riemann 解。
forge では運動量・エネルギ束を質量流束 $\dot m$ で書く形式に整理しているが、
等価に書き直すと

$$
\mathbf{F}^{\text{Roe}}_f \cdot \mathbf{S}
= \tfrac{1}{2}(\mathbf{F}_L + \mathbf{F}_R)\cdot\mathbf{S}
- \tfrac{1}{2}\, S\,\widetilde{R}\,|\widetilde{\Lambda}|\,\widetilde{L}\,(\mathbf{Q}_R - \mathbf{Q}_L),
\qquad S = |\mathbf{S}|.
$$

#### Roe 平均状態

$$
\widetilde{\rho} = \sqrt{\rho_L \rho_R},\qquad
\widetilde{\phi} = \frac{\sqrt{\rho_L}\,\phi_L + \sqrt{\rho_R}\,\phi_R}{\sqrt{\rho_L} + \sqrt{\rho_R}}
\quad (\phi \in \{u, v, w, H_t\}),
$$

$$
\widetilde{c} = \sqrt{(\gamma-1)\big(\widetilde{H}_t - \tfrac{1}{2}|\widetilde{\mathbf u}|^2\big)},
\qquad
\widetilde{U} = \widetilde{u} n_x + \widetilde{v} n_y + \widetilde{w} n_z .
$$

#### 状態量とフラックスのベクトル順序 (forge 実装規約)

forge は保存量・固有ベクトル・差分ベクトルすべてで

$$
\mathbf{Q} = (\rho,\; \rho E,\; \rho u,\; \rho v,\; \rho w)^\top
$$

の順序を採る (一般的教科書の $(\rho, \rho u, \rho v, \rho w, \rho E)$ とは
エネルギ位置が違うので注意)。固有値の並びは

$$
\widetilde{\Lambda} = \mathrm{diag}\big(|\widetilde{U} - \widetilde{c}|,\;
|\widetilde{U}|,\; |\widetilde{U}|,\; |\widetilde{U}|,\;
|\widetilde{U} + \widetilde{c}|\big),
$$

すなわち列 / 行は (遅い音波, せん断 ×3, 速い音波) の順。

#### 補助量

理想気体の圧力 $P = (\gamma-1)\big[\rho E - \tfrac{1}{2}\rho|\mathbf u|^2\big]$ を
保存変数で微分した量を定義する。

$$
\begin{aligned}
P_\rho   &= \tfrac{1}{2}(\gamma-1)|\widetilde{\mathbf u}|^2, &
P_{\rho E} &= \gamma - 1, \\
P_{\rho u} &= -(\gamma-1)\widetilde{u}, &
P_{\rho v} &= -(\gamma-1)\widetilde{v}, &
P_{\rho w} &= -(\gamma-1)\widetilde{w}.
\end{aligned}
$$

理想気体では $z \equiv |\widetilde{\mathbf u}|^2 - P_\rho/P_{\rho E} = \tfrac{1}{2}|\widetilde{\mathbf u}|^2$,
$\varphi \equiv P_\rho - \widetilde{c}^2$。

#### 右固有ベクトル $\widetilde{R}$

列 1 = $\widetilde{U}-\widetilde{c}$, 列 2–4 = $\widetilde{U}$, 列 5 = $\widetilde{U}+\widetilde{c}$。
行は $\mathbf{Q}$ 順 $(\rho, \rho E, \rho u, \rho v, \rho w)$。

$$
\widetilde{R} =
\begin{pmatrix}
1 & n_x & n_y & n_z & 1 \\[2pt]
\widetilde{H}_t - \widetilde{c}\widetilde{U}
 & z n_x + \widetilde{v} n_z - \widetilde{w} n_y
 & z n_y + \widetilde{w} n_x - \widetilde{u} n_z
 & z n_z + \widetilde{u} n_y - \widetilde{v} n_x
 & \widetilde{H}_t + \widetilde{c}\widetilde{U} \\[2pt]
\widetilde{u} - \widetilde{c} n_x
 & \widetilde{u} n_x
 & \widetilde{u} n_y - n_z
 & \widetilde{u} n_z + n_y
 & \widetilde{u} + \widetilde{c} n_x \\[2pt]
\widetilde{v} - \widetilde{c} n_y
 & \widetilde{v} n_x + n_z
 & \widetilde{v} n_y
 & \widetilde{v} n_z - n_x
 & \widetilde{v} + \widetilde{c} n_y \\[2pt]
\widetilde{w} - \widetilde{c} n_z
 & \widetilde{w} n_x - n_y
 & \widetilde{w} n_y + n_x
 & \widetilde{w} n_z
 & \widetilde{w} + \widetilde{c} n_z
\end{pmatrix}.
$$

#### 左固有ベクトル $\widetilde{L} = \widetilde{R}^{-1}$

実装では各行を一度ベタ書きし、最後にまとめて $1/\widetilde{c}^2$ で割る。

$$
\widetilde{L}_{1*} = \tfrac{1}{2\widetilde{c}^2}\!\left(
P_\rho + \widetilde{c}\widetilde{U},\;
P_{\rho E},\;
P_{\rho u} - \widetilde{c} n_x,\;
P_{\rho v} - \widetilde{c} n_y,\;
P_{\rho w} - \widetilde{c} n_z
\right),
$$

$$
\widetilde{L}_{2*} = \tfrac{1}{\widetilde{c}^2}\!\left(
-\varphi n_x - \widetilde{c}^2(\widetilde{v} n_z - \widetilde{w} n_y),\;
-P_{\rho E} n_x,\;
-P_{\rho u} n_x,\;
-P_{\rho v} n_x + \widetilde{c}^2 n_z,\;
-P_{\rho w} n_x - \widetilde{c}^2 n_y
\right),
$$

$$
\widetilde{L}_{3*} = \tfrac{1}{\widetilde{c}^2}\!\left(
-\varphi n_y - \widetilde{c}^2(\widetilde{w} n_x - \widetilde{u} n_z),\;
-P_{\rho E} n_y,\;
-P_{\rho u} n_y - \widetilde{c}^2 n_z,\;
-P_{\rho v} n_y,\;
-P_{\rho w} n_y + \widetilde{c}^2 n_x
\right),
$$

$$
\widetilde{L}_{4*} = \tfrac{1}{\widetilde{c}^2}\!\left(
-\varphi n_z - \widetilde{c}^2(\widetilde{u} n_y - \widetilde{v} n_x),\;
-P_{\rho E} n_z,\;
-P_{\rho u} n_z + \widetilde{c}^2 n_y,\;
-P_{\rho v} n_z - \widetilde{c}^2 n_x,\;
-P_{\rho w} n_z
\right),
$$

$$
\widetilde{L}_{5*} = \tfrac{1}{2\widetilde{c}^2}\!\left(
P_\rho - \widetilde{c}\widetilde{U},\;
P_{\rho E},\;
P_{\rho u} + \widetilde{c} n_x,\;
P_{\rho v} + \widetilde{c} n_y,\;
P_{\rho w} + \widetilde{c} n_z
\right).
$$

#### 散逸構成

$\Delta \mathbf{Q} = \mathbf{Q}_R - \mathbf{Q}_L$ を $\widetilde{L}$ で特性座標 $\Delta \mathbf{W}$ に
変換 → 固有値の絶対値を掛けて → $\widetilde{R}$ で物理座標に戻す。

$$
\Delta \mathbf{W} = \widetilde{L}\, \Delta \mathbf{Q},\qquad
\Delta \mathbf{F}^{\text{diss}} = \widetilde{R}\,(|\widetilde{\Lambda}|\,\Delta \mathbf{W}).
$$

forge の最終形 (面積 $S$ 込みで残差に積算)。質量・エネルギ・運動量の順:

$$
\begin{aligned}
\mathrm{Res}_{\rho}    &= \dot m - \tfrac{S}{2}\,\Delta F^{\text{diss}}_1, \\
\mathrm{Res}_{\rho E}  &= \dot m\,\tfrac{1}{2}(H_{t,L}+H_{t,R}) - \tfrac{S}{2}\,\Delta F^{\text{diss}}_2, \\
\mathrm{Res}_{\rho u}  &= \dot m\,\tfrac{1}{2}(u_L+u_R) + \tfrac{1}{2}(P_L+P_R) S_x - \tfrac{S}{2}\,\Delta F^{\text{diss}}_3, \\
\mathrm{Res}_{\rho v}  &= \dot m\,\tfrac{1}{2}(v_L+v_R) + \tfrac{1}{2}(P_L+P_R) S_y - \tfrac{S}{2}\,\Delta F^{\text{diss}}_4, \\
\mathrm{Res}_{\rho w}  &= \dot m\,\tfrac{1}{2}(w_L+w_R) + \tfrac{1}{2}(P_L+P_R) S_z - \tfrac{S}{2}\,\Delta F^{\text{diss}}_5,
\end{aligned}
$$

$$
\dot m = \tfrac{1}{2}\big[\rho_L (U_L) + \rho_R (U_R)\big]\, S,
\qquad U_{L/R} = \mathbf{u}_{L/R}\cdot\hat{\mathbf n}.
$$

#### エントロピーフィックス (Harten 形)

$\widetilde{U} \approx \pm \widetilde{c}$ で発生する非物理な expansion shock を防ぐため、
小さな固有値を放物線で持ち上げる。forge では速度比依存型を採用:

$$
\eta = 0.1\!\left(\frac{|\widetilde{U}|}{\widetilde{c}} + 1\right),
\qquad
|\widetilde{\lambda}| \leftarrow
\begin{cases}
|\widetilde{\lambda}| & (|\widetilde{\lambda}| > \eta),\\[2pt]
\dfrac{|\widetilde{\lambda}|^2 + \eta^2}{2\eta} & (|\widetilde{\lambda}| \le \eta).
\end{cases}
$$

van Leer 型と Harten 定数 $\eta$ 型もコード上に残っているが、現状は速度比依存型のみ有効。

#### 小値クランプ

$\rho_{L/R} \ge 10^{-8}$, $P_{L/R} \ge 10^{-8}$, $\widetilde{c}^2 \ge 10^{-8}$ で
除算保護する。

### SLAU

Shima–Kitamura の SLAU (Simple Low-dissipation AUSM)。低マッハ域の数値散逸を
小さく抑え、衝撃波の遷音速域でも単調性を保つ。

#### 入力量

$$
U_{L/R} = \mathbf{u}_{L/R}\cdot\hat{\mathbf n},\qquad
\widehat{c} = \tfrac{1}{2}(c_L + c_R)\quad\text{(算術平均、Roe 平均ではない)},
$$

$$
M_+ = U_L/\widehat{c},\qquad M_- = U_R/\widehat{c}.
$$

#### Mach 関数

$$
\beta_+(M) = \begin{cases}
\tfrac{1}{2}(M + |M|)/M & (|M| \ge 1),\\[2pt]
\tfrac{1}{4}(M+1)^2(2 - M) & (|M| < 1),
\end{cases}
\quad
\beta_-(M) = \begin{cases}
\tfrac{1}{2}(M - |M|)/M & (|M| \ge 1),\\[2pt]
\tfrac{1}{4}(M-1)^2(2 + M) & (|M| < 1).
\end{cases}
$$

#### 低マッハ補正

衝撃波付近の質量散逸を抑える因子 $g$ と、超音速近傍で散逸を切り替える $\chi$:

$$
g = -\max(\min(M_+, 0), -1) \cdot \min(\max(M_-, 0), 1),
$$

$$
|\widehat{V}_n| = \frac{\rho_L |U_L| + \rho_R |U_R|}{\rho_L + \rho_R},\qquad
|\widehat{V}_n^+| = (1-g)|\widehat{V}_n| + g\,|U_L|,\qquad
|\widehat{V}_n^-| = (1-g)|\widehat{V}_n| + g\,|U_R|,
$$

$$
\widehat{M} = \min\!\Big(1,\, \frac{\sqrt{\tfrac{1}{2}(|\mathbf u_L|^2 + |\mathbf u_R|^2)}}{\widehat{c}}\Big),
\qquad \chi = (1 - \widehat{M})^2.
$$

#### 質量流束と圧力束

$$
\dot m = \tfrac{S}{2}\Big[\rho_L (U_L + |\widehat{V}_n^+|) + \rho_R (U_R - |\widehat{V}_n^-|)\Big]
- \tfrac{\chi\, S}{2\,\widehat{c}}\,(P_R - P_L),
$$

$$
\tilde{p} = \tfrac{1}{2}(P_L + P_R)
 + \tfrac{1}{2}(\beta_+ - \beta_-)(P_L - P_R)
 + (1 - \chi)(\beta_+ + \beta_- - 1)\, \tfrac{1}{2}(P_L + P_R).
$$

質量散逸の主項 $\chi(P_R - P_L)/\widehat{c}$ は低マッハ ($\chi \to 1$) で
有効に効き、超音速 ($\chi \to 0$) では消える。

#### 風上化と残差

$\dot m$ の符号でエンタルピと運動量を風上選択する:

$$
\begin{aligned}
\mathrm{Res}_{\rho}   &= \dot m, \\
\mathrm{Res}_{\rho u} &= \tfrac{1}{2}(\dot m + |\dot m|) u_L + \tfrac{1}{2}(\dot m - |\dot m|) u_R + \tilde{p}\, S_x, \\
\mathrm{Res}_{\rho v} &= \tfrac{1}{2}(\dot m + |\dot m|) v_L + \tfrac{1}{2}(\dot m - |\dot m|) v_R + \tilde{p}\, S_y, \\
\mathrm{Res}_{\rho w} &= \tfrac{1}{2}(\dot m + |\dot m|) w_L + \tfrac{1}{2}(\dot m - |\dot m|) w_R + \tilde{p}\, S_z, \\
\mathrm{Res}_{\rho E} &= \tfrac{1}{2}(\dot m + |\dot m|) h_{t,L} + \tfrac{1}{2}(\dot m - |\dot m|) h_{t,R},
\end{aligned}
$$

ここで $h_t = \gamma P / [(\gamma-1)\rho] + \tfrac{1}{2}|\mathbf u|^2 = H_t$ (全エンタルピ)。

#### SLAU2 (圧力束の低マッハ改良)

SLAU の圧力束 $\tilde p$ の**第 3 項のみ**を Kitamura–Shima (2013) の SLAU2 形に置き換える。
質量流束 $\dot m$・風上化・残差組み立ては SLAU と完全に同一:

$$
\tilde{p}_{\text{SLAU2}} = \tfrac{1}{2}(P_L + P_R)
 + \tfrac{1}{2}(\beta_+ - \beta_-)(P_L - P_R)
 + (\beta_+ + \beta_- - 1)\,\sqrt{\tfrac{1}{2}\big(|\mathbf u_L|^2 + |\mathbf u_R|^2\big)}\;\bar\rho\,\widehat{c},
$$

ここで $\bar\rho = \tfrac{1}{2}(\rho_L + \rho_R)$、$\widehat c$ は界面音速 (上と同じ算術平均)。

SLAU の第 3 項 $(1-\chi)(\beta_++\beta_--1)\tfrac{1}{2}(P_L+P_R)$ は超音速近傍で
$(1-\chi)$ を介して平均圧力に比例するのに対し、SLAU2 の項は
$\bar\rho\,\widehat c\,|\mathbf V| \;(\propto M\cdot\gamma P)$ に比例する。Kitamura–Shima (2013)
の本来の狙いは、この圧力束の再定式化による**衝撃波ロバスト性の改善**(カーバンクル等の
抑制) と $(1-\chi)$ 依存の除去であり、超音速域では両式は同等で解は実質不変。

**低マッハ残差フロアとの関係 (注意)**。低マッハの圧力–速度カップリング
(チェッカーボード抑制) を主に担うのは、SLAU/SLAU2 で**共通**の質量流束の圧力散逸項
$-\dfrac{\chi}{\widehat c}(P_R-P_L)$ である。SLAU2 が変えるのは圧力束 (運動量側) の第 3 項
だけなので、**低マッハのエネルギー残差フロア自体は SLAU2 では解消しない**。実際、3D ノズル
($M_{\text{chamber}}\sim0.06$, 陰解法) で SLAU と SLAU2 を同条件比較したところ、`rms_roe`
フロア・チャンバー圧力の非物理ばらつきはほぼ同一であった
([`convection-slau2-lowmach.md`](../../.github/plans/convection-slau2-lowmach.md) §9)。
低マッハフロアの根治には時間項の**前処理** (Weiss–Smith 等。質量流束の散逸スケールごと
是正する) が必要。`solver: SLAU2` は衝撃波ロバスト性向けの選択肢として有効 (既定は `SLAU`)。

これは時間項の低マッハ**前処理** (Weiss–Smith 等) とは別物で、空間 (フラックス) 側の
散逸スケールのみを是正する軽量な対処である。`solver: SLAU2` で選択する (既定は `SLAU`)。

### HLLE

Harten–Lax–van Leer–Einfeldt 形 HLL。中央波 (接触不連続) を捨て、最大・最小波速で
左右領域を区切る単一中間状態を取る。

#### 波速見積もり

Roe 平均 $\widetilde{U}, \widetilde{c}$ と左右の $U_{L/R}, c_{L/R}$ から

$$
S_L = \min(U_L - c_L,\; U_R - c_R,\; \widetilde{U} - \widetilde{c}),
\qquad
S_R = \max(U_L + c_L,\; U_R + c_R,\; \widetilde{U} + \widetilde{c}).
$$

#### 数値フラックス

$$
\mathbf{F}^{\text{HLLE}} =
\begin{cases}
\mathbf{F}_L & (S_L \ge 0), \\[2pt]
\mathbf{F}_R & (S_R \le 0), \\[2pt]
\dfrac{S_R\,\mathbf{F}_L - S_L\,\mathbf{F}_R + S_L S_R\,(\mathbf{Q}_R - \mathbf{Q}_L)}{S_R - S_L} & (\text{中間}),
\end{cases}
$$

ここで $\mathbf{F}_{L/R} = (\rho U,\; \rho u U + P n_x,\; \rho v U + P n_y,\; \rho w U + P n_z,\; \rho H_t U)^\top$、
$\mathbf{Q}_{L/R} = (\rho,\; \rho u,\; \rho v,\; \rho w,\; \rho e)^\top$
(エネルギは $\rho e = P/(\gamma-1) + \tfrac{1}{2}\rho|\mathbf u|^2$、$H_t = (\rho e + P)/\rho$)。

最終的に面積 $S$ を掛けて両側セルへ符号反転加算する。

#### 特性

接触不連続は HLL 平均で塗りつぶされるため、Roe / SLAU よりも数値散逸が大きい。
代償として中間状態に追加の Riemann 解を必要とせず、左右の波速だけで定義できる
ため極めて堅牢で、強衝撃波 (`mach3 wind tunnel` 等) のラフィングに有効。

### KEEP 系 (中心型)

非散逸な KEEP 中心束 (運動エネルギ・エンタルピ保存) に、別スキーム (FVS / SLAU)
の散逸項をブレンドする構成。圧縮性 LES 用途で衝撃波のない領域の数値散逸を抑える。

### AUSM+ / AUSM+UP

質量流束と圧力束を分離する古典的 AUSM 系。AUSM+UP は低マッハ補正項を加える。

## 設定パラメータ

`solverConfig` (YAML) で次を指定する。

- `solver`: 数値フラックス名 (`SLAU` / `HLLE` / `ROE` ほか)。
- `convMethod`: 再構成次数 (`0`: 1 次, `1`: MUSCL 2 次, `2`: MUSCL 3 次, `-1`: ゴースト)。
- `limiter`: リミッタ種別 (実装上は MUSCL 2/3 次は無リミット相当、`-1` で minmod 系を選択)。

## 参考: 残差への寄与

各面で得た $\mathbf{F}^{\text{conv}}_f \cdot \mathbf{S}_f$ を `res_*` (運動量・エネルギ残差)
に符号反転で加算し、左セル `ic0` から減算、右セル `ic1` に加算する。
これが時間積分モジュール ([`docs/time_integration/`](../time_integration/)) に
渡される RHS のうち対流寄与となる。

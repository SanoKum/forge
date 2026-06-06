# 軸対称定式化 — 理論

forge は 3 次元直交座標で書かれた密度ベース有限体積ソルバだが、軸対称な
流れ (周方向 $\theta$ 依存性なし、$u_\theta = 0$) を **2D メッシュ** 上で解く
モードを `isAxisymmetric = 1` フラグで提供する。本ファイルでは離散化方針
(B 流儀: 幾何量を $r$ 重み付け) と、その結果として残るソース項・Roe フラックス
ヤコビアンの固有分解を整理する。実装側の対応は
[implementation.md](implementation.md) を参照。

本ファイルは solver core の理論を扱う。`case/23.axi_nozzle` の検証ノズル形状は
別問題であり、現在の generator は「Witoszynski 収束部 + sharp-throat planar MOC
+ 軸対称出口半径への相似補正」を土台にしている。喉部直後の滑らか化に 5 次多項式
ブレンドを使う場合でも、ここで扱う軸対称支配方程式や B 流儀の離散化自体は変わらない。

ただし、検証ノズルを **axisymmetric MOC** で作る場合は、geometry generator 側で
planar MOC の `K_\pm = \theta \pm \nu` 一定則の代わりに、軸対称の幾何発散に対応する
source-corrected compatibility

$$
\frac{d}{ds}(\theta + \nu) = +S, \qquad
\frac{d}{ds}(\theta - \nu) = -S, \qquad
S = \frac{\sin\mu\,\sin\theta}{r}
$$

を一次近似で積分する。これはノズル壁形状生成の近似モデルであり、solver core の
軸対称有限体積離散化そのものを変更するものではない。

## 座標系と前提

- 対称軸: x 軸固定、半径 $r = y$。任意軸対応は将来拡張とする。
- メッシュ: $(r, z) = (y, x)$ 平面の 2D 三角形 / 四角形セル。
  `case/23.axi_nozzle/run_slau_2d` がリファレンスケース。
- 周方向は $\theta$ 依存性ゼロ、$u_\theta \equiv 0$ を仮定する。
- 本節ではまず非粘性の軸対称保存形を導出するが、実装では粘性
  ($\tau_{\theta\theta}$ 由来) ソース項も追加済みである。陰解法 Jacobian の
  軸対称化は後続 plan に分離する。

## 支配方程式

### 円筒座標の保存形

$u_\theta = 0$, $\partial_\theta = 0$ の軸対称非粘性圧縮流は

$$
\partial_t \rho + \frac{1}{r}\partial_r(r\rho u_r) + \partial_z(\rho u_z) = 0,
$$

$$
\partial_t(\rho u_r) + \frac{1}{r}\partial_r(r\rho u_r^2)
+ \partial_z(\rho u_r u_z) = -\partial_r p,
$$

$$
\partial_t(\rho u_z) + \frac{1}{r}\partial_r(r\rho u_r u_z)
+ \partial_z(\rho u_z^2) = -\partial_z p,
$$

$$
\partial_t(\rho E) + \frac{1}{r}\partial_r(r\rho H u_r)
+ \partial_z(\rho H u_z) = 0.
$$

半径運動量の右辺の $-\partial_r p$ を保存形に直すには、

$$
-\partial_r p = -\frac{1}{r}\partial_r(rp) + \frac{p}{r}
$$

を使う。これにより半径運動量は

$$
\partial_t(\rho u_r) + \frac{1}{r}\partial_r\!\bigl(r(\rho u_r^2 + p)\bigr)
+ \partial_z(\rho u_r u_z) = \frac{p}{r}
$$

の "圧力をフラックスに含めた + $p/r$ ソース" 形になる。
**残った $p/r$ が Cartesian 座標には無い軸対称特有のソース項である。**

### B 流儀: $r$ 倍して planar セルで積分

forge は離散化方式として **B 流儀** = 「両辺に $r$ を掛けてから planar 2D セル
$A \subset (r, z)$ 平面で積分する」を採る。半径運動量に対して両辺 $r$ 倍:

$$
r\,\partial_t(\rho u_r) + \partial_r\!\bigl(r(\rho u_r^2 + p)\bigr)
+ \partial_z(r \rho u_r u_z) = p
$$

(z 偏微分の $r$ は外に出せる)。これを planar セル $A$ で $\iint dr\,dz$ 積分し、
発散項に Gauss を適用すると

$$
\underbrace{\int_A r\,dA}_{=:\,V}\,\partial_t(\rho u_r)
+ \oint_{\partial A}\!\bigl(r(\rho u_r^2 + p)\,n_r + r\rho u_r u_z\,n_z\bigr)\,dl
= \int_A p\,dA.
$$

連続式・軸方向運動量・エネルギー式も同じ手順で

$$
V\,\partial_t Q_i + \sum_{f \in \partial A} \mathbf{F}_f \cdot \mathbf{S}_f = R_i,
\qquad i \in \{\rho,\, \rho u_r,\, \rho u_z,\, \rho E\}
$$

を得る。離散幾何量と各方程式のソースは

$$
V \;=\; \bar r\,A_{\text{planar}}, \qquad
\mathbf{S}_f \;=\; \bar r_f\,dl_f\,\hat{\mathbf n}_f,
$$

| 方程式 | ソース $R_i$ |
| --- | --- |
| 連続 ($\rho$) | $0$ |
| 半径運動量 ($\rho u_r$) | $p_{\text{cell}} \cdot A_{\text{planar}}$ |
| 軸方向運動量 ($\rho u_z$) | $0$ |
| エネルギー ($\rho E$) | $0$ (非粘性) |

ここで $\bar r$, $\bar r_f$ はセル / 面の重心 $r$ 座標、$A_{\text{planar}}$ は
planar セルの 2D 面積 (m²)、$dl_f$ は planar 面 (= $(r,z)$ 平面の線分) の長さ。

### $r$ 重み付けが特異性を吸収する仕組み

「セル中心ソース」式の $p/r$ を planar セル積分するときの単位を追うと

$$
\int_A \frac{p}{r}\,dV \;=\; \int_A \frac{p}{r}\cdot r\,dr\,dz
\;=\; \int_A p\,dr\,dz \;=\; p\cdot A_{\text{planar}}
\quad [\text{Pa}\cdot\text{m}^2]
$$

となり、$dV = r\,dr\,dz$ の中の $r$ が $1/r$ を完全に打ち消す。
**したがって離散ソースに $1/r$ は現れず、軸 ($r \to 0$) 上での
特異性は B 流儀では構造的に発生しない。**

軸 ($r=0$) 上の face は $\bar r_f = 0$ により $\mathbf{S}_f = \mathbf{0}$ となり、
対流フラックスの寄与が自動的に 0 になる。すなわち **「軸 BC を明示的に
指定しなくても、軸面の flux 計算は無害に skip される」** のが B 流儀の
利点である。実装上は意図明示のため `kind: axis` という名称付きで slip
互換のゴーストを置くが、機能的には何もしないのと等価である。

## $2\pi$ の取扱い

完全に revolved した形 ($dV = 2\pi r\,dr\,dz$) で書いても、$V, \mathbf{S}_f, R_i$ の
すべてに $2\pi$ が乗るだけで時間更新 $\Delta Q = -\Delta t/V \cdot
(\text{flux} - R)$ から完全に約分する。forge では既存 3D ケースの `volume` 単位
(m³) と数値スケールが一致するよう **$2\pi$ を付けない**ことに統一する。

物理的な revolved 量 (mdot, 推力など) を出力する場合は **出力ステージで明示的に
$2\pi$ を乗じる**ことで対応する。これは出力側 (`output_*`) の責務とし、ソルバ
コアでは扱わない。

## 軸 ($r \to 0$) の数値挙動

軸上 (= $\bar r \to 0$) のセルでは $V \propto \bar r$ で 0 に近づく。CFL 時間刻み

$$
\Delta t_{\text{CFL}} = \mathrm{CFL} \cdot \frac{V}{\sum_f |\mathbf{F}\cdot\mathbf{S}_f|}
$$

は分子・分母とも $\bar r$ 倍されるため理論上は比が不変。しかし離散誤差で軸上
セルの $\Delta t$ が病的に小さくなることがあるため、実装では軸上セル限定で

$$
\bar r_{\text{eff}} \;=\; \max\!\bigl(\bar r,\; \varepsilon\,h_{\text{local}}\bigr),
\quad \varepsilon \approx 10^{-3}
$$

の下限フロアを `setDT_d` で適用する (詳細は
[implementation.md](implementation.md))。

なお半径運動量の保存変数 $\rho u_r$ 自体は軸上で 0 に収束すべきだが、これは
B 流儀の幾何重みと圧力ソースの組み合わせから自然に成立する性質であり、
明示的なクランプは Phase 1 では入れない。

## Roe フラックスヤコビアン (軸対称 4 波)

軸対称ではフラックス Jacobian の固有値が 4 個になる ($u_n - c$, $u_n$ (entropy),
$u_n$ (shear), $u_n + c$)。$u_\theta$ 方向のせん断波が落ちる分、3D の 5 波から
1 個減る。本節では `convMethod = ROE` での face flux 評価に用いる固有分解を
記す。実装は既存 3D Roe を流用し、軸対称固有の処理は追加しない (理由は
[implementation.md](implementation.md) §"Roe スキームの取扱い")。

ここでは **教科書順** の保存変数

$$
\mathbf{Q} = (\rho,\; \rho u_r,\; \rho u_z,\; \rho E)^\top
$$

を使う。forge の 3D Roe 実装は $(\rho, \rho E, \rho u, \rho v, \rho w)$ 順で
書かれているため、軸対称導出を実装に流用する際は行入れ替えが必要なことに
注意 (詳細は実装側で扱う)。

### 補助量

面法線 $\hat{\mathbf n} = (n_r, n_z)$、Roe 平均量を $\widetilde{(\cdot)}$ で示し、
以下の補助量を定義する。

$$
u_n \;=\; u_r n_r + u_z n_z,\qquad
u_t \;=\; -u_r n_z + u_z n_r,
$$

$$
q^2 \;=\; u_r^2 + u_z^2,\qquad
H \;=\; \frac{\rho E + p}{\rho},\qquad
c \;=\; \sqrt{(\gamma-1)\bigl(H - q^2/2\bigr)},
$$

$$
\beta \;=\; \frac{\gamma - 1}{c^2},\qquad
\alpha \;=\; \frac{\beta\,q^2}{2}.
$$

### 固有分解 $A_n = R\,\Lambda\,L$

$$
\Lambda \;=\; \mathrm{diag}\bigl(u_n - c,\; u_n,\; u_n,\; u_n + c\bigr).
$$

右固有ベクトル行列 $R$ (列が固有ベクトル、列順 = $\Lambda$ の固有値順、
行は $\mathbf{Q} = (\rho, \rho u_r, \rho u_z, \rho E)$ 順):

$$
R \;=\;
\begin{pmatrix}
1 & 1 & 0 & 1 \\
u_r - c\,n_r & u_r & -n_z & u_r + c\,n_r \\
u_z - c\,n_z & u_z & n_r & u_z + c\,n_z \\
H - c\,u_n & q^2/2 & u_t & H + c\,u_n
\end{pmatrix}.
$$

各列の物理意義:

| 列 | 固有値 | 波の種類 |
| --- | --- | --- |
| 1 | $u_n - c$ | 音波 (left-going) |
| 2 | $u_n$ | エントロピー波 |
| 3 | $u_n$ | せん断波 (面接線方向) |
| 4 | $u_n + c$ | 音波 (right-going) |

左固有ベクトル行列 $L = R^{-1}$ (行が左固有ベクトル、行順 = $\Lambda$ 順):

$$
L \;=\;
\begin{pmatrix}
\tfrac{1}{2}\!\left(\alpha + \tfrac{u_n}{c}\right)
& \left(-\tfrac{1}{2c} - \tfrac{\beta u_n}{2}\right) n_r + \tfrac{\beta u_t}{2} n_z
& \left(-\tfrac{1}{2c} - \tfrac{\beta u_n}{2}\right) n_z - \tfrac{\beta u_t}{2} n_r
& \tfrac{\beta}{2}
\\[4pt]
1 - \alpha
& \beta u_n n_r - \beta u_t n_z
& \beta u_n n_z + \beta u_t n_r
& -\beta
\\[4pt]
-u_t
& -n_z
& n_r
& 0
\\[4pt]
\tfrac{1}{2}\!\left(\alpha - \tfrac{u_n}{c}\right)
& \left(\tfrac{1}{2c} - \tfrac{\beta u_n}{2}\right) n_r + \tfrac{\beta u_t}{2} n_z
& \left(\tfrac{1}{2c} - \tfrac{\beta u_n}{2}\right) n_z - \tfrac{\beta u_t}{2} n_r
& \tfrac{\beta}{2}
\end{pmatrix}.
$$

行の物理意義は $R$ の列と対応する (1: 音波-, 2: entropy, 3: shear, 4: 音波+)。

数値フラックスは Roe 散逸付き中央差分

$$
\mathbf F^{\text{Roe}}_f \cdot \mathbf S
\;=\; \tfrac{1}{2}(\mathbf F_L + \mathbf F_R)\cdot \mathbf S
- \tfrac{1}{2}\,S\,\widetilde R\,|\widetilde \Lambda|\,\widetilde L\,
(\mathbf Q_R - \mathbf Q_L),\quad S = |\mathbf S|
$$

として書ける。3D の Roe 表現 ([`docs/convection/theory.md`](../convection/theory.md))
と異なるのは、(i) 固有値が 1 個少ない、(ii) せん断波が 1 個に縮退する、
(iii) $w$ 方向成分が落ちる、の 3 点のみ。

## Phase 1 のスコープと将来拡張

| 項目 | Phase 1 (本理論) | 将来 plan |
| --- | --- | --- |
| 連続 / 運動量 / エネルギー (非粘性) | 対象 | — |
| 半径運動量の圧力ソース $p\cdot A_{\text{planar}}$ | 対象 | — |
| 軸 ($r=0$) BC | 自動退化 + slip 互換のラベル | 任意軸 (BCで指定) |
| 粘性ソース $\tau_{\theta\theta}$ 由来 ($-2\mu u_r/r$ 等) | **非対象** | Phase 2 で追加 |
| 陰解法 Jacobian の軸対称化 | **非対象** | Phase 3 |
| 任意軸対応 (axis_dir 指定) | **非対象** | 別 plan |

## 参考

- Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*,
  Ch.11 (Approximate Riemann Solvers, 軸対称版の固有分解は §11.3 系)。
- Hirsch, *Numerical Computation of Internal and External Flows*,
  Vol.2 §16-3-4 (Axisymmetric formulation)。

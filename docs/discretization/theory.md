# 離散化レイアウト (cell-centered / node-centered) 理論

forge は密度ベース圧縮性 Navier–Stokes を有限体積法 (FVM) で解く。本ドキュメントは、
解変数を **どの幾何要素に乗せるか** = 制御体積 (control volume, CV) の取り方に関する理論を扱う。

- **cell-centered FVM** (現行): primal メッシュのセルを CV とし、解をセル重心に代表させる。
- **node-centered FVM** (本整備の追加対象): primal メッシュの**ノードを CV 中心**とし、各ノードの
  周囲に **中点双対 (median dual)** で CV を構成する。vertex-centered FVM とも呼ぶ。

両者は同じ保存則の積分形を異なる CV 上で離散化したものであり、**数値スキーム
(対流・勾配・粘性・陰解法) の定式化は CV と面 (CV 間の界面) の抽象だけで書ける**。
forge ではこの事実を利用して両対応を実現する (→ [implementation.md](implementation.md))。

## 1. 有限体積法の積分形 (CV 非依存)

保存量 $\mathbf{U}$ について、任意の CV $\Omega_i$ 上で

$$
\frac{d}{dt}\int_{\Omega_i}\mathbf{U}\,dV
+ \oint_{\partial\Omega_i}\big(\mathbf{F}_c - \mathbf{F}_v\big)\cdot\mathbf{n}\,dS
= \int_{\Omega_i}\mathbf{S}\,dV
$$

が成り立つ。$\mathbf{F}_c$ は対流フラックス、$\mathbf{F}_v$ は粘性フラックス、$\mathbf{S}$ はソース項。
離散化では

$$
V_i\frac{d\mathbf{U}_i}{dt}
= -\sum_{f\in\partial\Omega_i}\big(\mathbf{F}_c - \mathbf{F}_v\big)_f\cdot\mathbf{n}_f\,S_f
+ V_i\mathbf{S}_i
$$

となる。ここで必要な幾何量は **CV 体積 $V_i$**、**CV 代表点 $\mathbf{x}_i$**、各 CV を囲む
**面 $f$ の法線 $\mathbf{n}_f S_f$ と面重心 $\mathbf{x}_f$**、そして **面が繋ぐ 2 つの CV** だけである。
この抽象に cell / node の区別は現れない。差異はすべて「CV と面をどう作るか」に押し込められる。

## 2. cell-centered の CV と面

- CV = primal セル。$V_i$ = セル体積、$\mathbf{x}_i$ = セル重心。
- 面 = primal セル間の界面 (内部面) と境界面。面は 2 つのセルを繋ぐ。
- 境界はゴーストセル法で扱う (境界面の外側に鏡像 CV を置く)。

## 3. node-centered の CV と面 (中点双対)

primal メッシュは変えずに、その上に**双対メッシュ**を構成する。

### 3.1 双対 CV (median dual)
ノード $p$ の CV $\Omega_p$ は、$p$ に接続する各 primal セルを
**「セル重心 — 面重心 — エッジ中点 — ノード」** の小片に細分し、$p$ に帰属する小片を集めたもの。

- 2D: 各セル内で、$p$ から出る 2 本のエッジ中点とセル重心を結ぶ四辺形が $p$ の取り分。
- 3D: セル重心・面重心・エッジ中点・ノードから成る四面体片の和。
- 性質: $\sum_p V_p = \sum_{\text{cell}} V_{\text{cell}}$ (双対体積の総和は領域体積に一致)。

### 3.2 双対面 (edge ↔ dual face)
primal の各**エッジ** $(p,q)$ にちょうど 1 枚の双対面 $f_{pq}$ が対応する。
$f_{pq}$ はエッジ $(p,q)$ を囲むセル重心・面重心を結ぶパッチ群であり、
$\Omega_p$ と $\Omega_q$ の界面をなす。法線ベクトルはパッチの面積ベクトル和

$$
\mathbf{n}_{pq}S_{pq} = \sum_{\text{patch}}\mathbf{S}_{\text{patch}}, \qquad
\mathbf{x}_{pq} = \frac{\sum_{\text{patch}} \mathbf{S}_{\text{patch}}\,\mathbf{x}_{\text{patch}}}{\sum_{\text{patch}} \mathbf{S}_{\text{patch}}}
$$

で与える。よって **1 エッジ = 1 面 = 2 CV を繋ぐ** という cell-centered と同型の接続が得られる。

### 3.3 閉性条件 (consistency)
各内部ノード $p$ について、$p$ を囲む全双対面の符号付き面積ベクトル和が消える:

$$
\sum_{q\in N(p)} \mathbf{n}_{pq}S_{pq} = \mathbf{0}.
$$

これは「一定場で勾配がゼロ」「一定流で対流フラックスが保存」を保証する離散整合条件で、
median dual は構成上これを厳密に満たす。実装では前処理で数値的に検証する (→ implementation.md)。

### 3.4 境界 CV と半割双対面
境界ノードは CV 中心が**境界面上に乗る**ため、ゴースト鏡像が定義できない。
代わりに各 primal 境界面を構成ノードへ細分した **半割双対面 (half median-dual face)** を境界ノードに
割り当て、そこで境界フラックスを弱形式で課す。境界面積の分配は保存的でなければならない:

$$
\sum_{p\in f_b} S_{b,p} = S_{f_b}.
$$

## 4. cell-centered と node-centered の比較

| 観点 | cell-centered | node-centered (median dual) |
| --- | --- | --- |
| CV | primal セル | ノード周りの双対セル |
| 自由度数 (tet メッシュ) | $\approx 5.5\,N_{\text{node}}$ | $N_{\text{node}}$ (約 1/5.5) |
| 面ループ単位 | primal 面 | primal エッジ |
| 境界 | ゴーストセル | 半割双対面 + 弱形式 |
| tet 非構造での効率 | 自由度過多 | 自由度少で有利 |
| 軸上特異性 (軸対称) | セル CV が軸から離れる | 境界ノード CV を半割で扱える |

tet 非構造メッシュでは node-centered の自由度が大幅に少なく、計算コスト削減が見込める
(本整備の主動機)。一方 hex メッシュでは自由度比が逆転し得るため、両対応として残す。

## 5. 精度・安定性に関する注意

- **勾配**: Green–Gauss を双対 CV 上で適用する。§3.3 の閉性が一次精度の前提。
- **対流の高次再構成 (MUSCL)**: CV 重心の勾配と CV 重心→面重心ベクトルで外挿する形式は
  双対 CV 上でもそのまま意味を持つ。limiter も CV 値基準で不変。
- **粘性 over-relaxed 補正**: 双対では面法線とエッジ方向が概ね揃い、非直交補正項が小さく安定しやすい。
- **混合要素双対**: セル重心がノード平均近似のとき双対体積も近似を継ぐ。高歪み・非凸要素では
  負体積に注意 (前処理でチェック)。

## 参考

- median-dual / edge-based FVM は非構造ソルバ (SU2, Fluent の node-based など) で広く使われる定式化。
- 実装上の対応とデータ構造は [implementation.md](implementation.md)、
  実装計画は [`.github/plans/discretization-median-dual.md`](../../.github/plans/discretization-median-dual.md) を参照。

# 勾配再構成 — 理論

forge では各セル中心におけるセル平均量 $\phi_c$ から、勾配 $\nabla \phi_c$ を
**Green-Gauss 法 (cell-centered, 線形補間)** で求める。
対流フラックス再構成・粘性フラックス・速度発散 $\nabla \cdot \mathbf{u}$ の評価に用いる。

## 支配関係

セル $c$ (体積 $V_c$、面集合 $\partial c$) について発散定理から

$$
\int_{V_c} \nabla \phi \, dV = \oint_{\partial c} \phi \, \mathbf{n}\, dS .
$$

セル中心勾配を体積平均で近似し、面ごとに代表値 $\phi_f$ を与えると

$$
\nabla \phi_c \;\approx\; \frac{1}{V_c} \sum_{f \in \partial c} \phi_f \, \mathbf{S}_f ,
\qquad
\mathbf{S}_f = \mathbf{n}_f \, A_f ,
$$

ここで $\mathbf{S}_f$ は外向き面ベクトル (面積込み)。

## 面値の補間

forge では面 $f$ を共有する 2 セル $c_1, c_2$ に対し、面ごとの幾何重み $f_x \in [0,1]$
(メッシュ生成時に算出) を用いた線形補間を採用する。

$$
\phi_f = f_x \, \phi_{c_1} + (1 - f_x) \, \phi_{c_2}.
$$

$f_x$ は面中心からセル $c_2$ 中心までの距離を $c_1$–$c_2$ 中心間距離で割った値で、
等間隔メッシュでは $f_x = 1/2$ となる。
歪んだメッシュにおける線形精度を保つための重みであり、本実装では
スキュー補正・最小二乗による高次補正は行っていない。

## 面ループによる加算則

両側のセルへ符号反転して累積する離散版を採用する。
内部面 (両側に実セルが存在) について

$$
\nabla \phi_{c_1} \mathrel{+}= \phi_f \, \mathbf{S}_f, \qquad
\nabla \phi_{c_2} \mathrel{+}= -\phi_f \, \mathbf{S}_f,
$$

最後にセルごとに $V_c$ で除する。境界面については、ghost セルを $c_2$ 側として
扱い、ghost セル値は境界条件モジュール (`docs/boundary/`) が事前に与える。

## 計算する量

各セル中心で次の勾配と発散を同時に計算する。

| 量 | 用途 |
| --- | --- |
| $\nabla U_x, \nabla U_y, \nabla U_z$ | 粘性応力テンソル, 対流再構成 |
| $\nabla T$ | 熱伝導フラックス |
| $\nabla P$ | 対流フラックス再構成 (KEEP, SLAU, AUSM 系) |
| $\nabla \rho$ (GPU 経路のみ) | 密度ベース再構成 |
| $\nabla \cdot \mathbf{u}$ | 衝撃波・人工粘性指標, 圧縮性指標 |

## 既知の制約

- 制限関数 (limiter) は適用していない。スカラー量も含め純粋な Green-Gauss。
- 高次再構成 (二次以上) と最小二乗勾配は未実装。
- 境界面の寄与はゴーストセル値を `c_2` として扱うことで吸収する設計。
  GPU カーネル `calcGradient_b_d` は実装されているが現状コメントアウトされており、
  境界寄与は内部面ループ + ゴーストセル値で表現される。

## 参考

- 実装詳細とソース対応は [implementation.md](implementation.md) を参照。
- アーキテクチャ全体の中での位置付けは
  [`docs/architecture/overview.md`](../architecture/overview.md) を参照。

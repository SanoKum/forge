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

## 6. 境界条件: ゴースト法 vs 弱形式 (node-centered)

### 6.1 なぜゴースト法が node-centered で問題になるか
cell-centered では境界面の外側に鏡像ゴースト CV を置き、対流フラックスも Green–Gauss 勾配も
「ゴーストを隣接 CV とする内部面」として一括処理できる。node-centered では境界ノードが
**境界面上に乗る**ため鏡像ゴーストが退化し (§3.4)、ゴースト位置を人為的にずらす便宜が要る。
この便宜のもとでは:

- **対流フラックス**は境界面でゴースト経由 1 次精度になり、低マッハでは妥当だが、
- **Green–Gauss 勾配**は境界面値を `0.5(φ_node + φ_ghost)` で評価し、境界ノードの半割 CV と
  あいまいなゴースト幾何が混ざる。低マッハでは緩く問題ないが、**高マッハ (強膨張・衝撃) の壁近傍で
  境界ノード勾配が不正確になり 2 次 MUSCL 再構成が過大 → 発散**する (M2 で bump Mach1.65 の
  下壁後縁 x≈1.98 で実証)。

### 6.2 弱形式境界 (weak-form boundary)
ゴーストを介さず、境界 CV の開いた面 (半割双対面) に対し**境界フラックス・境界面値を直接課す**:

- **対流**: 境界半割面 $f_b$ を通る数値フラックスを境界状態 $Q_b$ と内部状態から評価し残差へ直接加える
  $$ \mathbf{R}_p \mathrel{+}= \hat{\mathbf{F}}(\mathbf{Q}_i, \mathbf{Q}_b, \mathbf{n}_b)\, S_b. $$
- **勾配 (Green–Gauss の閉曲面積分の境界寄与)**: 境界 CV の勾配は内部双対面に加え**境界半割面の寄与
  $\phi_b \mathbf{S}_b$ を必ず含めねば閉じない** (発散定理):
  $$ \nabla\phi_p = \frac{1}{V_p}\Big( \sum_{\text{内部面}} \phi_f \mathbf{S}_f + \sum_{\text{境界半割面}} \phi_b \mathbf{S}_b \Big). $$
  ここで $\phi_b$ は境界面での値 (BC 種別で決まる: slip 壁なら接線速度・内部 $\rho,P$、超音速流入なら
  自由流、超音速流出なら内部外挿)。**ゴーストの平均ではなく境界値そのもの**を使う点が肝。

$Q_b$ / $\phi_b$ は各 BC 種別で定まる量で、既存のゴースト状態計算 (slip 鏡像・inlet 規定・outlet 外挿) と
同じ物理を「面の値」として与えるもの。実装では既存ロジックを流用しつつ、勾配・フラックスへの
寄与だけを弱形式に置き換える (→ [implementation.md](implementation.md) §7)。

### 6.3 壁の強形式 (Dirichlet) と壁ゴースト撤廃

粘性壁 (no-slip) は弱形式フラックスより**強形式 (Dirichlet) の方が node-centered では自然かつ堅牢**である。
理由は §3.4 の通り**境界ノードが壁面そのものの上に乗る**ため、その点の速度は物理的に厳密にゼロでなければ
ならない:

$$ \mathbf{u}_p = \mathbf{0} \qquad (p \in \text{wall}). $$

cell-centered ではセル中心が壁から半セル内側にあり、壁面での $\mathbf{u}=0$ は鏡像ゴースト
($\mathbf{u}_{\text{ghost}}=-\mathbf{u}_{\text{interior}}$) で**面平均**を 0 にして表す。これはセル中心速度が非ゼロ
(近壁速度) であることと整合する。しかし node-centered で同じ鏡像を使うと、面平均は 0 になるが
**壁ノード自身の速度は 0 に固定されない** (残差で更新される)。これは壁ノードが壁上にあるという幾何と矛盾する。

**壁ゴーストは不要**:
- **断熱壁を通る対流フラックスは恒等的に 0**。$\mathbf{u}\cdot\mathbf{n}=0$ より質量・エネルギー流束 0、運動量は
  Dirichlet で固定 (壁反力が圧力を吸収)、断熱で熱流束 0。よって壁半割面の寄与は丸ごと不要で、CV を
  「壁側で閉じない」まま扱ってよい (落とした圧力起因の偽運動量は運動量 Dirichlet が吸収する)。
- **壁せん断の物理は内部双対面が担う**。近壁の内部ノード $q$ が感じる粘性せん断 $\mu\,\partial u_t/\partial n$ は、
  $q$ と壁ノード $p$ ($\mathbf{u}_p=0$) を繋ぐ**通常の内部双対面** $f_{pq}$ の粘性フラックスとして自然に現れる。
  したがって壁ゴーストを撤廃しても近壁せん断は失われない (むしろ鏡像で実効距離が 2 倍になる歪みが消える)。

> **実装上の補足 (検証で更新)**: 残差射影 Dirichlet は壁の**圧力寄与も落とすため不正**で発散した。また壁ゴーストを
> 完全に消すと壁ノードの Green-Gauss 勾配が閉じず悪化した。よって実際の採用形は: **壁ゴーストは勾配/粘性の
> ミラー用に残し、対流のみ壁境界値 (bvar 速度=0) で pressure-only に弱形式化**する (no-slip の本質=対流の運動量は
> 圧力のみ・せん断は内部双対面が担う、を満たす)。さらに**近壁極小 CV の粘性は陰解法では十分 implicit 化されず
> explicit (cfl≤0.1) が安定**。詳細・検証は [implementation.md](implementation.md) §7.2。

#### 6.3.1 陰解法での Dirichlet 整合 (Jacobian 行 decouple)

壁 $\mathbf{u}_p=\mathbf{0}$ を**残差射影だけ** ($\mathbf{res}_{\rho\mathbf{u}}=0$) で課すと、explicit では $\Delta(\rho\mathbf{u})=0$ になり厳密に効くが、**implicit (block-DPLUR) では不十分**である。block-DPLUR の対角 5×5 ブロックは壁ノードの運動量を連続・エネルギーと連成したままなので、線形 solve は残差 0 でも $\Delta(\rho\mathbf{u})\neq0$ を返し、**壁速度が drift する** (実測: 壁 $|\rho u_x|$ が時間発展で増大、状態射影 `enforceWallNoSlip` が毎ステージ余分な KE を剥ぐ持続的エネルギーシンク化)。

SU2 (同じ vertex-centered median-dual) はこれを **`Jacobian.DeleteValsRowi`** で解決する: 壁ノードの**運動量行を単位行に置換** ($\Delta(\rho u)=\Delta(\rho v)=\Delta(\rho w)=0$ を solve から直接得る)。連続・エネルギー行は保存式のまま残し、$u=0$ より $\rho e = \rho e_\text{int}$、圧力は EOS が復元する (CPG: $P=(\gamma-1)\rho e$、TP: $\rho e_\text{int}\to T\to P=\rho R_\text{mix}T$)。decouple は運動量行のみを触り thermo 行に依らないので **CPG/TP 共通** に動く。forge では既存の軸 `axis_flag` 行 decouple 機構を**壁の運動量3行へ拡張**して実装する。計画: [`discretization-node-wall-implicit-dirichlet.md`](../../design/accepted/discretization-node-wall-implicit-dirichlet.md)。

> 軸 (§7.1) の `roUy` 単行 decouple は r=0 の極小 CV が多方程式的に ill-posed なため単独では発散したが、**壁は r=0 特異点を持たない**ので、運動量3行 decouple が drift 機構を直接除去できる。

### 6.4 コーナー所有優先 (boundary∩boundary)

2 つの境界に同時に乗るノード (例: 壁∩出口の出口リップ) は、ゴースト法では**境界ごとに別々のゴースト**が
立ち、同一 CV の残差に矛盾する条件 (壁の $\mathbf{u}=0$ と出口の超音速外挿) が同時に課されて発散する。

これを避けるため、各境界ノードを**優先度に従い 1 つの境界条件にのみ所有させる**:

$$ \text{wall} > \text{inlet} > \text{outlet} > \text{slip/axis}. $$

物理的に、壁とその他境界の交線上の点は**壁点**である (リップは固体壁)。よって壁が最優先で所有し、当該ノードは
壁 Dirichlet として扱い、他境界 (出口等) の条件・半割面からは除外する。これによりコーナーの矛盾が解消する。

> 段階導入: まず壁を強形式化＋壁ゴースト撤廃 (Phase 1)、続いて流入出/slip/axis のゴーストも半割面弱形式へ
> 置き換え node モードを完全 ghost-free にする (Phase 2)。

## 参考

- median-dual / edge-based FVM は非構造ソルバ (SU2, Fluent の node-based など) で広く使われる定式化。
- 実装上の対応とデータ構造は [implementation.md](implementation.md)、
  実装計画は [`.github/plans/discretization-median-dual.md`](../../design/active/discretization-median-dual.md) を参照。

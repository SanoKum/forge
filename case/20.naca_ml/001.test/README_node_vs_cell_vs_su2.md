# NACA2412 node vs cell vs SU2 比較 (Euler / 層流 / SST)

forge の **node-centered (median-dual) と cell-centered** 離散化を、同一メッシュ・同一 BC で
**Euler・層流・SST** の 3 物理について比較し、中立基準として **SU2 v8.5.0** を併置した検証。

## 結論サマリ (2026-06-26 更新)

**当初 node が「物理に鈍感で数値散逸が大きい」ように見えたのは、node 安定化で入れた
`bndFirstOrder=1`(境界隣接 CV の再構成勾配をゼロ化＝壁を1次精度化) が原因だった。**
`bndFirstOrder=0` にすると:

- **Euler**: node が cell・SU2 と一致 (CL 0.889/0.888/0.885、CD 0.0047/0.0044/0.0048、吸込みピーク Cp=−2.0 が重なる)。
- **物理感度が回復**: node の CL が Euler 0.889 / 層流 0.882 / SST 0.429 と物理ごとに変化
  (bnd1 では全物理 0.82 固定だった)。境界層も物理に応じて形成 (mid-chord 第1オフ壁 |U|:
  euler 123 / 層流 98 / SST ~1)。

→ node-centered (median-dual) 自体は健全 (SU2 も同方式)。`bndFirstOrder=1` の壁1次化が
精度・物理感度を殺していた。亜音速では `bndFirstOrder=1` は不要・有害。
本筋の修正は SU2 流の「壁勾配を slip 補正 (Blazek 8.40)＋弱形式 BC 流束」化
([discretization-node-boundary-ghostless.md])。

## 条件

- 形状: NACA2412、単位翼弦 c=1、2D 単層押し出しメッシュ (34691 cells、forge `naca.h5` = SU2 `naca.su2` で同一)。
- 流れ: 入口 `inlet_Pressure_dir` (Pt=107853.4 Pa, Tt=293.337 K, 方向 (1,0,0))、出口 `outlet_statPress` Ps=101325 Pa。
  → **M∞ ≈ 0.30、α=0**(メッシュ実効 AoA により有限揚力)、Re_c ≈ 6.9×10⁶、q∞=6383 Pa。
- top/bot=slip、side1/side2=2D 対称面、翼面=Euler は slip / 粘性は no-slip。
- forge: explicit RK3 (timeIntegration=3)、MUSCL+Venkatakrishnan、visc=1.8e-5 定数 (層流/SST)。
- SU2: ROE+MUSCL+Venkat、EULER_IMPLICIT、MU_CONSTANT=1.8e-5、SST は TI=0.8%/μt/μ=67 (forge の k=1/ω=1000 相当)。
- CL/CD は翼面圧力積分 (forge: 壁面 `res_wall_5` の CONNE 面、SU2: `surface.csv` の保存量→P)。q∞ 正規化は両者一致。

## 計算 run 一覧

| run | discretization | 物理 | 主結果 (CL / CD) | 状態 |
| --- | --- | --- | --- | --- |
| `run_cmp_cell_euler` | cell | Euler | CL 0.888 / CD 0.0044 | still-converging |
| `run_cmp_node_euler` | node | Euler | CL 0.820 / CD 0.0210 | plateau |
| `run_cmp_cell_lam`   | cell | 層流  | CL 0.869 / CD 0.0058 | NOT CONV (非定常) |
| `run_cmp_node_lam`   | node | 層流  | CL 0.820 / CD 0.0210 | NOT CONV (非定常) |
| `run_cmp_node_sst`   | node | SST (wt=0) | CL 0.819 / CD 0.0210 | plateau |
| `run_cmp_cell_sst_wt1` | cell | SST (wt=1) | CL 0.794 / CD 0.0075 | DIVERGED (ω 壁 NaN, flow は完走) |
| `run_cmp_node_euler_bnd2` | node (bnd0) | Euler | CL 0.889 / CD 0.0046 | **cell/SU2 と一致** |
| `run_cmp_node_lam_bnd0` | node (bnd0) | 層流 | CL 0.882 / CD 0.0047 | NOT CONV (非定常) |
| `run_cmp_node_sst_bnd0` | node (bnd0) | SST | CL 0.429 / CD 0.0368 | NOT CONV (厚 BL) |
| `run_su2_euler` | SU2 | Euler | CL 0.885 / CD 0.0048 | conv rms-7.8 |
| `run_su2_lam`   | SU2 | 層流  | CL 0.472 / CD 0.0364 | conv rms-6.8 |
| `run_su2_sst`   | SU2 | SST   | CL 0.602 / CD 0.0290 | conv rms-7.1 |

図: `cp_euler.png` / `cp_lam.png` / `cp_sst.png` (node/cell/SU2 の Cp 分布)。

## 主要所見

### 1. Euler — node の差は `bndFirstOrder` が原因 (最も明快)

亜音速 Euler は **CD=0 が厳密解** (d'Alembert) なので、CD は数値散逸の指標になる。

| | CL | CD (≈0 が理想) | 前縁吸込みピーク Cp |
| --- | --- | --- | --- |
| forge cell | 0.888 | 0.0044 | −2.0 |
| SU2        | 0.885 | 0.0048 | −2.0 |
| **forge node `bndFirstOrder=0`** | **0.889** | **0.0046** | **−2.0** |
| forge node `bndFirstOrder=1` (旧) | 0.820 | 0.0210 | −1.5 |

- **forge cell ≈ SU2 ≈ node(bnd0)** (Cp 分布が重なる、CL/CD 一致)。図 `cp_euler_bnd.png`。
- **`bndFirstOrder=1` のときだけ** node は前縁吸込みピークを鈍化 (Cp −1.5) させ spurious CD 4–5 倍。
  原因は [calcGradient_d.cu:788](../../solver_density_cuda/cuda_forge/calcGradient_d.cu#L788) が
  **境界隣接 CV の再構成勾配を全ゼロ化**＝翼面接線方向の圧力勾配を消して1次精度化するため。
  → node-centered 自体は健全。`bndFirstOrder=1` (亜音速では不要) が精度を殺していた。

### 2. 層流 (Re 6.9×10⁶) — 本質的に非定常、forge 未収束

- **SU2** は層流分離を捕捉し定常化 (吸込みピーク −1.0 へ崩壊、**CL 0.885→0.472**)。
- **forge cell/node** は explicit で**未収束 (非定常プラトー)**、分離を捉えず Cp は Euler 的のまま。
  → 高 Re 層流は層流 BL が分離して非定常になる領域で、explicit forge では定常比較が取りにくい
    (ユーザー予想どおり「粘性計算は暴れる」)。node/cell の優劣判定には不適。

### 3. node の「物理に鈍感」も `bndFirstOrder=1` が原因 (解決)

`bndFirstOrder=1` では forge node の CL は **Euler 0.820・層流 0.820・SST 0.819 とほぼ不変**だった。
壁の1次化散逸が物理粘性 BL を埋もれさせ、node_lam の第1オフ壁速度が node_euler と同一 (123 m/s, BL 無し) になっていた。

**`bndFirstOrder=0` で物理感度が回復**:
| physics | node bnd1 (旧) | node bnd0 (新) | cell | SU2 |
| --- | --- | --- | --- | --- |
| Euler | 0.820 | 0.889 | 0.888 | 0.885 |
| 層流 | 0.820 | 0.882 | 0.869 | 0.472 |
| SST | 0.819 | 0.429 | 0.794 | 0.602 |

mid-chord 第1オフ壁 |U|: euler 123 / 層流 98 / SST ~1 と**物理ごとに BL が形成**。
(SST bnd0 は CL 0.429＝BL 過厚で振れ; 層流/SST とも Re6.9e6 で非定常・未収束のため値は plateau。)

### 4. node SST 壁処理のロバスト性が cell と逆 (要対処)

- **node SST**: `wallTreatmentSST=1` (壁関数) は **ω が壁で発散** (omega 先頭、cold IC step~1)。
  `wallTreatmentSST=0` (wall-resolved) のみ安定 (本比較は wt=0)。既知の node 壁関数バグ系
  ([turbulence-node-sst-wallfunction], [diffusion-node-wall-viscous-distance])。
- **cell SST**: `wallTreatmentSST=1` で安定、`wt=0` は逆に ω 即発散 (y+~50 で wall-resolved 不整合)。
- → **同一 wt で両者を回せず**、node は wt=0・cell は wt=1 を強制 (壁処理ロバスト性が逆転)。
  公平な SST 直接比較が難しく、これ自体が node の制約。
- 注意 (実装上の落とし穴): SST を **k/ω を持たない層流場から seed すると ω 即 NaN**。
  cold IC (naca.h5 の freestream 初期化) か、k/ω を持つ場から restart すること。

## 結論

- **cell-centered は SU2 と整合** (Euler で Cp・CL・CD 一致)。
- **node-centered (median-dual) 自体も健全**: `bndFirstOrder=0` にすれば Euler で cell・SU2 と一致し
  (CL 0.889/0.888/0.885)、物理感度も回復する。当初の「node は散逸大・物理に鈍感」は
  **私が安定化で入れた `bndFirstOrder=1` (壁1次化) のアーティファクト**だった。
- **対処**: 亜音速で `bndFirstOrder=1` は不要・有害。既定や運用を見直す。本筋は SU2 流の
  「壁勾配を slip 補正 (Blazek 8.40)＋弱形式 BC 流束」化 ([discretization-node-boundary-ghostless.md])。
- 層流は Re 6.9e6 で非定常、SST は wt 非対称＋過厚 BL で、定常 node-vs-cell 比較は **Euler が最も信頼できる**
  (そこで node=cell=SU2 が確認できた)。

## 追補: node 境界半割面拡散 skip の回帰 (2026-06-27)

plan [diffusion-node-boundary-real-distance.md](../../../.github/plans/diffusion-node-boundary-real-distance.md)
で node 境界半割面の k/ω 拡散を ∇·S 弱形式→**skip** 化 (commit af5b98d)。本 NACA-ML node SST でも
平板 case/26 と同じ差分回帰を実施:

- **方法**: 収束済み node SST 場 `run_cmp_node_sst_wt1_fsinit/res_120000` を同一 restart として、**旧 ∇·S
  バイナリ (`run_node_sst_gradS_verify`)** と **新 skip バイナリ (`run_node_sst_skip_verify`)** で +5000 step
  継続 (wt=1, MUSCL, RK3)。
- **結果**: **CL・CD が完全一致** (両者 CD=0.0283 / CL=0.312、同一簡易圧力積分器で 5 桁一致)。全流れ場も
  skip vs ∇·S で relL2 ≤ 2e-4 (ro 9.6e-6, P 1.1e-5, Ux 9.3e-5)。唯一 k が relL2 8.9e-4 と僅差で、差は
  翼面近傍 (wall_dist<1mm に 87%, max|Δk|=3.97) に局在し far/slip 域は 0% → 平板同様 **空力係数は無影響**。
- NaN なし (5000 step、全 rms 列非有限なし)。両 run の `check_convergence.py` VERDICT は **NOT CONVERGED
  (plateau)** = 収束済み場からの継続でプラトー床 (base 自身が block-DPLUR プラトー、新規収束ではない)。
  → ここでの主張は「収束」ではなく「skip 化が収束場の CL/CD を動かさない (差分回帰)」。
- 図: `run_node_sst_skip_verify/residual_history.png`。**結論: skip 化は NACA-ML node SST でも空力係数不変で安全**。

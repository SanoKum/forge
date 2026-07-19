# node 境界半割面 拡散の実距離 over-relax 化 (∇·S 弱形式の置換)

## メタ

- **area**: `diffusion`
- **status**: `done`
- **related_docs**:
  - `methods/diffusion.md`
  - `methods/diffusion.md`
- **related_plans**:
  - `.github/plans/discretization-node-boundary-ghostless.md` (親: node 境界 ghostless)
  - `.github/plans/architecture-node-centroid-value-position.md` (centCoords=node 統一)
  - `.github/plans/diffusion-node-wall-viscous-distance.md`
- **created**: `2026-06-26`
- **owner**: `CFD Dev`

## 1. 目的

node (median-dual) 境界半割面の拡散流束を、**ghost mirror (dcc≈0) を使う退化形 / 暫定の ∇φ·S 弱形式**から、
**境界ノードと法線内部ノードの実距離による over-relax 形**へ統一する。退化割り算を原理的に排し、
境界種別 (Dirichlet/Neumann) を正しく反映する。

## 2. 背景

node 境界では境界ノード W が境界面上に乗り、ghost mirror が dcc≈0 に退化 → Fick `ρD(Y1-Y0)/dcc·δ` が
0/0。対策として現状 2 形が混在:
- **退化形** (旧): ghost dcc をそのまま使用 → garbage (species 拡散 `species_diffusion_d` が該当だったが
  本セッションで暫定 ∇·S 化)。
- **∇φ·S 弱形式** (暫定): cell 勾配を流用し `ρD(∇φ·S)`。退化は回避するが (a) cell 勾配の境界閉包に依存、
  (b) Neumann でも厳密 0 にならない、という弱点。k/ω 拡散 (`scalar_diffusion_first_order_d:84`) と
  species 拡散 (本セッション) が該当。

ユーザー方針: **「壁面ノードと内部ノードの実距離を測って over-relax」**で内部面と同形に統一したい。

## 3. ⚠ 設計の核心論点 (実装前に確定が必要)

境界半割面の拡散流束が **BC 種別でどうあるべきか**を先に確定する。候補:

| 扱い | Neumann (slip/断熱/zero-grad) | Dirichlet (壁固定値/inlet) |
| --- | --- | --- |
| (a) ∇φ·S 弱形式 (現状) | cell 勾配次第で微小非0 (誤差) | cell 勾配の境界投影 |
| (b) 実距離 `ρD(φ_W-φ_I)/d_n` (ユーザー案) | **(φ_W-φ_I)/d ≠ 0 → 過剰付与 (誤)** | 壁法線フラックス (正) |
| (c) 境界半割面を 0 にスキップ | 0 (Neumann 厳密・正) | W はピンで上書き → 内部面が担うので 0 でよい? |

**未確定の問い**: 境界半割面フラックスは「内部 W↔I 双対面が既に実距離で担う拡散」と**別物として加えるべきか**、
それとも **Dirichlet は W ピン＋内部面で足り、Neumann は 0** で、**半割面は本来スキップが正しい**のか。
- 1D 検証: 壁 Dirichlet では W=固定値、内部面 W↔1 が `ρD(φ_W-φ_1)/d` を node1 に運ぶ (実距離・正)。
  → 半割面を 0 にしても node1 は正しい。Neumann も内部面＋半割面0で zero-grad に収束。**(c) で足りる可能性**。
- ただし k/ω 拡散で ∇·S 化が omega 安定に効いた実績があり、(c) との差異要検証 (omega ピン箇所では
  半割面フラックスは W に入り上書きされるはずだが、実測で挙動確認)。

→ **(b) を一律適用すると Neumann で過剰付与の懸念**。ユーザー案は壁 (Dirichlet 的) には妥当だが、
slip/出口 (Neumann) には不適。

### 暫定結論 (2026-06-27 議論): **(c) skip が正しく最もシンプル**

検証で詰めた論理:
- **Dirichlet** (壁 u=0 / ω ピン / 入口 / 等温壁 T=T_wall): 境界ノードはピン or 値固定 → 半割面に何を足しても
  上書きされ無意味。**内部ノード I は内部双対面 W↔I が実距離で拡散を運ぶ** (ROE_d/viscousFlux_d の主ループで計算、
  下記検証で確認)。
- **Neumann** (slip / 断熱 / zero-grad 出口): 物理的に $\partial_n\phi=0$ ＝半割面フラックス 0。**skip がそのまま正解**。
- 例外は**指定非ゼロ Neumann 流束 (壁熱流束 q≠0 を陽に課す) のみ**。forge の壁は断熱/Dirichlet/等温(=Dirichlet)
  なので不要。等温壁の熱流束は viscousFlux 側 (scalar 拡散ではない)。
- → 実距離 over-relax (b) は「半割面に新規フラックスを足す」案だが、実距離拡散は**内部面 W↔I が既に担う**ので、
  半割面は「足さない (skip)」で足りる。

### 補足検証 (本セッションで確認済)
- **接線 W↔W' 面 と 境界→内部 W↔I 面 は両方とも ROE_d/SLAU_d/viscousFlux_d の主ループ (`convPlaneBound=nNormalPlanes`,
  `normal_halo_planes_d[0:convPlaneBound]`) で計算されている** (bc.iPlanes でない内部 plane)。
  `convectiveFlux_boundary_d` は境界半割面 (node↔ghost, 外向き BC 流束) 専用。
- node 双対 CV の**閉性残差 1.8e-7** (全境界ノード) → 全面 present・完全。slip 接線輸送の欠落は無い (杞憂だった)。

## 4. 次の一手 (実装＋検証)

**(c) skip を実装**: `scalar_diffusion_first_order_d` と `species_diffusion_d` の node 境界半割面分岐
(現状 ∇φ·S 弱形式) を **`return` (skip, 残差に何も足さない)** に置換。粘性運動量/エネルギーは
viscousFlux_d が既に境界を別処理 (主ループは [0:nNormalPlanes] で境界除外) なので対象外。

検証: 平板 `case/26.flat_plate_sst` (壁解像 node SST) で **Cf が現状 (∇·S) と不変**であること (=skip が正しい・
弱形式の壁寄与は元々ピンで上書きされていた) を確認。必要なら (a)/(c) を並べて Cf/δ99/log則を SU2 と比較。
OK なら species へ展開。cell は isNode ガードでビット不変。

## 5. スコープ

- **やる**: 境界半割面拡散の正しい扱いを検証で確定 → 共通ヘルパ (境界ノード→法線内部ノード＋距離＋
  over-relax 係数、representative point は壁関数 `bestI` を一般化) → k/ω で実証 → species/viscous/energy へ展開。
- **やらない**: 対流の ghostless 化 (本セッションで別途実施済)。flow(ρ,ρu,ρe) の bvar 弱形式 (既に正しい)。
- **cell 不変**: 全変更 `isNode` ガードでビット不変。

## 6. 段階

1. (検証) flat_plate で (a)/(b)/(c) の Cf/BL 比較 → 正しい形を確定。
2. 共通ヘルパ実装 (境界半割面の representative 内部ノード＋距離＋over-relax)。
3. k/ω 拡散を確定形に置換 (∇·S 撤去)、平板で SU2 一致確認。
4. species・粘性運動量・エネルギー拡散へ展開。各段で cell ビット不変＋node 検証。

## 変更ログ

- 2026-06-26: draft 作成。実装前に §3 の核心論点 (Neumann で実距離は過剰付与/半割面スキップが正しいか)
  を §4 の平板検証で確定する方針。
- 2026-06-27: **(c) skip を実装・検証完了 (status: done)**。
  - **実装**: `scalar_diffusion_first_order_d` (k/ω) と `species_diffusion_d` の node 境界半割面分岐
    (`isNode!=0 && (ic0>=nCells || ic1>=nCells)`) を ∇φ·S / ∇Y·S 弱形式から **`return` (skip)** に置換。
    勾配引数 (`dphidx`/`dYdx` 等) は未使用化するがシグネチャ・呼び出し側は維持 (cell churn 最小化)。
  - **検証** (平板 `case/26.flat_plate_sst`, node SST MUSCL, 収束場 run_node_sst_muscl_cont/res_90000 から +5000 step):
    - 旧 ∇·S バイナリ (`run_node_gradS_verify`) と 新 skip バイナリ (`run_node_skip_verify`) を**同一 restart から
      同ステップ継続**して apples-to-apples 比較。
    - **Cf / u_τ / δ99 は 3 station (x≈0.30/0.60/0.89) で完全一致** (skip vs ∇·S 差 0.00%、両者とも元収束場
      ref90k と 0.01% 以内)。**Cf 不変を確認** → §4 の判定 PASS。
    - 唯一の場差は k (relL2 6.5%) で、**全て前縁上流 x<0 のスリップ壁・y≈0 域に局在** (near-wall y<1mm に 100%)。
      そこでは ∇·S 弱形式が Neumann スリップで漏れて k を ~10 に過剰生成していたのに対し、**skip は ∂k/∂n=0 を
      厳密に満たし k≈3 (=内部値) を保持** → §3(c) の予測どおり skip がむしろ物理的に正しい。平板 BL (x∈[0,1]) は無影響。
    - NaN/発散なし (5000 step、rms 全列で非有限なし)。
  - **cell ビット不変**: 変更は全て `isNode!=0` 分岐内。cell SST (run_0001) を新旧バイナリ 2 step で比較すると
    `ro 1.19e-7 / roUx 7.6e-6 / roOmega 4.0` 差だが、**同一新バイナリ 2 回の atomicAdd 非決定論の床と完全一致** →
    変更由来でなく atomicAdd ノイズ。cell は床内でビット不変を実証。
  - **残課題**: species (`species_diffusion_d`) も同方針で skip 化したが node+species の検証ケースは未整備
    (k/ω パターン準拠だが未実測)。将来エネルギー拡散へ skip 拡張時は等温壁 (指定非ゼロ熱流束 Neumann) に注意
    (forge の現状は断熱/Dirichlet なので不要、熱流束は viscousFlux 側)。

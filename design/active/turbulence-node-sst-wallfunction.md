# node-centered SST 壁関数の代表点修正 (ic=壁ノード u=0 トラップ)

## メタ

- **area**: `turbulence`
- **status**: `in_progress` (SU2 調査完了→方針 A 確定、実装中 2026-06-24)
- **related_docs**:
  - `docs/turbulence/implementation.md`
  - `.github/plans/diffusion-node-wall-viscous-distance.md` (§11 twall, 同一トラップ)
- **created**: `2026-06-24`

## 1. 確定済みの問題 (コード + 実測)

`ransWallFunction_d.cu` (`computeWallFrictionSST_d`, `wallTreatmentSST=1`) は `ic=bplane_cell[ib]` から
`ro[ic]`/`vis_lam[ic]`/`wall_dist[ic]`/`Ux,Uy,Uz[ic]` を読み `wf_pk[ic]` に書く。**node では ic=壁ノード W**
で no-slip Dirichlet により `Ux=Uy=Uz=0`・`wall_dist≈0`。よって ([L62,75,84-89](../../solver_density_cuda/cuda_forge/ransWallFunction_d.cu#L62)):

```cpp
const geom_int ic = bplane_cell[ib];   // node: 壁ノード W
const flow_float uc = Ux[ic];          // = 0 (Dirichlet)
... Ut = |U_c-(U_c·n)n| = 0 ...
if (Ut <= kSmall) { utau_b=0; ypls_b=0; wf_pk[ic]=0; return; }
```

→ **node の全壁ノードで `utau=0, ypls=0, wf_pk=0` に退化** (実測: `run_node_wallstress_on/res_wall_4_30.h5`
で utau/ypls の nonzero frac = 0.00)。これは [[node-sst-wallfunction-utau-zero]]・twall≡0 と**同一の
ic=壁ノード u=0 トラップ**。cell は ic=内部セル (u≠0) で正常。

**症状 (case/36)**: 壁生産退化 + ω 壁 BC 縮退で node 近壁 μt が異常 (μt/μ~4-5千がコアへ侵入) → BL 厚化
→ 発散ダクトの超音速コア加速失敗 (M 1.69 vs cell/SU2 1.82) → 擬似衝撃波が ~80mm 上流 (node 46mm vs
cell 142 / SU2 132)。背圧一致 (Ps1.90) でも node shock 46mm に漸近 (背圧交絡でない)。平板 (case26,
`wallTreatmentSST=0`) は壁関数不使用ゆえ node/cell 一致 (Cf±1%)。詳細: diffusion plan §11.5 / 分析応答。

## 2. 修正設計の論点 (ユーザー memo 2026-06-24, SU2 調査前の暫定)

`Ut` だけ内部ノードから取るのは**不十分**。`wall_dist[W]≈0` なので Reichardt 逆解きの `y` が不整合。
少なくとも次を**一貫して内部代表点/壁法線距離から**取る:

- `Ut` (接線速度), `rho`, `nu=vis_lam/rho`, `y` (壁法線距離)

### 2.1 複数内部エッジの代表点選択 (twall とは別方針)

twall は全内部エッジの粘性力を**集約** (保存的, bb90036)。だが SST 壁関数は「壁から距離 y の代表点速度」
モデルなので**単純な総和/平均は危険**。暫定推奨:

- `bplane` ごとに壁法線 `n` を使う。
- `W=bplane_cell[ib]` に接続する内部ノード候補 `I` を走査。
- `dn=(x_I-x_W)·n>0` かつ**最小**の候補を代表点に (= 第一オフ壁点; `wall_y_eff` と同趣旨)。
- 接線ズレ大の候補は除外/低優先。
- `Ut,rho,nu,y` をその代表 `I` と `dn` で一貫計算。`utau_b[ib],ypls_b[ib]` は bplane 値で書く。

### 2.2 `wf_pk` overwrite 問題

`wf_pk[ic]=...` は直接代入。node では同一 W に複数 bplane/wall face が対応し得 **非決定的 overwrite**。
候補: (a) W に面積重み平均で集約, (b) 最大値, (c) 適用先を代表内部点 I に変更, (d) bplane 変数化。
**最小変更の低リスク案**: 参照値は内部点から、適用先は現行どおり W、複数寄与は**面積重み平均**。
注意: (c) `wf_pk[I]` 案は `ransBoundary_d.cu` が node で `omega[W]` を直接ピン留めしている件と整合要。

## 3. SU2 調査結果 (完了 2026-06-24, `.external/su2-src` 12eb826)

SU2 の vertex-centered wall function は**壁頂点自身の速度を使わず、各壁 vertex の `Normal_Neighbor`
(壁法線と cos 最大の隣接内部点) を代表点**にする。
- `FindNormal_Neighbor` ([CPhysicalGeometry.cpp:7660](../../.external/su2-src/Common/src/geometry/CPhysicalGeometry.cpp#L7660)): 隣接点のうち壁法線との cos 最大を選ぶ。
- `CNSSolver::SetTau_Wall_WF` ([CNSSolver.cpp:798](../../.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp#L798)): `Point_Normal` の速度 `VelTangMod` と距離 `WallDistMod` で U_Tau/YPlus/Tau_Wall。
- `CTurbSSTSolver.cpp:523`: SST 壁関数は流れ側 YPlus/UTau を使い Normal_Neighbor の k/omega を更新。
- **複数 edge の面積平均ではなく代表内部点 1 点**。wall_dist=0 の壁頂点値は代表距離に使わない。

### 3.1 確定方針 (A, SU2 整合)

`compute_wall_friction_sst_d` を cell/node 分岐:
- **cell**: 現行どおり `irep=ic`(内部隣接セル), `y=wall_dist[ic]` → **ビット不変**。
- **node**: `W=bplane_cell[ib]`(壁ノード)。`W` の入射内部双対面 (`cell_planes` CSR, `ip<nNormalPlanes`,
  `wall_flag[I]==0`) を走査し、**壁内向き法線 `-n` との cos = `-(x_I-x_W)·n/|x_I-x_W|` が最大の I を代表点**に選ぶ
  (SU2 `Normal_Neighbor`)。`irep=I`、`y=dn=-(x_I-x_W)·n`(壁法線距離>0)。`Ut/rho/nu` を全て `I` から一貫取得。
  候補なしは退避 (utau=0)。
- 壁単位法線 `n` は `bplane_plane[ib]` の `sx/sy/sz/ss`。`utau_b[ib]`/`ypls_b[ib]` は bplane 値。
- `wf_pk` は現行どおり**壁ノード W に適用** (omega ピン留め先 W と一致, ransBoundary 再設計を避ける)。
  1 wall bcond 内では W↔半割面が 1:1 なので**bcond 内 race なし** (直接代入で可)。複数 wall bcond が共有する
  コーナーのみ overwrite し得るが、検証ケース (case36/26 とも単一 wall bcond) では発生しない。**複数 wall bcond
  コーナーの面積重み集約は将来の堅牢化** (本 fix の物理・検証には不要)。
- y-整合: utau の `y=dn`(代表点壁法線距離) は ransBoundary の omega 用 `wall_y_eff`(最近接隣接 wall_dist) と
  近似一致 (構造近壁メッシュでは代表点=最近接)。厳密一致は将来 refine。

## 3.2 第2の欠落＝壁関数 τ_w が運動量に課されていない (SU2 AddTauWall, 2026-06-24)

§3.1 の代表点修正 (utau/ypls 是正) は必要だったが**不十分**だった (case36 shock 46→65mm 止まり)。SU2 を
さらに対比して**真の主因**を特定:

- **SU2 の流儀** (`CNSSolver::SetTau_Wall_WF`, `flow_diffusion.cpp:AddTauWall`): 壁関数 ON/OFF に依らず
  **常に u=0 strong Dirichlet**。壁関数のときは「壁ノード↔第一内部点エッジの**解像粘性応力テンソルを
  スカラー τ_w/|τ_resolved| で再スケール**」して、壁エッジの運動量流束がモデル τ_w を運ぶ
  (`tau[i][j] *= τ_w/|τ_tangent|`, `xor` ガードで片端のみ壁ノードのエッジに作用)。Neumann 流束加算でも
  速度自由化でもない。OFF は `Tau_Wall=-1` で再スケールを切るだけ。
- **forge cell**: `viscousFlux_wall_d`(wallTreatment==1) がモデル τ_w を**内部セル ic** に課す → 正しい。
- **forge node (欠落)**: モデル τ_w は壁境界半割面 (W→ghost) → **壁ノード W (Dirichlet で破棄)** へ。W-I 内部
  エッジは**生の解像応力 μ_total·u_I/dcc のまま**で τ_w を一切知らない。**壁ノード残差は両者とも捨てるが
  (ユーザー確認)、SU2/cell は τ_w を「届く面」に乗せ、node だけ「捨てられる面」に乗せていた** (ic=壁ノード
  トラップの3例目)。

→ **修正 (採用, SU2 AddTauWall 厳密移植)**: node で per-node `Tau_Wall=ρu_τ²` を壁関数が格納し、
  `viscousFlux_d` で片端のみ壁ノードの W-I 内部双対面 (`Tau_Wall>0` の xor) の接線 traction を τ_w に
  再スケール (`scale=τ_w·area/|接線traction|`)。cell/非WF は `Tau_Wall≡-1` で無効 (ビット不変)。壁ノード
  残差は不変 (Dirichlet ゼロ化のまま; τ_w は W-I エッジ流束経由で内部ノード I に届く)。

**検証 (決定的, `run_0067_node_sst_tauwall_bp1p90`, Ps1.90)**: コア加速 **Mmax 1.69→1.92**、衝撃 **46→~171mm**
(cell142/SU2132 域へ; 上流固着・背圧逆応答が消失)、utau/ypls 物理 (y+→145)。cell SST 不変 (142/1.82,
`run_cell_tauwall_regr`)。**主病理解消**。残差: shock が cell/SU2 より ~30mm 下流 (node Mmax1.92>cell1.82=BL
やや薄=τ_w やや過小付与) は別途リファイン (τ_w 大きさ・代表点 ρ・y の精査)。図 `cmp_node_addtauwall.png`。

## 4. 検証チェックリスト (修正後)

- `utau`/`ypls`/`wf_pk` が node で**非ゼロ・物理的**になる (現状 0)。
- `mut/mu` の異常なコア侵入が改善。
- 中心線 Mach の上流加速が cell/SU2 に近づく (M→~1.82)。
- shock front が下流へ移動し SU2(132)/cell(142) に近づく。
- `check_convergence.py` の VERDICT を引用 (case36 は plateau→「準定常/shock 漸近」と表現、「収束」と書かない)。
- 残差 CSV + 最終 `res_*.h5` で NaN/Inf チェック。
- run 作成時は case README「計算 run 一覧」を同期。
- 修正前後を**同条件**で比較 (背圧・mesh・convMethod 一致; 背圧不一致 1.84 vs 1.90 の交絡に注意)。

## 変更ログ

- 2026-06-24: 診断確定 (コード + utau/ypls=0 実測) と設計論点 (ユーザー memo) を記録。実装は SU2 調査結果待ちで
  保留 (status=blocked)。case36 検証 run: `run_node_sst_bp1p90_matched` (Ps1.90 matched, shock→46mm 漸近)。
- 2026-06-24: SU2 調査完了 (§3, Normal_Neighbor 流) → 方針 A 確定。`compute_wall_friction_sst_d` を cell/node
  分岐実装 (node: 壁内向き法線 cos 最大の代表内部点 I から Ut/rho/nu/y を一貫取得; cell ビット不変)。隔離
  worktree で build clean・検証。**結果**: `utau`/`ypls` が node で nonzero・物理 (y+~98) に修正 (実測
  nonzero_frac 0→1)、peak μt/μ 91k→57k 是正、case36 衝撃が **46→~65mm と +20mm 下流改善** (`run_0054_node_sst_wffix_bp1p90`,
  Ps1.90, shock 漸近 ~65mm, ただし `check_convergence` は plateau=未収束)。場は wt=1 node のみ変化、cell/flat-plate(wt=0) 不変。
  **ただし case36 はまだ cell/SU2(132-142mm) に届かず (65mm, Mmax 1.69 vs 1.82)**。残差主因は壁関数でなく
  **median-dual のコア・チェッカーボード** (x=50mm 半径分布で軸上コア Mach が波打ち沈む 1.62-1.68; μt は cell
  以下なので BL ブロッケージでない) と判明 → `gradLSQ`/`nodeMidpointFx` 検証は別途。壁関数 fix 自体は独立バグ
  修正として commit。

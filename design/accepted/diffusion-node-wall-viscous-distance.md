# node-centered 粘性壁 flux の法線距離修正 (mirror-ghost 退化対策)

## メタ

- **area**: `diffusion`
- **status**: `done` (§11 node 壁摩擦応力 twall=内部双対面集約; §1-10 の mirror-ghost 法線距離は別軸で保留)
- **related_docs**:
  - `docs/diffusion/theory.md`
  - `docs/diffusion/implementation.md`
  - `docs/discretization/implementation.md`
- **created**: `2026-06-20`
- **owner**: `未定`

## 1. 目的

case 29 (`29.bell_vs_conical`) の **node-centered viscous** 解で壁近傍に出る振動/implicit
不安定を解消する。Euler の node 解はきれいで、cell viscous も安定なのに node viscous だけ
壁近傍が振動する。原因は粘性壁 flux `viscousFlux_wall_d` の**法線勾配評価が mirror ghost
距離 `dcc = |cc_ghost - cc|` に依存している**ことにある。

## 2. 切り分け (コード確認結果)

- ghost 中心生成 ([`mesh.cpp` L450-456](../../solver_density_cuda/mesh/mesh.cpp#L450)):
  `cc_ghost = cc + 2((pc-cc)·n) n`。
- **cell**: cell 中心は壁面から法線距離 $d_n=(pc-cc)\cdot n$ だけ内側 → `dcc = 2 d_n`。
  mirror (`Ux[ig]=-Ux[ic]`) と合わせ法線項は $(U_{wall}-U_{cell})/d_n$ ($U_{wall}=0$) に一致。正しい。
- **node (median-dual)**: 境界 CV 代表点 (=ノード) も半割面重心 `pc` も**壁面上**にある →
  $(pc-cc)\cdot n \approx 0$ で `cc_ghost ≈ cc`、**`dcc ≈ 0`** に退化。法線項
  `(Ux[ig]-Ux[ic])/dcc` が $0/0$ 的に爆発・偶奇振動。
- 勾配側 ([`calcGradient_d.cu` L719-750](../../solver_density_cuda/cuda_forge/calcGradient_d.cu#L719))
  は既にこの退化を認識し、壁半割面を cellgather から外して **bvar 弱形式**で閉じている。
  → 「境界処理が無い」のではなく、**粘性壁 flux の法線項だけが mirror ghost に取り残されている**。

### 2.1 診断実測 (case 29 node viscous, wall physID 3, 272 半割面, `FORGE_VISC_WALL_DIAG=1`)

| 量 | min | max | 備考 |
| --- | --- | --- | --- |
| `dn`=(pc-cc)·n (符号付) | -3.06e-4 | 3.05e-4 | **272 面中 133 面が負** (符号反転) |
| `|dn|` | 1.05e-8 | 3.06e-4 | **4〜5 桁ばらつく**・最小は ~0 |
| `dcc`=\|cc_ghost-cc\| | 2.13e-8 | 6.13e-4 | 構成上 `=2|dn|` |
| `dcc/(2|dn|)` | 0.998 | 1.009 | **常に ≈1** (cell の関係式は node でも成立) |
| 接線オフセット \|(pc-cc)-dn·n\| | — | 2.58e-3 | `pc-cc` の**主成分**(壁ノードは壁面上) |

→ 当初想定の「`dcc≈0` 一律退化」ではなく、より病的: `pc-cc` は接線成分が支配的 (最大 2.58e-3)
で、法線成分 `dn` はその微小・符号反転する残差にすぎない。`dcc=2|dn|` がこの残差なので、法線項が
掛ける `1/dcc` は隣接壁ノード間で **4〜5 桁** 変動する。これが近壁の偶奇振動・implicit 剛性の正体。
第一候補だった `(Uxb-Ux)/dn` も同じ `dn` を使うため救済にならない (採用せず)。

## 3. 設計方針

`cfg.discretization=="node"` かつ `nodeWallViscGradFlux==1` (既定) のとき、壁法線項を
**node セル勾配ベース** $\mu\,\nabla u_i\cdot\mathbf{S}$ に置換する。$\nabla u_i$ は calcGradient で
bvar 弱形式により壁面で閉じているため退化距離に依存しない。

- 法線項: `mu*((U[ig]-U[ic])/dcc)*sss` → `mu*(dU_i/dx*sx + dU_i/dy*sy + dU_i/dz*sz)`。
- 熱流束: `tc*((Ts[ig]-Ts[ic])/dcc)*sss` → `tc*(dTdx*sx+dTdy*sy+dTdz*sz)` (断熱壁では ∇T≈0)。
- **転置項・体積粘性項は cell/node 共通で不変** (もともと node セル勾配にフル `S` を内積する形)。
- cell モードは `nodeWallGradFlux=0` を渡し従来挙動をビット不変で維持。
- 第一候補だった「dn と Dirichlet 距離 `(Uxb-Ux[ic])/dn`」は、node では `dn≈0` のため
  同様に退化する (診断で確認) ので採らない。`dcc` 除算が無くなるので NaN floor は不要。

## 4. 切替フラグ

- `mesh.nodeWallViscGradFlux` (int, 既定 1)。0 で旧 mirror-ghost 挙動 (A/B 比較用)。
- 診断 env `FORGE_VISC_WALL_DIAG=1`: 壁半割面の $d_n$, `dcc`, `dcc/(2d_n)`, 接線オフセットを
  1 回だけ集計表示 (`viscousWallDiag_d`)。

## 5. 検証

- ベース (現行 mirror-ghost): `case/29.bell_vs_conical/run_dual_visc_conical_node_expl_long/`
- 参照: Euler node `run_dual_eul_conical_node_m3/`, cell viscous `run_dual_visc_conical_cell_m3/`
- 新 run: `case/29.bell_vs_conical/run_NNNN_*` (explicit RK3, cfl=0.1, convMethod=0, bndFirstOrder=1)。
- 指標: 壁近傍 Ux/Uy/P/dUxdy/twall の x 方向振動低減、残差全列低下、`P≤Pt`/`ro>0`/`T>0`、
  出口リップ近傍の局所 NaN/極端 wall shear 消失、推力/出口プロファイルが Euler/cell/SU2 から
  大きく崩れないこと。

## 6. 検証結果と所見 (2026-06-20)

case 29 node viscous (explicit RK3, cfl=0.1, convMethod=0, bndFirstOrder=1) で A/B 検証。
run path: `case/29.bell_vs_conical/run_0019〜run_0025`。

### 6.1 診断 (確定)

- `FORGE_VISC_WALL_DIAG=1` で壁半割面 272 を集計 → §2.1 の通り `|dn|∈[1e-8,3e-4]` が 4〜5 桁
  ばらつき符号反転。`dcc=2|dn|` がそのまま振動源。**仮説どおり退化を実証**。

### 6.2 修正候補と結果 (すべて発散 or 無効果)

| run | 方式 | 結果 |
| --- | --- | --- |
| run_0020 | 旧 mirror-ghost (`flag=0`) | **安定**。twall_x TV/range=14.8, 符号反転 172, スパイク -1.77e4。T max 2967K(>Tt 1500K), 中心線 Mach は Euler/cell と一致(バルク健全)。基準。 |
| run_0019 | 純勾配 ∇u·S | **発散** (step ~13250, roe 先頭)。断熱壁から ∇T·S で熱漏れ。 |
| run_0021 | 純勾配 + 断熱熱流束=0 | step20k 安定だが **壁せん断 0.37Pa (cell ~4840Pa の 1e4 倍過小)** = no-slip 喪失(ほぼ slip 壁)。120k で別途発散。 |
| run_0022 | 滑らか距離 2·vol/sss | **発散** (step ~9750, roe 先頭)。ドラッグ弱すぎ。 |
| run_0024 | 滑らか距離 + 転置項除去 | **発散** (step ~9760)。転置項は無関係と確認。 |
| run_0023 | nodeWallDirichlet (残差射影) | **発散** (step ~22700, roe 先頭)。タスク既知の「壁全域発散」を再現。 |
| run_0025 | floor `max(dcc, 0.05·vol/sss)` | **安定だが twall ノイズ残存** (TV/range=15.8, 反転 180)。c 小=安定だが元同様ノイズ。 |
| run_0026 | クリーン Dirichlet (毎ステージ enforceWallNoSlip: roe から KE 除去 + roU=0 + P/T を内部エネルギーから再計算) | **全域崩壊** (step<10000)。残差射影のみ版 (run_0023) と同様に発散し、エネルギー一貫処理を入れても改善せず。 |
| run_0027 | クリーン Dirichlet + **壁ノード全保存量残差を凍結** (update から除外: res_ro/roU/roe=0) | **発散が ~3 倍遅延** (発散 ~22.6k 行 → ~68.5k 行)。壁ノード (特に壁∩出口コーナー) の**質量/エネルギー積算が主要な発散ドライバ**と確認 (ユーザー仮説の裏付け)。ただし完全凍結でも残る緩やかな第二不安定があり、全壁凍結は壁物理を固定するため正解でない。 |

### 6.2.1 コーナー (壁∩出口) が主因 (ユーザー指摘・検証)

> **【2026-06-24 訂正】下記の「質量/エネルギーが溜まる (accumulation/outflow できない)」という機構説明は
> §9.8 の細分出力 build-up + slip テストで否定された。** 実体は壁∩出口コーナーでの**出口静圧 BC の数値不安定
> (成長する圧力振動)** で、質量蓄積でも no-slip BL でもない (slip 壁でも同所で同様に発散)。詳細・訂正は §9.8。
> run_0027 の「全保存量凍結で3倍遅延」も、凍結が振動を一時鈍らせただけと再解釈する。

壁ノードが出口と一致するコーナーでは、壁優先所有で `wall_flag=1` のため運動量は射影されるが、
質量/エネルギー残差は射影されず、かつコーナーの出口方向 半割面が**壁ミラーゴースト (流出しない)**
になるため、質量/エネルギーが溜まり **roe が先頭発散**する。run_0027 で全保存量を凍結すると発散が
3 倍遅延したことがこれを裏付ける。これは既知の**「壁∩出口で矛盾ゴースト 2 個 → exit-lip 発散
(case/29 viscous node)」** (plan `discretization-node-boundary-ghostless.md`) と同一の問題で、
正攻法は同 plan Phase 2 の**コーナー面ごと BC** (コーナーの出口面は出口流束、壁面は no-slip) である。

### 6.3 結論

- **強 no-slip ドラッグ (小 dcc) が安定性に必須かつ twall 振動の源**であり、`viscousFlux_wall_d`
  内の局所的な距離操作では両立できない。弱めると (勾配/滑らか距離/floor 大/Dirichlet) no-slip を
  失い **energy が先頭発散**、強いままだと振動が残る。
- **バルク解は健全** (中心線 Mach が Euler/cell と一致)。問題は近壁 twall 出力と (implicit/高 CFL での) 安定性に局在。
- node 壁ノードは壁面上に乗るため、壁せん断は**壁ノード単体では作れず内部双対面 (壁ノード↔内部ノード)
  が担う**。よって正攻法は「壁ノード u=0 を一貫した質量/エネルギー処理つきで固定」か「内部隣接ノード
  距離を使うコンパクト法線作用素」で、いずれも `viscousFlux_wall_d` の距離一行では収まらない**アーキ
  テクチャ変更**。floor が物理距離 (vol/sss) でも発散することから、壁ノードの粘性 timestep/スペクトル
  半径の見積りも疑わしい (別軸の調査候補)。

→ 既定は `nodeWallViscGradFlux=0` (元挙動を保持・回帰なし)。診断カーネルとフラグ/floor は切り分け・
  将来の正攻法用に残置。正攻法はスコープ要相談 (status を done にせず保留)。

## 7. SU2 比較に基づく再設計 (2026-06-20, `.external/su2-src`)

SU2 は forge node と同じ**頂点中心 median-dual**。その no-slip / 粘性壁実装を確認した結論:

**SU2 の流儀 (`CNSSolver::BC_HeatFlux_Wall_Generic`, `CIntegration::Space_Integration`,
`CAvgGrad_Base::CorrectGradient`, `CEulerVariable::SetVelocity_Old`):**

1. **ゴーストを使わない**。各境界マーカーが自分の半割面を積分しノードに加算 (弱形式)。
   **コーナーノードは複数マーカーの頂点リストに入る**ため、出口の流出 (質量/エネルギー) と壁の
   no-slip を**両方**受け取る。
2. **BC は2パス**: ①弱 (inlet/outlet/farfield, 残差へ流束加算) → ②**壁を最後に**強制適用。
   壁が出口の後なのでコーナーで壁 no-slip が速度を上書きするが、出口の質量/エネルギー流出は残る。
3. **no-slip = 強 Dirichlet**: 壁ノードで運動量残差=0 (`LinSysRes(vel)=0`)・truncation error=0・
   (implicit) `Jacobian.DeleteValsRowi` で速度行を単位化。`SetVelocity_Old` は `Solution_Old=ρ·u_target`。
   ρu を凍結 (ρ 変化による微小 u ドリフトは許容)。**壁ミラーゴースト粘性 flux は存在しない**。
   壁せん断は内部エッジが担う。
4. **断熱壁 = エネルギー弱流束 0**: `Res_Visc = Wall_HeatFlux·Area = 0`。壁面熱流束を一切加えない。
5. **内部粘性 flux = 平均勾配 + CorrectGradient** (エッジ整合補正: 面平均勾配のエッジ射影を
   コンパクト差分 (φ_j-φ_i)/|r_ij| に一致させる over-relaxed)。forge は同等の over-relaxed を保有。

**forge が発散する真因 (SU2 との差分)**: forge node は**ゴースト+単一所有**でコーナー (壁∩出口) が
壁ミラー (非流出) になり質量/エネルギーが溜まる。SU2 はコーナーが出口面を持つので流出する。
→ 「壁ノード全凍結で発散3倍遅延」(§6.2.1 run_0027) と整合。

### 7.1 再設計 (SU2 丸パクリ方針)

本問題は「粘性壁 flux の距離」ではなく**境界クロージャ (ゴースト→弱形式マルチマーカー)** が本体。
よって本 plan の主眼を [discretization-node-boundary-ghostless.md](../active/discretization-node-boundary-ghostless.md)
Phase 2 と統合し、次を実装する:

1. **node 境界をゴーストレス弱形式に**: 各マーカーが半割面流束をノードに加算。コーナーは所有を
   単一化せず**全マーカーから寄与**を受ける (出口辺=流出流束, 壁辺=no-slip)。
2. **BC 適用順を 弱→壁 に**: 対流境界 (inlet/outlet) を先に残差へ加算し、**壁 no-slip を最後に**
   強制 (運動量残差=0)。コーナーで壁が速度を上書き、出口の質量/エネルギー流出は残す。
3. **壁ミラーゴースト粘性 flux を廃止** (`viscousFlux_wall_d` を node では呼ばない)。せん断は内部
   双対面 (壁ノード u=0 ↔ 内部ノード) が担う。断熱はエネルギー弱流束 0、等温は規定温度で熱流束。
4. **強 no-slip の ρ ドリフト対策**: SU2 `SetVelocity_Old` 同様、更新後に壁ノードの ρu を ρ·0=0 へ
   再設定 (explicit でも u=0 を厳密化)。本 plan の `enforceWallNoSlip` (KE 除去+roU=0 毎ステージ)
   はこの一部。残差ゼロだけでは ρ ドリフトで u≠0 になるため状態再設定が必須。
5. implicit (block-DPLUR) では SU2 同様 Jacobian の速度行を単位化して強 Dirichlet と整合させる。

twall 出力ノイズ (dcc 退化) は副次問題。壁ミラーゴースト粘性 flux を廃止すれば twall も内部エッジ
勾配ベースで滑らかに出せる (本 plan の診断フラグはその検証に流用)。

## 8. 具体的実装仕様 (SU2 化, Explore マップ確定 2026-06-20)

Explore で現状を確定: **forge node 境界は既にほぼ弱形式 (bvar ベース)**。ghost は非壁対流とミラー
勾配にのみ使用。勾配 (`calcGradient_b_d`)・壁対流 (`convectiveFlux_boundary_d`, pressure-only)・粘性壁は
すべて per-bcond の bvar/半割面を消費する。**唯一の構造的欠陥は gmshReader の半割面 emit が単一所有**で
ある点 (`buildMedianDual`, gmshReader.hpp 行 1482-1536):

```cpp
const int ow = ownerBc[N];          // ← 単一所有: ノード N の半割面を「N の所有 bcond」へ
auto& h = halfByOwner[ow][N];       //   コーナー(壁∩出口)は壁所有 → 出口の流出面を持たない
```

### 8.1 変更点 (最小・SU2 整合)

1. **gmshReader.hpp `buildMedianDual` emit をマルチマーカー化** (本丸):
   - `halfByOwner[ownerBc[N]][N]` → `halfByOwner[ib][N]` (辺の属する bcond `ib` へ集約)。
   - これでコーナーノード N は壁 bcond の half-face (壁辺由来) と出口 bcond の half-face (出口辺由来) の
     **両方**を得る。各 bcond の `dualBcondNodes`/`iPlanes`/bvar 配列に N が両方入る。
   - `ownerBc` は wall_flag 構築 (どのノードが壁か) には引き続き使用。
   - 注意: コーナーは 2 つの境界 plane (壁・出口) を持ち、各々 ghost を生成する (mesh.cpp は plane ごとに
     ghost を作るので自動)。勾配弱形式は両 half-face の bvar·S を加算 (正しい)。

2. **mesh.cpp**: 変更最小。plane ごとに ghost 生成する既存ロジックがマルチマーカーでもそのまま動く
   (コーナーが 2 ghost)。`wall_flag_d` は壁 bcond の iCells から構築 (コーナー含む, 既存)。

3. **convectiveFlux_d.cu**: `convPlaneBound` は壁 plane を除外し、壁 plane は `convectiveFlux_boundary_d`
   (pressure-only) で処理 (既存)。コーナーの出口 plane は非壁なので主ループ+ghost で**流出流束**が入る
   (既存ロジックがマルチマーカーで自動的にコーナーへ流出を与える)。← これが発散の直接修正。

4. **BC 適用順 + 強 no-slip を最後に** (boundaryCond.cpp `applyBconds` / main.cpp `assembleResidual`):
   - 残差積算後、**最後に**壁 no-slip を強制: `enforceWallNoSlip` (roe から KE 除去 + roU=0) を毎ステージ
     更新後に適用し、`zeroWallMomentumResidual` で運動量残差も 0。コーナーで壁が速度を上書きしつつ、
     出口由来の質量/エネルギー流出 (res_ro/res_roe) は残す (SU2 の2パス順序と等価)。
   - SU2 `SetVelocity_Old` 相当: ρ ドリフト対策に更新後 roU=0 を再設定 (enforceWallNoSlip が担う)。

5. **壁ミラーゴースト粘性 flux を node で無効化**: `viscousFlux_wall_d` を node では呼ばない (せん断は
   内部双対面が担う)。断熱はエネルギー寄与 0。twall 出力は内部エッジ勾配ベースに変更 (診断フラグ流用)。

6. **implicit (block-DPLUR)**: 壁ノードの速度行を単位化 (SU2 `DeleteValsRowi` 相当) して強 Dirichlet と
   整合 (Phase 2 後半)。

### 8.1.1 実装試行の結果 (2026-06-20, run_0028/0029)

§8.1-1 (emit を `ownerBc[N]`→`ib` でマルチマーカー化) を実装し conical.msh を再変換:
- **幾何は健全**: `buildMedianDual` の体積閉性 relErr=8e-8、境界閉性 normalized=4.7e-6、各マーカーの
  半割面積=primalArea。コーナーが壁・出口の両半割面を正しく分割保持。
- **しかしソルバが step 0 で発散** (roe NaN、nodeWallDirichlet=0 の元壁処理でも同じ)。
  → **emit 変更だけでは不十分**。ソルバ側の境界 plane/ghost/CV ループが「ノードあたり単一 bcond /
  単一境界 plane / 単一 ghost」を暗黙前提しており、コーナーが 2 境界 plane + 2 ghost を持つと破綻する
  (setMeshMap_d の plane 順序・対流主ループの ghost 経路・bvar 配列・CV 体積閉包のいずれか)。
- 結論: **full ghostless は emit のマルチマーカー化に加え、ソルバ境界ループのマルチ bcond 対応**
  (1 ノードが複数境界 plane/ghost を持てるデータ構造とループ) が必須。これは
  [discretization-node-boundary-ghostless.md](../active/discretization-node-boundary-ghostless.md) Phase 2 の本体。
- 共有作業ツリーが並行セッションで頻繁にビルド不能 (`boundaryCond_d.cu`/`calcWallDistance_kdtree`) に
  なり、反復デバッグが困難。**専有ツリー (worktree or 別クローン) での一括実装+検証**を推奨。
- emit のマルチマーカー変更は単体ではソルバを壊すため working tree から撤回した (再実装時は
  ソルバ側マルチ bcond 対応とセットで入れる)。

### 8.1.2 詳細切り分け (worktree, 2026-06-20)

隔離 worktree (HEAD + multi-marker emit + HighFive) で再ビルド・再変換し詳細に切り分けた:
- **multi-marker mesh + 元の壁処理 (nodeWallDirichlet=0) でも発散**: 先頭 5 step は max cfl≈0.58 で正常
  だが、30000 step かけて **max cfl→0 で場が崩壊**。残差 rms は step 0 outer_begin から NaN
  (roUx/roe)、ただし **detectNaN (state 検査) は未発火** → NaN は実セルでなく**余分なコーナー ghost**
  側に出ている (rms 集計が ghost を含む)。
- ghost 計上自体は整合 (`nCells_ghst=nBPlanes`=878, multi-marker のコーナー重複を含む)。
  問題は**コーナーが持つ 2 つの境界 plane/ghost の state 設定・flux ルーティング**で、ソルバの境界
  ループ (対流主ループ ghost 経路・BC kernel の ghost 書込み・rms 集計範囲) が「1 ノード=1 境界 plane/
  ghost」を暗黙前提しているため、コーナー ghost が正しく扱われず残差 NaN + 緩慢な崩壊を招く。

### 8.1.2.1 【重要・訂正】multi-marker は無実、fresh conical.msh が独立に発散 (2026-06-20)

worktree で **single-owner (baseline, 私の変更なし) でも同じ fresh conical.msh で同一に発散** することを
確認 (step0 rms NaN + cfl→0)。さらに **旧 nozzle.h5 (27472 ノード, 272 壁) は現 HEAD worktree binary で
正常** (cfl 0.50, rms_ro 2e-5 低下, NaN なし)。

→ **§8.1.1 の「multi-marker emit がソルバを壊す」は誤り**。真相:
1. **multi-marker emit は発散を起こさない** (single-owner と挙動同一)。安全。
2. **fresh conical.msh (34138 ノード, 細かい) が現 HEAD の node viscous で独立に不安定** (convert の体積/
   閉性チェックは通るのに step0 で 1 セル残差 NaN→緩慢崩壊)。これは粗い旧メッシュ (27472) では出ない
   **メッシュ解像度/品質依存の別問題**。run_0028/0029 の発散はこれに**交絡**していた。
3. よって corner/multi-marker 修正の検証には **動く再変換可能メッシュが必須**だが、唯一の再変換元
   conical.msh が壊れている (旧 27472 メッシュの素 .msh は不明/変化)。粗い conical.msh 再生成 (gmsh) か
   fresh mesh 不安定の原因究明が前提。
4. 旧 (動く) メッシュでの既知: 元壁=安定 (twall ノイズ, run_0020)、nodeWallDirichlet 単一所有=発散
   ~22k (run_0023, コーナー非流出)。multi-marker を旧メッシュに適用できれば ~22k 発散が直る可能性が
   高い (コーナー流出) が、旧メッシュは pre-converted h5 で再変換不可のため未検証。

### 8.1.2.2 粗メッシュでの検証 (run_0030, 2026-06-20)

gmsh で粗い conical (ny=40, 13858 ノード) を再生成し、動く baseline を確保:
- **粗 single-owner (run_dbg_coarse, Dirichlet なし)**: 安定 (cfl 0.1, rms_ro 1.16e-6)。
- **粗 multi-marker + nodeWallDirichlet=1 (run_0030)**: **40000 step 安定** (NaN なし, ro>0,
  **P max=4.003e6≈Pt** [旧メッシュ 4.78e6=19%超過から改善], T max=1827 [Tt 1500 を 22%超過, 旧 2967=
  98%超過から改善だが残存])。
- **ただし粗メッシュは single-owner と multi-marker で per-bcond ノード数が同一** (inlet41/outlet41/
  wall338/axis338)。すなわち**このメッシュでは wall∩出口コーナーが分割されず multi-marker が実質無効**。
  → run_0030 の安定性は「粗メッシュ + Dirichlet」由来で、**コーナー修正の純粋検証にはなっていない**。
- 含意: (1) nodeWallDirichlet (強 no-slip) は**メッシュ依存**で、粗メッシュでは安定だが旧 27472
  (run_0023) では ~22k 発散。(2) コーナー修正を検証するには **wall∩出口コーナーを持ち、かつ node
  viscous で安定なメッシュ**が必要だが、fine (コーナーあり) は独立に発散、coarse はコーナー分割が
  出ない。multi-marker emit が fine では +4 corner を作るが (874→878)、その fine 自体が発散するため
  検証不可。
- T overshoot 残存 (1827>Tt) は近壁/コーナーの粘性発熱の処理が未完なことを示す (壁ミラーゴースト粘性
  flux の廃止 + 内部エッジせん断への一本化が必要; §8.1.3 パスB)。

### 8.1.2.3 fresh-fine conical.msh 不安定の原因究明 (タスク1, run_dbg_fine, 2026-06-20)

fine conical.msh (ny=80, 34138 ノード) の step0 発散原因を特定:
- 発散は **rms_roUx=inf @ step0 outer_begin** = **初回粘性 flux 評価で Inf**。detectNaN(state)未発火と整合。
- **dual CV 体積分布**: min=3.46e-11, mean=9.24e-8 (**min/mean=3.7e-4**)。最小級セルは**全て同一 vol
  3.456e-11**(系統的, ランダムでない)で、位置は **x≈0.18–0.97mm, r≈10.0–10.16mm = スロート直下
  arc 部 (Rfd=0.382Rt) の近壁帯**。
- 機構: スロート直下の**強曲率 arc + 細かい radial** が小体積・歪み near-wall CV を作り、その内部双対面
  の **dcc が極小→ over-relaxed 粘性項 (Ux1-Ux0)/dcc·delta が Inf**。これは**壁 flux の dcc 退化と同類
  の現象が内部スリバセルで起きたもの**。粗 (ny=40) はスリバを作らないので安定。
- **結論: fresh-fine の不安定はメッシュ品質 (スロート部スリバ) 起因で、ソルバ/コーナー/壁 flux 修正とは
  独立**。run_0028/0029 の発散はこれに交絡していた (§8.1.2.1)。

### 8.1.2.4 検証メッシュ確保のブロッカー (コーナー共有)

- 現 `make_nozzle.py` 生成メッシュは **single-owner と multi-marker で per-bcond 数が同一** (粗 ny=40 で
  inlet41/outlet41 不変) = **wall∩境界コーナーが分割されない (+0)**。一方 committed conical.msh は
  multi-marker で +4 corner (874→878) = コーナー共有あり。→ 現生成法ではコーナー検証メッシュが作れず、
  committed (コーナーあり) は throat スリバで発散。**「コーナー共有 かつ スリバ無し」のメッシュが未確保**。
- 次手: (a) make_nozzle.py の throat_dn arc 部の radial/axial spacing を調整しスリバを除去 (min/mean を
  ~1e-2 以上に)、(b) Physical Curve のコーナー共有を確認/修正、の双方を満たすメッシュを生成してから
  single+Dirichlet (発散?) vs multi+Dirichlet (安定?) でコーナー修正を実証する。

### 8.1.2.3.1 スリバ除去 (メッシュ修正, 2026-06-20, 完了)

§8.1.2.3 のスリバ原因を**メッシュ生成側で根治** (floor で隠さない方針 §8.1.2.5 に整合)。
- 原因の定量: 軸方向 dx が chamber 1.11mm / **throat_dn 0.0225mm (nx=44, 軸長~1mm)** / cone 0.42mm
  (nx=200)。throat_dn が cone の **~19 倍細かい**過剰分割で、スロート下流に薄い near-wall CV を作っていた。
- 修正 (`mesh/make_nozzle.py`): `NX_ZONE['throat_dn']` **44→10** (dx~0.099mm)、発散部の軸クラスタ比
  `--prog-x` 既定 **1.0→1.012** (cone をスロート側=細→出口=粗に段差レス接続)。`.msh` は gitignore
  生成物なので **make_nozzle.py がソース変更**。`mesh/conical.msh`/`bell.msh` を ny=80 で再生成。
- 結果: **vol min/mean 3.7e-4 → 4.8e-2 (130 倍均一化)**、負体積 0。**ny=80 fine conical が node viscous
  で安定** (NaN なし, cfl 0.1, rms_ro 低下) — 以前は step0 発散。スリバ起因の不安定は解消。
- 残: このメッシュが wall∩出口コーナーを共有するか (multi-marker 効果) は未確認。コーナー修正の純粋
  検証 (§8.1.2.4) はコーナー共有の確認後。

### 8.1.2.5 方針確定: dcc フロアは不採用 (2026-06-20, ユーザー指示)

**原則 (ユーザー指示): dcc フロアを物理 flux の通常計算に使わない。** 退化した極小 dcc が生じる時点で
メッシュ(または境界クロージャ)が間違っており、マスクして物理 flux を改変するのは誤り。直すべきは:
- **内部スリバ (§8.1.2.3)**: メッシュ生成側で除去 (throat 下流アークの spacing)。floor で隠さない。
- **壁 mirror-ghost の退化 dcc (§2.1)**: 退化した距離を作るクロージャ自体が誤り。**ゴーストレス弱形式**
  (SU2 流) に置換し、そもそも mirror-ghost 距離を使わない (=退化 dcc を生まない)。floor で隠さない。

→ 実装決定:
- `viscousFlux_wall_d` に入れた **dcc フロア (nodeWallViscGradFlux/nodeWallDistFloorCoef) は撤回**し、
  物理 flux は raw dcc のみで計算する元の形に戻した。config フラグは未使用 (deprecated, 既定で無影響)。
- 診断カーネル `viscousWallDiag_d` (env `FORGE_VISC_WALL_DIAG`, flux を変えず dcc/dn を計測) は温存。
- 根治は (a) スリバ無し+コーナー共有メッシュの生成、(b) ゴーストレス弱形式壁/境界クロージャ (§8.1.3 パスB)
  の2本に一本化する。

### 8.1.3 実装の2パス (次セッション向け)

full ghostless を完遂するには、emit マルチマーカー化に加えて次のどちらかが必要:

- **パスA (ghost 拡張)**: コーナーが複数境界 plane/ghost を持てるよう、(1) 各境界 plane の ghost state を
  対応 bcond kernel が確実に書く (corner の wall ghost=mirror, outlet ghost=外挿)、(2) rms 集計を実セル
  [0,nCells) に限定、(3) 対流主ループが corner の outlet plane を ghost で正しく流出させる。ghost を
  維持する分、既存コードからの差分は小さめだが corner の二重 ghost 整合が要注意。
- **パスB (真のゴーストレス, SU2 厳密)**: 境界対流を**ghost 無しの弱形式**に置換 (calcGradient_b_d と
  同型: 各半割面で interior 状態+bvar から流束を作りノード残差へ加算)。outlet/inlet も weak 化。ghost を
  全廃しコーナーは複数 bcond の weak 流束を自然に受ける。SU2 に最も忠実だが convectiveFlux の境界経路を
  全面書換え。

いずれも **専有ツリーで一括実装+段階検証** (まず Euler node で回帰なし→viscous node で発散解消→twall
平滑→中心線/出口が Euler/cell/SU2 一致) する。共有ツリーは並行セッションで頻繁にビルド不能になるため
不可。

### 8.2 リスクと進め方

- §8.1-1 (emit マルチマーカー化) は **node モード全体に影響**する破壊的変更。半端な適用は全 node ケースを
  壊すため、**安定した専有ツリーで一括実装+検証**する必要がある (現 main は他セッションの
  `boundaryCond_d.cu` 編集でビルド不能・共有汚染中)。
- 検証: case 29 node viscous (explicit RK3, cfl 0.1) で①NaN/発散解消 ②twall 平滑 ③中心線/出口が
  Euler/cell/SU2 一致、を確認。既存 Euler node (run_dual_eul_conical_node_m3) が回帰しないことも確認。
- 親 plan [discretization-node-boundary-ghostless.md](../active/discretization-node-boundary-ghostless.md) Phase 2 と統合。

## 9. 自走検証の結果 (2026-06-20, run_0030〜0034)

B (multi-marker emit + nodeWallDirichlet 強 no-slip) を node 壁の主軸に確定し、自走で (1)-(6) を検証。

### 9.1 コーナー修正 A/B (uniform 均一メッシュ, 30k)
- A `run_0031_corner_single_dir` (single-owner, コーナー非共有: outlet=80, wall∩outlet=0):
  解は正常・収束するが **出口リップ・コーナー1セルで max cfl=2.92e25 (一定・退化)**。
- B `run_0032_corner_multi_dir` (multi-marker, コーナー共有: outlet=81, wall∩outlet=1):
  **max cfl=0.1 正常**。100k で **check_convergence PASS (収束)**。
- → multi-marker は良いメッシュ上では**解を変えず、退化コーナー cfl/dt を正常化するロバスト性改善**
  (cfl 律速 dt/implicit で dt→0 を防ぐ潜在地雷を除去)。

### 9.2 他 node ケース回帰 (PASS)
- case 24 (laminar_channel_bl, node viscous): multi/single とも NaN なし・cfl 1.0・rms_ro 同等 (壊れない)。
- case 08 (bump, node Euler): multi/single とも NaN なし・cfl 0.5 (multi はむしろ残差低)。

### 9.3 壁寄せ (radial clustering) メッシュ
- **誤検出の訂正**: 当初 prog_r=0.3-0.7 で median-dual 破綻 (非正体積・非閉) と報告したが、これは
  **gmsh Progression が等比で ny=80 セルに対し比^79 ≈ 10⁻⁴¹ という異常グレーディングを作るため**で、
  node スキームの限界ではない。
- **正しい壁寄せ (Bump 0.4, case 24 同方式) は健全**: 変換 OK (非正0・閉性OK・vol min/mean 2.6e-2)、
  node viscous で安定 (`run_0033_wallclust_bump_node`, cfl 0.1, rms_ro 1.6e-6)。**壁寄せは node で問題なく動く**。

### 9.4 T overshoot は近壁「解像不足」アーチファクト
- uniform 均一メッシュ: Tmax=1689 (Tt 1500 を **+13%** 超過, 非物理)。cell viscous も +24%, 旧 node +98%。
- **壁寄せ (Bump 0.4, BL 解像) で Tmax=1508 (+0.6%, ほぼ解消)**。
- → T overshoot は scheme バグでなく**近壁の steep 層を均一粗格子が解像できず数値散逸で過熱**するもの。
  **壁寄せ (BL 解像) が物理的な正解**で、これで消える。cell でも同根 (要 BL 解像)。

### 9.5 CFL 上限 (壁寄せ Bump0.4 メッシュ, node vs cell)

| 解法 | node-B (multi+Dirichlet) | cell-center |
| --- | --- | --- |
| **explicit RK3** | **~0.7** (0.7安定/1.0発散) `run_0034_cfl*_wc` | **~0.7** (0.7安定/1.0発散step35) `run_0036_cell_expl_*` |
| **implicit (blockDPLUR, nInner5, warm-start)** | **~5-10** (5優/10marginal/20発散) `run_0035_imp_ws_*` | **~5-10** (5優/10marginal/20発散) `run_0036_cell_imp_*` |

- **CFL 上限は node と cell でほぼ同一** → node 化 (B) による CFL ロバスト性ペナルティなし。
- explicit ~0.7 は壁寄せ近壁セルの粘性 dt 制限が律速。implicit で ~7-14倍。
- **implicit は warm-start 必須**: 冷態 isentropic IC からは node/cell とも step ~10 で発散
  (理想 1D IC + 起動過渡が implicit 線形化を超える)。「isentropic IC → explicit で発達 → implicit」が実用手順。
- 全代表 run に流れ場 h5 出力済 (outStepInterval=2000)、nbad=0・P≤Pt・ro>0 を確認。
- cell も壁寄せで **T overshoot 解消 (Tmax 1499.9, +0%)** = §9.4 の解像不足説を再確認 (scheme 非依存)。

### 9.6 SU2-lam クロスチェック (PASS)
- 中心線 Mach が **SU2-lam と ~2-5% 一致** (出口: node-B 3.80 / node-WC 3.77 / SU2 3.71 / Euler 3.85)。
- **no-slip BL 正常**: node viscous は壁セル |U|=0 (欠損あり)、Euler は壁でも ~1450 (slip)。
- 推力: node-B 1889 N, node-WC 1907 N (~1% 一致, 妥当)。

### 9.7 結論
- **B (multi-marker + Dirichlet) を node 壁の主軸に採用**。コーナー退化 cfl を解消し、SU2-lam と一致、
  no-slip BL を正しく形成。回帰なし。
- **T overshoot の根治 = 壁寄せ (BL 解像)**。scheme 側の追加修正は不要。
- 残: WC の長時間収束、implicit CFL、既定値 (nodeWallDirichlet を node 既定 ON 化) の確定、cell 側の
  BL 解像での T overshoot 再確認。

### 9.8 【症状再発】case36 実コーナーで roe 先頭発散が再発 (2026-06-24)

case36 node 層流 1次 (Ps1.70, `run_node_lam1st_bp1p70`) が **step7185 で発散**。detectNaN ダンプ
(`res_nan_7186.h5`) で局在化: **壁∩出口コーナー** (x=673-690mm=出口端, y=21.9-22.3mm=壁, wall_dist≈0)
の 30 セルで ro/roe/roUx/roUy が NaN (P/T 有限)、max cfl→4.3e11。**§6.2.1 の roe 先頭発散と同一**。

**重要 — マルチマーカーは効いているのに再発**: case36 node mesh は**コーナー2ノード (x=690,y=±22.27) が
outlet・wall の両 bcond に所属** (multi-marker emit `ow=ib` は動作)。それでも発散 → **§8.1.1 の予言通り
「emit だけでは不十分、ソルバ側のマルチ bcond 対応が必須」**。§9.1 の検証は**実コーナーの無い均一/粗メッシュ**
(§8.1.2.4: per-bcond 数が single/multi で同一=コーナー非共有) だったため、**実コーナーでの破綻は未検出だった**。
→ コーナー修正は emit のみ適用で**ソルバ側 (コーナーの出口 outflow) は未完**のまま、case36 の実コーナーで再発。

**【機構訂正 2026-06-24, 細分出力 build-up + slip テストで確定】**: §6.2.1/上記の「質量/エネルギー蓄積
(outflow できず溜まる)」は**誤り**。実体は**壁∩出口コーナーの成長する圧力振動 (数値不安定)**:
- 細分出力 (50step毎, `run_nodelam_fineout`): コーナー P が **規定 Ps=1700 を中心に振幅増大して振動** (std
  30→333 kPa で約10倍)、ρ も 21↔6.6 と振動し発散。**蓄積 (単調増) ではない**。コーナー速度は u=0 のまま
  (壁面に沿う流量は無い=ユーザー指摘どおり)。先に near-wall 出口ノードが 1700→1000 へ落ち、続いてコーナーが
  振動成長。図 `corner_divergence_buildup.png`。
- **slip テスト (`run_nodelam_slipwall_cfl1`, 壁を slip 化)**: それでも**ほぼ同じ step3176 で発散・コーナー P
  std309kPa で振動** → **no-slip BL は無関係**。トリガは**出口静圧 BC が壁と交わるコーナーでの Ps 規定**そのもの
  (slip/no-slip 不問)。
- → **真因 = 壁∩出口コーナーでの outlet_statPress (Ps 規定) の数値不安定**。高 CFL (陰解法) で圧力振動が成長、
  低 CFL で減衰。質量流束/蓄積/no-slip BL はいずれも犯人でない。

**CFL 律速 (実測, run_nodelam_*)**: cfl_pseudo=1.0 で step7185 発散、**0.3/0.1/陽解法(0.2) は 15000 step 安定**
(NaN なし)。当面 node viscous/laminar は **cfl_pseudo≤0.3** で運用可。安定後の解は過散逸 (Mmax1.40, 衝撃入口
寄り)。3者比較 `cmp_laminar_3way.png`: SU2 clean 155mm/1.92, cell 入口+train/1.64, node cfl0.3 入口/1.40。

**根治方針 (修正)**: 「コーナーの outflow を直す」ではなく**「壁∩出口コーナーでの出口静圧 BC を数値安定化」**。
候補: (a) コーナーで出口 Ps 規定を**弱形式化＋緩和** (圧力の急な強制を避ける)、(b) 出口を下流へ延長して壁∩
出口コーナーを除く (幾何)、(c) コーナーノードの出口/壁 BC の優先順位・整合を見直す。**CFL 低減 (cfl≤0.3) は
対症療法だが当面有効**。SU2 が同条件で振動しないのは固定CFL＋弱形式出口の差 (本 plan の SU2 対比と整合)。

## 10. スロート dPdx スプリアスピーク = メッシュ段差起因 (2026-06-20, ユーザー指摘, run_0037)

ユーザーがスロート付近の dPdx 不連続を指摘。軸中心ラインで数値検証:
- **音速点 (M=1) は x≈1.7mm、一方 dPdx ピークは x≈0.10mm** で一致せず。x≈0.10 は throat_up→throat_dn
  ゾーン境界 (x=0) 直後で M はまだ ~0.87。
- 当初「forge dPdx が FD(P) と一致するから物理的」と判断したが**誤り**。勾配カーネルは健全 (FD と一致)
  だが、**P 場自体が dx 段差で kink を持つ** (対流スキームが dx 不連続で偽の圧力勾配を作る)。
- **決定的検証 (run_0037)**: throat_up/throat_dn の dx を揃える (throat_up=30→60, throat_dn=10→8, 共に
  dx~0.12mm, 段差解消) と、**x≈0.10 の dPdx 偽ピーク (-2.76e8) が消滅**し dPdx は throat 域で平滑
  (-1.74→-1.81e8) に。ピークは x≈1.1 (M~0.95, 音速点側) へ移り大きさも -2.06e8 に低下。
  → **dPdx 偽ピークは throat の dx 段差起因と確定** (ユーザー指摘が正しい)。
- 修正: `make_nozzle.py` NX_ZONE `throat_up 30→60, throat_dn 10→8` (dx~0.12 で揃える)。`mesh/conical.msh`/
  `bell.msh` 再生成 (26892 ノード, throat dx 0.118→0.125 連続, vol min/mean 5.6e-2)。
- 残: x≈1.1 (Nc, throat_dn→cone 近傍) の弱いピーク (-2.06e8) が物理 (transonic) か微小残差かは要追確認。

## 11. 専用カーネル `viscousFlux_wall_node_d` (壁隣接内部双対面のせん断, 2026-06-24)

ユーザー要望: node モード専用の壁粘性カーネルを新作し、(a) 壁ノード→内部ノードの map 上でループ、
(b) せん断 flux を**内部ノードの残差のみ**に加算 (壁ノードには入れない)、(c) **over-relaxed** を
ちゃんと組む、(d) 面勾配は壁ノードと内部ノードの**平均 (dual mesh なので 1/2)**、(e) ghost (`ig`) は
使わず壁値は `Uxb` (no-slip=0) を使う。

### 11.1 `ic` の意味 (確認・確定)

`ic = bplane_cell[ib]` は半割面所有 CV のインデックス。**cell**: 壁隣接の内部セル (`planes[pln].iCells[0]`,
[`gmshReader.hpp` makeMesh])。**node**: **壁ノードそのもの** (`bconds[ib].iCells = dualBcondNodes[ib]`,
[`gmshReader.hpp:1718` replacePrimalWithDual])。→ **node の `ic` は内部ノードではなく壁面上の壁ノード**。
ユーザーの不安は的中。よって `dUxdx[ic]`=壁ノード勾配、`Ux[ig]`=壁ミラーゴーストで、node では不適切。

### 11.2 【訂正】二重計上ではない — 残差は内部カーネルが担い、壁カーネルの実役は twall/y+

当初 §11 草案で「新カーネルが I に加算すると二重計上」と書いたが**不正確だった (ユーザー指摘で訂正)**。
正しい切り分け (コード確認):

- median-dual で**壁ノード W ↔ 内部ノード I を結ぶ primal edge は「内部双対面」**であり、`viscousFlux_d`
  ([`viscousFlux_d.cu:124-167`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu#L124)) が全内部双対面を
  **保存形** (`+flux→ic0, -flux→ic1`) で処理。W-I 面では `Ux[W]≈0` (Dirichlet) なので `(Ux[I]-0)/dcc` が
  no-slip 法線せん断そのもの → **内部ノード I の運動量せん断は既に完全に計算済み**。
- 一方 `viscousFlux_wall_d` (境界半割面) は `res_*[ic]` の **`ic`=壁ノード W** に加算する。node で
  `nodeWallDirichlet=1` のとき壁ノード運動量残差は `zeroWallMomentumResidual_d` /
  `timeIntegration_d.cu:1256` (壁3行 decouple) / `enforceWallNoSlip_d` で**破棄**される。
  → 壁カーネルの残差加算は **I への二重計上ではなく「壁ノードへ捨てている」だけ** (無害だが無効)。
- 実消費を確認: **`twall_x/y/z_b` は `viscousFlux_wall_d` が書くのみで他に消費なし=壁摩擦応力の出力**。
  **`ypls_b` は `ransWallFunction_d.cu:55,150` が消費=SST 壁関数へ供給**。

→ **node での `viscousFlux_wall_d` の実効的役割は `twall` (壁摩擦応力) と `y+` の算出だけ**。そしてそれが
**退化 dcc (§2.1) と壁ノード勾配 `∇U[W]·S` で不正確** (§6.2 で twall ノイズ・符号反転を実測)。これが
真の問題で、ユーザーの「壁の摩擦応力を出すためにやってる」が的確。

### 11.3 設計 (twall/y+ 専用算出, 残差は不変)

**狙い**: node の壁摩擦応力 `twall` と `y+` を、退化 dcc/壁ノード勾配ではなく**壁ノードに接続する内部双対面
の集約**で正しく出す。**運動量・エネルギー残差には一切触らない** (内部カーネルが担うため)。
→ `viscousFlux_d` のマスクは**不要**、二重計上の懸念は無し、安定性への破壊的影響も無し (草案より安全)。

物理的定義 (保存的): 壁ノード W の CV に内部双対面から入る粘性運動量 flux の総和 = 流体が W の CV に及ぼす
粘性力。W は Dirichlet (u=0) なので、これが**壁が流体に及ぼすせん断力**に等しい。よって
`twall(W) = (Σ_{W の内部双対面} 粘性運動量 flux) / wallArea(W)` の接線成分。= ユーザー選択「全隣接エッジ集約」。

1. **map (新規, ホスト構築)**: 各壁ノードに接続する内部 primal edge 一覧 `wallAdjPlane[]` と `wallIsC0[]`
   (どちらの端点が壁ノードか)。幾何 (`sx/ss/fx/dcc`) は既存 plane 配列を流用。壁ノード→その壁半割面
   (面積・法線) の対応も保持 (`twall` の除算面積・y+ の壁法線距離用)。
2. **`viscousFlux_wall_node_d` (新カーネル)**: `wallAdjPlane` 上でループ。各面で
   - 壁端 `iw`, 内部端 `ii` を `wallIsC0` から決定。
   - **over-relaxed 法線せん断** = `mu*((U[ii]-Uwb)/dcc)*delta`、`Uwb`=`Uxb/Uyb/Uzb`=0 (ghost 不使用)。
     `delta=dcc*ss²/D_safe` は内部カーネルと同式 (mirror-ghost dcc でなく W-I 物理距離)。
   - 面勾配 = **1/2(∇[iw]+∇[ii])** (ユーザー指定 dual 平均) を非直交補正 `·k`・転置項 `·S` に使用。
   - 得た粘性運動量 flux を**壁ノード W ごとに集約** (`atomicAdd` で per-wall-node の twall アキュムレータへ)。
   - **`res_*` (運動量/エネルギー) には加算しない**。
3. **後段で twall/y+ を確定**: 集約した壁ノード粘性力を壁半割面積で割り接線成分→`twall_*_b`。
   `utau=sqrt(|twall|/ro)`、`y+` は**内部ノードまでの壁法線距離** (退化 dcc でなく) で算出し `ypls_b` へ。
4. **node では `viscousFlux_wall_d` の twall/y+ 出力を本カーネルで置換** (運動量残差加算は元々無効なので
   呼ばない or no-op 化)。cell モードは従来どおり `viscousFlux_wall_d` (ビット不変)。
5. **config フラグ** `nodeWallStressEdgeKernel` (既定 1)。0 で旧 (壁ノード勾配ベース twall) に戻し A/B 比較。

### 11.3.1 二重計上の確定切り分け (nNormalPlanes vs nPlanes) と cell/node 非対称性

**`nNormalPlanes == nPlanes` なら境界半割面も内部カーネルが処理し `viscousFlux_wall_d` と二重計上になるが、
node ではそうならない** ([`gmshReader.hpp:1685-1687`](../../solver_density_cuda/mesh/gmshReader.hpp#L1685)):
`nNormalPlanes = nDualInternal` (内部双対面), `nPlanes = nDualInternal + nBHalf` (+境界半割面)。plane 配列は
**内部双対面 `[0,nDualInternal)` + 境界半割面 `[nDualInternal,nPlanes)`** の順。よって:

- `viscousFlux_d` は `ip < nNormalPlanes` = **内部双対面のみ** (境界半割面は触らない)。
- `viscousFlux_wall_d` は `bc.map_bplane_plane_d` = **境界半割面のみ**。
- → **両者の plane 集合は重ならず、この2カーネル間の二重計上は無い** (cell/node 共通)。

**ゆえに二重計上が起きるのは「新カーネルが内部 W-I 面 (`[0,nNormalPlanes)` 内) を再処理して `res_*` に
加算する」場合のみ** → 新カーネルは twall/y+ 専用 (res 不加算) にすれば回避 (§11.3)。

**cell/node の非対称性 (ユーザー指摘「これがややこしい」):**

| | wall 境界面の所属 CV | 内部ノード/セルへの壁せん断 | `viscousFlux_wall_d` の res 加算 |
| --- | --- | --- | --- |
| **cell** | 壁隣接の**内部セル** ic の CV | これ以外に供給源なし (内部カーネルは境界面を skip) | **必須・正しい** (ic に加算; 残す) |
| **node** | **壁ノード W** の CV (W は Dirichlet) | 内部双対面 W-I (`viscousFlux_d`) が供給 | **無効** (W に捨てる); twall/y+ だけ実役 |

→ **cell では `viscousFlux_wall_d` の res 加算が壁せん断の唯一の供給源なので必須**(ユーザー確認の通り)、
**node では res 加算は破棄され実役は twall/y+** という、同一カーネルが離散化で役割を変える非対称が本質。
実装方針: **cell は `viscousFlux_wall_d` をそのまま (res+twall+y+)**、**node はそれを呼ばず
`viscousFlux_wall_node_d` (twall/y+ のみ, res 不変)** に振り分ける。

### 11.4 決定 (ユーザー確認済 2026-06-24) と実装方針

**決定: 案 (A) を採用。** 内部ノード I の運動量せん断は `viscousFlux_d` の内部 W-I 面が担い (既に正しく
加算済・二重計上回避のため新カーネルは res に足さない)、新カーネル `viscousFlux_wall_node_d` は
**`twall`/`y+` の算出専用**。マスク不要・残差/場は不変。(B 置換案は不採用: マスクが §8.2 の地雷・1/2 平均が
`fx` 重みより粗く数値利得が薄いため。)

**実装をローリスク化する判断**: mesh 構築側に「壁ノード→内部ノード」専用マップを**新設しない**。代わりに
既存 `msh.wall_flag_d` (壁ノード=1) で内部面の壁隣接を判定し集約する。これは「全隣接エッジ集約」と
**同一結果**で、§8.1.1/§8.2 が警告する gmshReader/mesh.cpp 改変 (発散多発の地雷原) を回避できる。

実装構成 (node のみ・cell は `viscousFlux_wall_d` を従来どおり呼ぶ):

1. **per-node 粘性力アキュムレータ** `wallStressAccum[nCells*3]` (device, 既存 scratch か新規確保)。
   壁摩擦応力の定義は保存的: 壁ノード W の CV に内部双対面から入る粘性運動量 flux の総和
   = 流体が及ぼす粘性力 = (W は u=0 ゆえ) **壁が流体に及ぼすせん断力**。
2. **kernel-1 `accumWallStress_node_d`** (内部面 `ip<nNormalPlanes` をループ):
   `wall_flag[ic0]^wall_flag[ic1]==1` (片側のみ壁) の面で、壁端を `Uwb=0` として `viscousFlux_d` と
   同形の over-relaxed 粘性運動量 flux を計算し、**壁端ノード W の `wallStressAccum` に `atomicAdd`**
   (面勾配は §11.3 のユーザー指定 1/2 平均)。`res_*` には一切触れない。
3. **kernel-2 `finalizeWallStress_node_d`** (境界半割面 `ib` をループ): `W=bplane_cell[ib]`、
   `twall_*_b[ib] = wallStressAccum[W]/wallArea(W)` の接線成分、`utau=sqrt(|twall|/ro[W])`、
   `y+` は壁法線方向の**内部ノード距離** (退化 dcc でなく) で算出し `ypls_b[ib]` へ。
4. **config フラグ** `nodeWallStressEdgeKernel` (既定 1)。0 で旧 (壁ノード勾配 ∇U[W]·S ベース twall) に
   戻し A/B 比較可能。
5. cell モードは一切不変 (`viscousFlux_wall_d` が res+twall+y+ を従来どおり; §11.3.1 の通り cell は res 必須)。

順序: docs (diffusion theory/impl に「node 壁摩擦応力=内部双対面集約」を追記) → 実装 (上記 1-4) →
検証: case 29/36 node viscous で ①`twall` 平滑化 (TV/range・符号反転が §6.2 比で改善)、②`y+` 妥当
(壁関数が破綻しない)、③運動量場・中心線が**回帰なし** (res 不変ゆえ場は不変のはず) を確認。

### 11.5 実装・検証結果 (2026-06-24, 隔離 worktree)

`solverConfig.{hpp,cpp}` (フラグ `nodeWallStressEdgeKernel` 既定1)、`cuda_forge/viscousFlux_d.cu`
(カーネル `viscousFlux_wall_node_d` + wrapper 配線)、`docs/diffusion/implementation.md` を実装。native build
(RelWithDebInfo→Release, arch86) **クリーン (exit0)**。検証 run (main tree, worktree binary):
`case/36.passive_pseudoshock_control/run_node_wallstress_{off,on,off2}/` (node SST, wallTreatmentSST=1,
nodeWallDirichlet=1, 既存発達場から 30 step)。

1. **無 NaN・完走**: off/on とも exit0、NaN なし、dt ほぼ一致 (4.0936e-9)。
2. **場は不変 (回帰なし)**: OFF vs ON の保存量 maxabsdiff (ro 2.2e-4, roe 22, roOmega 96) は
   **OFF vs OFF2 (同一設定2回=GPU atomicAdd 非決定性ノイズ床: ro 2.5e-4, roe 18, roOmega 64) と同オーダー**
   → 新カーネルは解に影響しない (twall 出力専用) を実証。
3. **twall 改善 (決定的)**: 旧 (`nodeWallStressEdgeKernel=0`) の node twall は **恒等的に 0** (壊れていた)。
   原因: `wallTreatmentSST=1` の modeled τ_w が**壁ノードの接線速度** Ut を使うが、node は ic=壁ノードで
   u=0 (Dirichlet) → Ut=0 → τ≡0 (§11.1 の ic 取り違えそのもの)。新カーネルは内部ノード速度を使い、
   **物理的な壁せん断** (twall_x median 6.0e3 Pa, cf≈0.0046, max 7.2e4 Pa@近shock, TV/range 1.62・
   符号反転2=平滑) を出力。`case/36/compare_node_twall_old_vs_new.png` 参照。
4. **y+ は不変** (新カーネルは ypls に触れない; SST 壁関数入力を変えず場を不変に保つため)。本実装の対象外。

→ **status: done。** 残: cell 側との twall 符号規約整合 (現状 magnitude は物理, 向きは「流体→壁の traction」)、
長時間 run での twall 時系列、案 B (内部ノード運動量も専用化) は不採用維持。

## 変更ログ

- 2026-06-20: 診断カーネル (`viscousWallDiag_d`, env `FORGE_VISC_WALL_DIAG`)、config フラグ
  (`nodeWallViscGradFlux`, `nodeWallDistFloorCoef`)、`viscousFlux_wall_d` の dcc floor 分岐を実装。
  退化を実証 (§6.1) したが、局所距離修正は安定性と非両立 (§6.2-6.3)。既定は元挙動保持。正攻法保留。
- 2026-06-24: §11 追記。ユーザー要望の `viscousFlux_wall_node_d` (壁ノード→内部ノード edge map 上で
  せん断を内部ノードに加算) を設計。**最重要発見: W-I 内部双対面は既に `viscousFlux_d` が処理済み
  → 新カーネルは追加でなく内部カーネルからの抜き取り(置換)が必須 (二重計上回避)**。置換アーキ採用可否を
  ユーザー確認後に docs 更新+実装へ。
- 2026-06-24: 設計を案 A に確定 (ユーザー確認)。内部ノード運動量は `viscousFlux_d` が担い、新カーネルは
  twall/y+ 専用 (res 不変・マスク不要)。さらに「壁ノード残差は破棄・内部ノードは既加算」を確認し、
  実役は twall 出力のみと判明 (§11.2)。`nodeWallStressEdgeKernel` フラグ + `viscousFlux_wall_node_d`
  (壁半割面ごとに wall_flag で入射内部面を集約, Uxb=0・1/2勾配・over-relaxed) を実装。
  **検証 (§11.5): 場は不変 (非決定性ノイズ床内)、旧 node twall は wallTreatmentSST+壁ノード u=0 で恒等0
  だったのを物理的な壁せん断 (cf≈0.0046) に修正。status=done。** 隔離 worktree で実装/ビルド/検証。

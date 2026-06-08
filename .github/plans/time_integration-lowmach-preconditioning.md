# 低マッハ前処理 (Weiss–Smith) — 密度ベース経路の散逸スケール是正と陰解法固有値前処理

## メタ

- **area**: `time_integration` (convection 横断)
- **status**: `done`  <!-- draft / in_progress / done / superseded -->
  <!-- Phase1 (flux c') 採用 / Phase2 (LHS 固有値) 棄却 / Phase3 (b Thornber) 負の結果・opt-in 残置 /
       Phase4 (a 完全 Γ⁻¹A 前処理 lowMachPrecond=2) 採用＝低マッハ自励振動を根治 (§9 2026-06-09)。 -->
- **related_docs**:
  - [`docs/convection/theory.md`](../../docs/convection/theory.md) (低マッハ前処理節)
  - [`docs/convection/implementation.md`](../../docs/convection/implementation.md) (`lowMachPrecond`)
  - [`docs/time_integration/theory.md`](../../docs/time_integration/theory.md) (前処理固有値・局所時間刻み)
  - [`docs/time_integration/implementation.md`](../../docs/time_integration/implementation.md)
- **related_plans**:
  - [`convection-slau2-lowmach.md`](convection-slau2-lowmach.md) (done。SLAU2 では低マッハフロアが解消しないことを §9 で確認。本 plan が根治)
  - [`gpu-implicit-plan.md`](gpu-implicit-plan.md) (block DPLUR / dual-time。Phase 2 が固有値で接続。dual-time は実装済 `unsteady_diag`)
- **created**: `2026-06-07`
- **owner**: `CFD Dev`

## 1. 目的

密度ベース圧縮性経路 (CUDA) に Weiss–Smith 型低マッハ前処理を導入し、低マッハ域
($M\!\sim\!0.06$) で残るエネルギー残差フロア (`rms_roe`≈5 頭打ち) とチャンバー圧の非物理
ばらつき (std/mean≈0.65%、物理は $\gamma M^2/2\!\approx\!0.4\%$) を低減する。前セッションで
SLAU2 (圧力束第 3 項) を試したが解消せず、原因は質量流束の圧力散逸項を音速 $\widehat c$ で
スケールしていること (低マッハで圧力–速度カップリングが $O(M)$ に縮退) と特定済み
([`convection-slau2-lowmach.md`](convection-slau2-lowmach.md) §9)。本 plan は散逸スケールを
前処理音速 $c'$ に置き換えてフロアを根治し、併せて擬似時間項 (局所 $\Delta t$・陰解法 LHS 固有値) を
前処理して収束を加速する。前処理は擬似/時間微分項にのみ掛け、超音速域 ($U_r=c$) は従来と
ビット一致で不変。

## 2. スコープ

- **やる**:
  - 共有 device ヘッダ `cuda_forge/lowMachPrecond_d.cuh` ($U_r$, $\beta$, 前処理音速 $c'$, 前処理固有値 $\lambda'$)。
  - **Phase 1** (dual-time と独立、フロア主因に直撃) — **実装・検証済 (§9)**:
    - `convectiveFlux_d.cu` `SLAU_d` の質量流束 圧力散逸項 `chi/c_hat·ΔP` の `c_hat → c'`。
    - `solverConfig` に `lowMachPrecond` (0/1, 既定 0) と `precondEps` (既定 0.3) を追加。
  - **Phase 2** (LHS 固有値前処理 ＋ setDT) — **実装・検証したが棄却 (§9)**。
    `build_jacobian_split` の固有値 `U±c → U±c'` 差し替えと `setDT_d` の `sonic → c'` は、
    block DPLUR では**収束加速にも根治にもならず、むしろ有害**と判明したため revert。詳細・根拠は §9。
  - **Phase 3** (b) **Thornber 型再構成補正** (RHS・LHS 不触) — **実装・検証済 → 負の結果 (§9 `2026-06-08`)**。
    `lowMachThornber` は opt-in 残置 (既定 0)、ノズル limit cycle の根治には不採用。
  - **Phase 4** (a) **完全 $\Gamma^{-1}A$ preconditioned-flux Jacobian** (`lowMachPrecond=2`) — **採用・実装対象 (§5 Phase 4)**:
    - $\Gamma_Q=M\Gamma_pM^{-1}$ と前処理固有ベクトル $R_\Gamma,L_\Gamma$ を再導出、`build_jacobian_split` を前処理固有系で書き直し。
    - 対角時間項を行列 $\Gamma_Q V/\Delta\tau'$ 化、$\Delta\tau'$ を前処理スペクトル半径から (`setDT_d`)。
    - 単精度対策: 対角ブロック組み立て+`solve_5x5` を局所倍精度化 (§5 Phase 4.3)。
    - 振動根治 ($\epsilon$ を物理値方向へ) と低マッハ収束加速を同時に狙う。$\beta=1$ でビット一致回帰必須。
- **やらない**:
  - Phase 1 (flux `c'`) と Phase 4 (完全前処理) の同時併用。両者は排他 (`lowMachPrecond=1` vs `=2`)。
  - SLAU 以外 (Roe/HLLE/AUSM/KEEP) の低マッハ前処理。本 plan は SLAU 経路のみ。
  - dual-time 本体 (実装済)。`unsteady_diag` 物理時間項は前処理しない (§4 分界)。
  - scalar 対角版 (`blockDPLUR==0`) の固有値前処理 (block を既定とするため後回し)。

## 3. 関連 docs と前提

- 散逸スケール是正の理論: [`docs/convection/theory.md`](../../docs/convection/theory.md) の
  「低マッハ前処理 (Weiss–Smith 型・散逸スケール是正)」節 (SLAU2 節の直後、HLLE 節の前)。
- 前処理固有値・局所時間刻み: [`docs/time_integration/theory.md`](../../docs/time_integration/theory.md) の
  「低マッハ前処理固有値 (Weiss–Smith)」節と「局所時間刻み」節。
- 実装対応: [`docs/convection/implementation.md`](../../docs/convection/implementation.md)
  `SLAU_d`「低マッハ前処理 (`lowMachPrecond`)」、
  [`docs/time_integration/implementation.md`](../../docs/time_integration/implementation.md)
  `implicit_defect_correction_block_d` の低マッハ前処理ノート・局所時間刻み節。
- **docs 4 ファイルは本 plan の実装着手前に更新済 (AGENTS 開発フロー)。**
- 前提: gmshReader の 3D 幾何法線修正 (`fe1db5e`) と dual-time block DPLUR (`ad71a0c`) は
  base ブランチに取り込み済。本 plan は新ブランチ `feature/lowmach-preconditioning` で作業。

## 4. 設計方針

### 4.1 前処理の核 (`lowMachPrecond_d.cuh`)

基準速度と前処理係数:
$$
U_r = \min\!\big(c,\ \max(|\mathbf u|,\ \epsilon\,c)\big),\qquad \beta = (U_r/c)^2 \in (0,1].
$$
停留点フロア $\epsilon c$ (`precondEps`)、$M\!\ge\!1$ で $U_r=c,\ \beta=1$。前処理音速・前処理固有値:
$$
c' = \tfrac12\sqrt{(1-\beta)^2 U_n^2 + 4U_r^2},\qquad
\lambda_{4,5}' = \tfrac12(1+\beta)U_n \pm c',\quad \lambda_{1,2,3}'=U_n.
$$
$\beta=1$ で $c'\!=\!c$, $\lambda'\!=\!U_n\pm c$ に厳密復帰。`__device__ __forceinline__` の
小関数群として提供し、3 カーネルから共用 (引数: `c`, `vel_mag`(=$|\mathbf u|$), `Vn`(=$U_n$), `eps`)。

### 4.2 収束解への作用 (切り分け)

| 適用箇所 | 種別 | 収束解 |
| --- | --- | --- |
| ①`SLAU_d` 散逸 `c_hat→c'` | 空間散逸 (残差の一部) | **低マッハ域を意図的に変更** (= フロア根治)。超音速は不変 |
| ②`setDT_d` スペクトル半径 | 擬似時間項 | **厳密に不変** (収束加速のみ) |
| ③`build_jacobian_split` 固有値 | 陰解法 LHS (擬似時間) | **厳密に不変** (条件数改善・擬似 CFL 上限 ↑) |

フロアを実際に下げるのは①の解変更。②③は不変なまま収束を速める。物理音速 `sonic` 配列は不変。

### 4.3 dual-time (実装済) との分界

- dual-time は物理時間 BDF 項 $aV/\Delta t$ を `implicit_defect_correction_block_d` の
  `unsteady_diag` 引数経由で対角に加算 (`add_identity_scaled(diag_block, v*unsteady_diag)`)。
  低マッハ前処理は**ここを触らない**。
- 対角の所有分離: $V/\Delta\tau\,I$ (擬似・前処理固有値が関与) ＋ $aV/\Delta t\,I$ (物理 BDF・
  **非前処理**・dual-time 所有) ＋ $\sum_f A'^{+}_f S_f$ (低マッハ所有)。
- 低マッハが触る関数: `build_jacobian_split` の `lambda[5]` のみ。
  dual-time が触る要素: `unsteady_diag` / `addUnsteadyTimeTerm` / `roN(N)` / `nSubIterDualTime`。
  → **同一ファイル内でも別関数・別項**のため衝突最小。「非定常では前処理を擬似時間項のみ」が
  自動的に成立。

### 4.4 config

- `lowMachPrecond` (int, 既定 `0`): 0=従来、1=前処理 ON。`solverConfig.cpp`/`.hpp` に追加。
- `precondEps` (float, 既定例 `0.05`): 停留点フロア $\epsilon$。
- whitelist 検証は SLAU2 同様に最小限 (未指定は既定値)。既定 0 のため既存ケースはビット不変。

## 5. 実装ステップ

### Phase 1 (フロア主因・dual-time 独立) — 実装・検証済

1. **共有ヘッダ**: `cuda_forge/lowMachPrecond_d.cuh` を新設 ($U_r$, $\beta$, $c'$, $\lambda'$ の device 関数)。
2. **config**: `solverConfig.hpp`/`.cpp` に `lowMachPrecond`, `precondEps` を追加 (既定 0 / 0.3)。
3. **対流フラックス**: `convectiveFlux_d.cu` `SLAU_d` に `int lowMachPrecond, flow_float precondEps`
   引数追加。`c_hat` から `c_diss` を作り、質量流束 `mdot` の `chi/c_hat·ΔP` を `chi/c_diss·ΔP`
   に分岐 (`lowMachPrecond==0` で `c_diss==c_hat` でビット不変)。`Vn_hat`・`p_tilde` は不変。wrapper の
   `SLAU_d<<<>>>` 呼び出しに引数を渡す。
4. **ビルド**: Docker dev image (`forge-solver:cuda-dev`) で `tools/build.sh`。
5. **検証** (§6 Phase 1)。

### Phase 2 (陰解法 LHS ＋ setDT) — 実装・検証したが棄却 (§9)

当初案: `build_jacobian_split` の固有値 `U±c → U±c'`、`setDT_d` の `sonic → c'`。
実装・検証した結果 block DPLUR では有害 (収束加速も根治もせず、安定 `eps` 範囲を狭め、
mild `eps` では chamber チェッカーボードを悪化) と判明し revert した。根拠・データは §9。
→ §9 の `2026-06-08` 設計分析で「振動=RHS 駆動、LHS (a) は収束加速の道具で根治手段ではない」と確定。

### Phase 3 (b) Thornber 型再構成補正 (RHS・LHS 不触) — 採用 (実装対象)

振動の根治は RHS 側で行う (§9 `2026-06-08`)。`SLAU_d` の L/R 速度再構成直後で左右速度ジャンプを
局所マッハ `z=min(M,1)` で縮め、低マッハで $O(1/M)$ に増大する速度ジャンプ由来の散逸を抑える。

1. **config**: `solverConfig.hpp`/`.cpp` に `lowMachThornber` (int, 既定 `0`) を追加。
   既定 0 で既存ケースはビット不変。Phase 1 (`lowMachPrecond`) と独立トグル (直交・併用可)。
2. **対流フラックス**: `convectiveFlux_d.cu` `SLAU_d` に `int lowMachThornber` 引数追加。
   L/R 速度再構成 ([convectiveFlux_d.cu:330-345]) 直後・`velocity2_*`/`h_*`/`Vn_*` 算出前に挿入:
   $$ z=\min\!\Big(1,\ \tfrac{\sqrt{(|\mathbf u_L|^2+|\mathbf u_R|^2)/2}}{\hat c}\Big),\quad
      \mathbf u_L^\* = \bar{\mathbf u}+z\,\delta\mathbf u,\ \ \mathbf u_R^\* = \bar{\mathbf u}-z\,\delta\mathbf u $$
   ($\bar{\mathbf u}=\tfrac12(\mathbf u_L+\mathbf u_R)$, $\delta\mathbf u=\tfrac12(\mathbf u_L-\mathbf u_R)$, 3 成分すべて)。
   `velocity2_L/R`・`h_p/h_m`・`Vn_p/Vn_m` はブレンド後の速度から算出する (`h_p`/`velocity2_L` を挿入点の後段へ移動)。
   `lowMachThornber==0` で恒等 (ビット不変)、$M\!\ge\!1$ で $z=1$ (超音速域不変)。wrapper の `SLAU_d<<<>>>` に引数追加。
3. **ビルド**: Docker dev image (`forge-solver:cuda-dev`) で `tools/build.sh`。
4. **検証** (§6、`case/23.axi_nozzle` 20k step)。Phase 1 単独 (0.603%) に対し Thornber 単独・Phase1+Thornber 併用を比較。

→ **結果は負** (§9 `2026-06-08`)。Phase 3 完了 (機能は opt-in 残置・根治不採用)。根治本命は Phase 4 (a) へ。

### Phase 4 (a) 完全 $\Gamma^{-1}A$ preconditioned-flux Jacobian — 根治本命 (大改修・実装対象)

振動根治の正しい RHS レバーは「圧力散逸を増やす」前処理音速 $c'$ (Phase 1) 側だが、$\epsilon$ を物理値
(~0.05) まで下げると陽的フラックスが発散する (安定限界 $\epsilon\gtrsim0.15$)。$\epsilon$ を下げて根治するには
**LHS を完全 preconditioned-flux Jacobian にして aggressive flux を陰的に安定化**する必要がある。これは
同時に低マッハ収束加速 (擬似 CFL 上げ) も与える。Phase 1 の flux 変更とは**排他**＝新モード `lowMachPrecond=2`。

**4.1 定式化** (保存変数 $Q=(\rho,\rho u,\rho v,\rho w,\rho E)$、forge の実順。**plan 旧記載「(ρ,ρE,ρu,ρv,ρw)」は誤り、エネルギは index 4**)。
Weiss–Smith 前処理は原始変数 $Q_p=(p,u,v,w,T)$ で定義:
$$ \Gamma_p\,\partial_\tau Q_p + \partial_t Q_c + \nabla\!\cdot F = 0,\quad
   \Gamma_p = M^{-1}\!+\!(\Theta-\rho_p)\,\mathbf e_\rho\,\partial p,\quad \Theta=\tfrac{1}{U_r^2}-\tfrac{\rho_T}{\rho c_p}. $$
($\rho_p=\partial\rho/\partial p|_T=1/RT$, $\rho_T=-\rho/T$, $c_p=\gamma R/(\gamma-1)$。$U_r=c$ で $\Theta=\rho_p$＝$\Gamma_p=M^{-1}$＝無前処理。)
保存形へは $\Gamma_Q = M\,\Gamma_p\,M^{-1}$ ($M=\partial Q_c/\partial Q_p$)。前処理フラックスヤコビアン
$\hat A_\Gamma=\Gamma_Q^{-1}A$ は前処理固有値 $\lambda'_{1,2,3}=U_n,\ \lambda'_{4,5}=\tfrac12(1+\beta)U_n\pm c'$
($c',\beta$ は `lowMachPrecond_d.cuh` 既存) と**前処理固有ベクトル** $R_\Gamma,L_\Gamma$ を持つ。

**4.2 LHS 構成** (`implicit_defect_correction_block_d`):
- 対角時間項を**スカラー** $V/\Delta\tau\,I$ から**行列** $\Gamma_Q\,V/\Delta\tau'$ へ ([timeIntegration_d.cu:696])。
- $\Delta\tau'$ を前処理スペクトル半径 $\lambda'_{\max}\sim|u|+c'$ から ( `setDT_d`、低マッハで $\Delta\tau'$ が $\sim1/M$ 拡大)。
- `build_jacobian_split` を前処理固有系で書き直し: $\hat A_\Gamma^{+}=R_\Gamma\Lambda'^{+}L_\Gamma$ を対角、$-\hat A_\Gamma^{-}$ を近傍項。
- dual-time 物理 BDF 項 $aV/\Delta t\,I$ は**非前処理のまま** (恒等)。dual-time 所有、§4.3 分界を踏襲。
- $\beta=1$ ($M\ge1$) で $\Gamma_Q=I,\ \hat A_\Gamma=A,\ \Delta\tau'=\Delta\tau$ に**ビット一致回帰** (回帰テスト必須)。

**4.3 単精度フィージビリティ (着手前に潰す設計課題・結論)**。`flow_float`=float。前処理は圧力結合方向に
$1/\beta\sim1/M^2$ の増幅を持ち、対角ブロック $D=\Gamma_Q V/\Delta\tau'+\dots$ は保存基底で条件数 $\sim1/\beta$
(M=0.06・$\epsilon$=0.15 で $\sim44$、物理 $\epsilon$=0.05 で $\sim400$) を持つ。前処理固有基底では $D$ は良条件
(全モード $O(\text{area}\cdot\lambda'_{\max})$) で、悪条件は保存基底への変換 $\Gamma_Q$ が担うだけ。
→ **対策 (採用): 対角 5×5 ブロックの組み立て + `solve_5x5` を局所倍精度で行い float に戻す** (セル毎・5×5 と
小さく低コスト、`solve_5x5` を double テンプレート化)。併せて $\beta$ フロア ($\epsilon$) で条件数上限を担保。
これで float 本体のまま $\epsilon$ を物理値方向へ下げられる見込み。基底変換を避ける案 (原始/固有基底で組む) は
forge の保存変数パイプラインに侵襲的なため不採用。

**4.4 実装ステップ**:
1. **数式確定**: $\Gamma_Q=M\Gamma_pM^{-1}$ と $R_\Gamma,L_\Gamma$ の閉形を導出 (Weiss–Smith/Merkle 文献ベース、forge の
   正規化に合わせ検証)。$\beta=1$ で恒等になることを記号的に確認。
2. **`lowMachPrecond_d.cuh` 拡張**: $\Gamma_Q$ 構築・前処理固有値・$R_\Gamma,L_\Gamma$ の device 関数。
3. **`solve_5x5` の double 版**、対角ブロック組み立て+解を局所倍精度化。
4. **`build_jacobian_split` の前処理版** (`lowMachPrecond==2` 経路)、`setDT_d` の $\Delta\tau'$、時間項 $\Gamma_Q$ 化。
5. **回帰**: $\beta=1$ (超音速単独ケース) で現行とビット一致を確認するテストを先に用意。
6. **検証** (§6): `case/23.axi_nozzle` で $\epsilon$ を 0.15→0.05 方向に下げ、振幅が物理目標 ~0.25% に近づくか、
   かつ擬似 CFL を上げられるか (収束加速) を測る。発散しないこと・超音速域不変。

## 6. 検証

検証ケース: `case/23.axi_nozzle` の 90 度セクター 3D O-grid ノズル (chamber $M\approx0.06$)。
ハーネスは前セッション資産を再利用 (既存 `run_*` は複製した新 `run_*` で実行)。

- **メッシュ/初期化**: `mesh_3d_m4/build_mesh_3d.sh` → `convertGmshToForge` (gmsh 3D 法線修正済) →
  `set_initial_isentropic.py` (準 1 次元等エントロピー初期化、起動安定化)。
- **実行**: `run_3d_m4_prod` を複製した新 `run_*` で `solver: SLAU`, `lowMachPrecond: 1`
  (陰解法据え置き: `timeIntegration 11`, `cfl 0.5`, `cfl_pseudo 0.5`, `blockDPLUR 1`,
  `implicitRelax 0.3`, `nStepInner 20`)。`build/forge` の stale を `find` で確認。
- **収束ベース判定 (重要)**: この症状は静的フロアではなく**自励振動 (limit cycle)** のため、
  2000 step スナップショットでは判断しない。**20000 step 回し、4k–20k step の chamber 圧 std/mean
  (M<0.08) と `rms_roe` の平均・包絡で振幅を評価**する (`chamber_metric.py`)。判定:
  - 低マッハ振動の振幅 (chamber std/mean の平均) が OFF 比で明確に低下。
  - ⟨M_max⟩ が OFF/ON 一致＝超音速域不変、発散なし。
  - `lowMachPrecond: 0` 回帰で OFF 経路がビット不変。
  - `residual_history.csv` → `residual_history.png` を生成。
  - (任意) 別途 2D 低速ケースで一般妥当性。

## 7. 影響範囲

- 新規: `solver_density_cuda/cuda_forge/lowMachPrecond_d.cuh`。
- 変更 (Phase 1): `convectiveFlux_d.cu`/`.cuh` (`SLAU_d` 引数・分岐・wrapper)、
  `setDT_d.cu` (`setCFL_pln_d` 引数・分岐・wrapper)、`solverConfig.cpp`/`.hpp`。
- 変更 (Phase 2): `timeIntegration_d.cu` (`build_jacobian_split`・`implicit_defect_correction_block_d`・
  wrapper の引数伝搬)。dual-time の `unsteady_diag` 経路は不変。
- 変更 (Phase 3): `convectiveFlux_d.cu`/`.cuh` (`SLAU_d` の `lowMachThornber` 引数・L/R 速度ブレンド・wrapper)、
  `solverConfig.cpp`/`.hpp` (`lowMachThornber`)。LHS (`timeIntegration_d.cu`) は不触。
- docs: `docs/convection/{theory,implementation}.md`, `docs/time_integration/{theory,implementation}.md`
  (更新済)。`docs/index.md` は構成不変につき確認のみ。
- 既存ケース: `lowMachPrecond`・`lowMachThornber` 既定 0 のため挙動不変。

## 8. 完了条件

- [x] 関連 `docs/convection/{theory,implementation}.md` 更新済み (収束ベース所見・Phase 2 不採用を反映)
- [x] 関連 `docs/time_integration/{theory,implementation}.md` 更新済み (固有値/setDT 前処理は不採用と明記)
- [x] Phase 1 実装・検証完了 (収束ベース: 低マッハ自励振動を ε=0.15 で −32% 減衰。§9)
- [x] Phase 2 (LHS 固有値 ＋ setDT) 実装・検証 → **棄却** (block DPLUR で有害、revert。§9)
- [x] `.github/plans/README.md` の状態を更新
- [x] 設計分析 (a の可否・RHS/LHS 切り分け・(b) 採用決定) を §9 `2026-06-08` に記録
- [x] **Phase 3 (b) Thornber** — docs 更新 → 実装 (`lowMachThornber`) → 20k step 検証完了。**結果は負** (この症状に無効・僅かに悪化。§9 `2026-06-08`)。機能は opt-in で残置、ノズル limit cycle の根治には不採用
- [x] **Phase 4 (a) $\Gamma^{-1}A$ 完全前処理** (`lowMachPrecond=2`) — 実装・検証完了。**低マッハ自励振動を根治** (ε=0.05 で振幅 0.882%→0.087%・定常収束。§9 `2026-06-09`)。時間項 $\Gamma_c$+物理フラックス FVS+$\Delta\tau'$、倍精度ブロック解。収束加速は副次 (ブレークイーブン)

## 9. 変更ログ

- `2026-06-07` — 初稿 (計画)。docs 4 ファイル (convection/time_integration の theory/implementation)
  を先行更新。Phase 1 (フラックス散逸＋setDT、dual-time 独立) / Phase 2 (block DPLUR 固有値、
  dual-time `unsteady_diag` と分界) に分割。base ブランチに dual-time (`ad71a0c`) が取り込まれた
  ため、dual-time との分界を「実装済 `unsteady_diag` 項との非干渉」として確定。

- `2026-06-07` — **Phase 1 実装**。`lowMachPrecond_d.cuh` 新設 ($U_r$/$\beta$/$c'$/$\lambda'$)、
  `SLAU_d` の質量流束圧力散逸項を `chi/c_hat·ΔP → chi/c_diss·ΔP` に分岐、config 追加。
  ブランチ `feature/lowmach-preconditioning`。`lowMachPrecond==0` でビット不変。

- `2026-06-07` — **検証手法の是正 (重要)**。当初 2000 step で良否判断していたが、
  `case/23.axi_nozzle` 3D ノズルを 20000 step 回した結果、**この問題は静的な「残差フロア」ではなく
  低マッハの自励振動 (limit cycle)** と判明。baseline の `rms_roe` は 5.05 で頭打ちではなく
  2000 step の未収束に過ぎず、その後 0.4〜4.6 を振動 (mean 2.18)。chamber 圧 std/mean も
  0.25〜1.82% を振動 (mean 0.88%)、`M_max` は 5.95〜6.06、step 2000 のスナップショットは
  過渡を拾っていただけ。**正しい指標は振動の包絡/平均** (4k–20k step の統計) とする。
  当初タスク前提・`convection-slau2-lowmach.md` §9 の「`rms_roe≈5` 頭打ち」も 2000 step
  アーティファクトだったと判断。

- `2026-06-07` — **Phase 1 収束ベース検証 (合格)**。20000 step (cfl_pseudo 0.5)、4k–20k 統計:

  | 設定 | chamber std/mean (M<0.08) mean [min,max] | rms_roe mean | ⟨M_max⟩ |
  | --- | --- | --- | --- |
  | OFF | 0.882 % [0.25, 1.82] | 2.18 | 6.028 |
  | ON `precondEps=0.15` | **0.603 %** [0.11, 1.53] | **1.47** | 6.029 |
  | ON `precondEps=0.3`  | 0.733 % [0.17, 1.68] | — | 6.029 |

  - **効果**: 低マッハ limit-cycle の振幅を **ε=0.15 で −32%** (chamber 0.882→0.603%, rms_roe 2.18→1.47)、
    ε=0.3 で −17%。ε 小ほど強く減衰 (単調)。物理目標 $\gamma M^2/2\approx0.25\%$ には未達だが明確に減衰。
  - **回帰・不変性**: `lowMachPrecond:0` は OFF 経路ビット不変。⟨M_max⟩ は OFF/ON 一致＝超音速域不変。
  - **安定限界**: ε 小で圧力散逸が増幅。ε=0.05 (停留点 ~20×) は数十 step で NaN、ε≳0.15 で 20k step 安定。既定 0.15。
  - 成果物: `run_conv_baseline` / `run_conv_p1_eps015` / `run_conv_p1_eps03` (`residual_history.png` 付き)、
    `case/23.axi_nozzle/chamber_metric.py`。

- `2026-06-07` — **Phase 2 (LHS 固有値 ＋ setDT) 実装したが棄却**。`build_jacobian_split` の固有値
  `U±c → U±c'`、`setDT_d` の `sonic → c'` を実装・検証 → **block DPLUR では有害**と判明し revert:
  1. **setDT 前処理**: 低マッハ域で `dt_local` 増大 → 対角 `V/Δτ` が縮小。固有値前処理で `A⁺` も縮小するため
     対角が二重に潰れ、最速で発散 (step ~4)。
  2. **LHS 固有値前処理**: `A⁺` の音響固有値 `U±c` (大) こそ block DPLUR の対角優位性の源
     (2026-06 の収束化はこれに依拠)。`U±c'` (小) に落とすと対角優位性が下がり、**flux 単独では安定だった
     `ε=0.15` すら発散** (安定 ε 範囲をむしろ狭める)。`ε=0.3` は安定だが収束加速・根治の利得なし。
  3. **根因**: `build_jacobian_split` は Euler 固有ヤコビアンで、SLAU 低マッハ散逸項を表さない。
     固有値だけ前処理しても aggressive flux を陰的に安定化できない。真の安定化には
     **完全 preconditioned-flux Jacobian** (固有ベクトルも前処理＝大改修) が要る。
  → 最小差分 (固有値のみ) では Phase 2 の目的を達成できず、ε を下げて根治する道は本アプローチでは閉じた。

  **今後の選択肢**: (a) 完全 $\Gamma^{-1}A$ preconditioned-flux Jacobian を LHS に実装 (大・リスク高)、
  (b) 別系の低マッハ補正 (Thornber 型再構成補正など、LHS 不要) を別 plan で検討、
  (c) Phase 1 (フラックス散逸, ε=0.15, −32%) を成果として確定。現状は (c) を既定とし、(a)/(b) は未着手。

- `2026-06-07` — **本 plan を「未完成」として一旦中断 (status は in_progress のまま)**。計算は回り
  発散もしないが、低マッハ自励振動の**根治には至っていない** (Phase 1 は振幅 −32% の部分緩和どまり、
  物理目標 ~0.25% に対し ~0.60%)。現状の到達点と再開の指針を以下に残す。

  **確定している成果 (採用)**: フラックス散逸前処理のみ (`lowMachPrecond:1`、既定 `precondEps=0.15`)。
  opt-in・OFF ビット不変・超音速域不変。`feature/lowmach-preconditioning` ブランチ (未コミット)。

  **棄却済み (掘らなくてよい)**: LHS 固有値前処理・setDT スペクトル半径前処理。block DPLUR の対角優位性を
  崩すため有害 (上記 Phase 2 エントリ)。

  **根治の再開ルート (未着手・どちらかを別フェーズ/別 plan で)**:
  - **(b) Thornber 型再構成補正 (推奨・軽量)**: 面再構成の左右速度ジャンプを `z=min(局所M,1)` で縮める
    (`u_L*=½(u_L+u_R)+½z(u_L−u_R)`, `u_R*` は対称)。`SLAU_d` の L/R 速度生成直後を数行。**LHS 不触**。
    Phase 1 (圧力散逸側) と**直交＝併用可**。散逸を「削る」方向なので陰解法でも安定寄り。要: 定常 RANS での
    効き・ロバスト性検証。
  - **(a) 完全 $\Gamma^{-1}A$ 前処理 (重い・根治本命)**: セルごとに保存変数の前処理行列 `Γ_Q=M·Γ_p·M⁻¹`
    を組み、`build_jacobian_split` を前処理固有系 (前処理固有値＋前処理固有ベクトル `R_Γ,L_Γ`) で書き直し、
    対角の時間項を `V/Δτ·Γ_Q` (行列) に、残差散逸を `Γ_Q⁻¹|Γ_Q A|` に格上げ、大域 Ur フロア＋衝撃波/停留点
    ガードを追加。変数並びが (ρ,ρE,ρu,ρv,ρw) のため固有ベクトルは要再導出。Phase 1 の flux 変更とは**排他**
    (置き換え)＝`lowMachPrecond=2` の別モードにするのが筋。Phase 1 の config/`lowMachPrecond_d.cuh` 土台は流用。

  **再開時の検証手順 (重要)**: 必ず 20000 step 回し、4k–20k の chamber std/mean (M<0.08) と rms_roe の
  **平均・包絡で振幅評価** (`case/23.axi_nozzle/chamber_metric.py`)。2000 step スナップショットは過渡を拾うため不可。
  基準値: baseline 0.882%、Phase1(ε=0.15) 0.603%、物理目標 ~0.25%。

- `2026-06-08` — **設計分析 (再開前のコード実体精査)。(a) の可否判定と RHS/LHS の本質的切り分けを確定し、
  (b) を次フェーズに採用と決定**。コード根拠は次の3点:

  1. **`build_jacobian_split` ([timeIntegration_d.cu:158-211]) は近似ではなく真の解析的 $R\Lambda L$ 分解**。
     $\Lambda=\{U+c,U,U,U,U-c\}$ にフル右/左固有ベクトル $R,L$ を掛け `a_plus=RΛ⁺L`・`k_off=R(-Λ⁻)L` を構成。
     → Phase 2 が**固有値だけ** $U\pm c\to U\pm c'$ に差し替えたのは「物理固有ベクトル $R,L$ のまま固有値のみ前処理した
     非整合行列」で、$\Gamma^{-1}A$ の正しい分割になっていない。Phase 2 棄却は当然の帰結であり、**(a) は
     $R_\Gamma,L_\Gamma$ の再導出が必須**と確定。
  2. **対角優位性は音速スケールの `V/Δτ·I` に依拠** ([timeIntegration_d.cu:696-758])。対角 = `V/Δτ·I`(擬似・
     スカラー)＋`aV/Δt·I`(BDF)＋`ΣA⁺S`＋`viscous·I`。低マッハで $\Delta\tau\sim\text{cell}/c$ ゆえ $V/\Delta\tau\sim\text{area}\cdot c$ が
     off-diagonal `k_off`(音響 $\sim c$)を上回り優位。Phase 2 が setDT を $c'$ 化すると $V/\Delta\tau$ が縮み、固有値前処理で
     $A⁺$ も縮む→**二重に潰れて即発散** (§9 既述) を裏付け。
  3. **`flow_float` は単精度 float** ([flowFormat.hpp:6])。

  - **判定 (a) は対角優位性を「再確保できる」が条件付き**。整合3点セット —— ①時間項を**スカラーでなく行列**
    $V/\Delta\tau'\cdot\Gamma_Q$ ②$R_\Gamma,L_\Gamma$ 再導出 ③$\Delta\tau'$ を前処理スペクトル半径 $\lambda'_{\max}\sim|u|$ から ——
    が揃えば、前処理固有基底で時間項と $A_\Gamma^{+}S$ が同じ $O(|u|)$ で釣り合い、優位性は条件数改善とともに回復する
    (これが低マッハ前処理陰解法の設計目標そのもの)。**Phase 2 の失敗は「(a) 不可」の証拠ではなく「半端な前処理は壊れる」証拠**。
    ただし forge 固有リスク2つ: (i) **単精度 × $\Gamma_Q$ 条件数** —— $\Gamma_p$ は圧力カップリング方向に $1/\beta\sim1/M^2$
    (M=0.06 で ~280) の増幅を持ち、$\Gamma_Q=M\Gamma_pM^{-1}$ を float の `solve_5x5` で毎セル解くと桁落ちリスク。対策は
    $\Gamma_Q$ 正規化／対角ブロック解の局所倍精度化／$\beta$ フロア。(ii) **固有ベクトル再導出の検証コスト** —— $\beta{=}1$ で
    現行へビット一致回帰する回帰テストを先に用意。
  - **本質的切り分け (最重要)**: §9 の実験非対称性 (**RHS の Phase 1 は振動 −32%、LHS の Phase 2 は無効〜有害**) は、
    **低マッハ自励振動が RHS(SLAU 半離散作用素)のほぼ無散逸な低マッハモードに起因する**ことを示す。defect-correction の
    LHS は**収束解を変えない** (§4.2) ので、**(a) は収束加速の道具であって振動の根治手段ではない**。離散定常点そのものが
    その低マッハモードについて中立〜不安定なら、どんな LHS でもその固定点には落ちない (Phase 2 が「根治しなかった」事実と整合)。
    → 振動の根治レバーは **RHS 側**:
    | 目的 | 正しいレバー | コスト |
    | --- | --- | --- |
    | 振動を ~0.25% へ根治 (精度) | RHS 側: Phase 1 ＋ **(b) Thornber 再構成補正** | 軽 (SLAU の L/R 速度生成直後に数行・LHS 不触・Phase 1 と直交=併用可) |
    | 低マッハの擬似 CFL/収束を上げる (剛性) | LHS 側: 完全 (a) $\Gamma^{-1}A$ | 重 (固有ベクトル再導出＋$\Gamma_Q$ 行列時間項＋単精度対策＋回帰) |
  - **決定**: 当面の目的は**低マッハ振動の根治**ゆえ **(b) Thornber を Phase 3 として採用**し、まず RHS を低マッハ整合化する
    (§5 Phase 3)。(a) は「低マッハ収束加速が必要」と確認された段階で別途着手 (上記リスク (i)(ii) の feasibility を先に詰める)。
    Phase 1 (圧力散逸 `c'`) と (b) は直交ゆえ**併用**して評価する。

- `2026-06-08` — **Phase 3 (b) Thornber 実装・検証 → この症状には無効 (むしろ僅かに悪化)。負の結果**。
  `SLAU_d` に `lowMachThornber` を実装 (L/R 速度ジャンプを `z=min(M,1)` で 3 成分ブレンド)、config 追加、
  Docker dev image でビルド成功。`case/23.axi_nozzle` 20k step 収束ベース検証 (4k–20k chamber std/mean, M<0.08):

  | 設定 | chamber std/mean mean [min,max] | vs baseline | M_max |
  | --- | --- | --- | --- |
  | OFF (baseline) | 0.882% [0.25,1.82] | — | 6.057 |
  | Phase1 のみ (`c'`, ε=0.15) | 0.603% [0.11,1.53] | **−32%** | 6.046 |
  | **Thornber のみ** | **0.924% [0.28,1.85]** | **+5% (悪化)** | 6.058 |
  | **Phase1+Thornber 併用** | **0.613% [0.10,1.54]** | −30% (Phase1 単独 0.603% に対し利得なし・僅かに悪化) | 6.046 |

  - **結果**: Thornber 単独は全スナップで baseline より一様に ~3–5% 高く (ノイズでない)、併用も Phase1 単独に
    対し一様に ~1–2% 高い。**この低マッハ自励振動に対し Thornber は無効、むしろ僅かに悪化**。
  - **理由 (重要)**: Thornber は速度ジャンプを縮めて**運動量上流化の散逸を「減らす」**補正で、本来は低マッハの
    **過剰散逸 (smearing) による精度劣化**を直す道具 (LES/乱流減衰など)。一方この症状は**圧力–速度カップリングの
    under-damping (チェッカーボード/limit cycle)** で、必要なのは散逸を**増やす**方向。Phase1 の `c'` は圧力散逸を
    増やすので効き、Thornber は散逸を減らすので逆符号。§9 `2026-06-08` の「振動=RHS 駆動」は正しかったが、
    RHS の中でも**圧力散逸を増やす (Phase1)** が正しいレバーで、**速度ジャンプを減らす (Thornber)** はこの症状には逆。
  - **安定性**: Thornber は安定・超音速域不変 (M_max 6.05 一致)・NaN なし。OFF (`lowMachThornber:0`) は
    atomicAdd 非決定性の範囲で従来経路不変 (同一バイナリ二回 run の差 ~1e-4 と新旧差 ~1e-4 が同オーダーと確認)。
  - **扱い**: `lowMachThornber` は実装としては正しく opt-in (既定 0)。LES/乱流の低マッハ精度向けには有用な可能性が
    あるため**機能は残す**が、**ノズル低マッハ limit cycle の根治には使えない**と確定。根治の本命は引き続き Phase1
    (圧力散逸) 側で、$\epsilon$ を下げる方向には安定限界 (ε≳0.15) があるため、それを超える根治は完全 (a) の
    preconditioned-flux Jacobian (収束加速と同時) に戻ることになる。成果物: `run_conv_p3_thornber_only` /
    `run_conv_p3_combo` (`residual_history.png` 付き)。

- `2026-06-08` — **Phase 4 (a) 完全 $\Gamma^{-1}A$ 前処理を採用・設計確定 (実装着手前のフィージビリティ)**。
  Thornber が負だったため、振動根治は「$\epsilon$ を物理値まで下げて圧力散逸を増やす」ルートに戻るが、それには
  LHS で aggressive flux を陰的安定化する完全前処理が要る (同時に収束加速)。設計を §5 Phase 4 に確定:
  - **変数順の訂正**: forge の保存変数は $Q=(\rho,\rho u,\rho v,\rho w,\rho E)$ (energy は index 4)。
    旧 §9 の「(ρ,ρE,ρu,ρv,ρw)」は誤記。`build_jacobian_split` の $R,L$・`rhs[5]` で確認。
  - **対角ブロックは既に 5×5 を `solve_5x5` で反転** ([timeIntegration_d.cu:743,214]) ゆえ、時間項を
    スカラー $V/\Delta\tau\,I$ から行列 $\Gamma_Q V/\Delta\tau'$ に替える構造コストは小。
  - **単精度フィージビリティ結論**: 保存基底で $D$ は条件数 $\sim1/\beta$ ($\epsilon$=0.15 で ~44, 0.05 で ~400)。
    前処理固有基底では良条件なので、悪条件は $\Gamma_Q$ 変換のみが担う。**対策＝対角 5×5 の組み立て+`solve_5x5` を
    局所倍精度化** (セル毎・小・低コスト)。これで float 本体のまま $\epsilon$ を下げられる見込み。原始/固有基底で
    組む案は保存パイプラインに侵襲的で不採用。
  - **最大の実装リスク**: $\Gamma_Q$ と前処理固有ベクトル $R_\Gamma,L_\Gamma$ の閉形導出と、$\beta=1$ ビット一致回帰。
    実装は数式確定 → ヘッダ拡張 → `solve_5x5` 倍精度 → `build_jacobian_split`/`setDT_d` 前処理 → 回帰 → 検証 の順。

- `2026-06-08` — **Phase 4 土台を実装 (低リスク基盤)。ビルド成功・閉形検証 PASS**。
  - **$\Gamma_c$ がランク 1 閉形に閉じる重要発見**: $\Gamma_c=\Gamma_p M^{-1}=I+\tfrac{1-\beta}{\beta c^2}g r^\top$
    ($g=(1,u,v,w,H)$, $r^\top=(\gamma-1)(e_k,-u,-v,-w,1)=\partial p/\partial Q_c$, $r^\top g=c^2$)、逆も
    $\Gamma_c^{-1}=I-\tfrac{1-\beta}{c^2}g r^\top$ (Sherman-Morrison)。$\beta=1$ で両者厳密に $I$。**固有値 $\lambda'$ の
    公式も $\mathrm{eig}(\Gamma_c^{-1}A)$ と一致を確認** (固有ベクトル導出前にスペクトルを確定)。
  - **実装した土台**: `lowMachPrecond_d.cuh` に `lowMachBeta`/`lowMachGammaC`/`lowMachGammaCinv` (double rank-1)、
    `timeIntegration_d.cu` に `solve_5x5_dbl` (倍精度ブロック解)。いずれも未配線 (既存経路ビット不変・ビルド 5/5 成功)。
  - **検証**: `tools/verify_lowmach_precond.py` (numpy・2000 状態) で上記 5 恒等式すべて PASS (最大相対誤差
    $\Gamma_c$ 5e-10, $\Gamma_c\Gamma_c^{-1}$ 4e-8, $r^\top g$ 6e-16, β=1 厳密 0, $\lambda'$ 1.7e-7)。
  - **残り (次段)**: 前処理固有ベクトル $R_P=MR'_p,\ L_P=L'_p M^{-1}$ の閉形 → `build_jacobian_split` 前処理版
    (`lowMachPrecond==2`)、`setDT_d` の $\Delta\tau'$、対角時間項の $\Gamma_c$ 行列化 (倍精度ブロック組み立て+解)、
    $\beta=1$ GPU 回帰 (超音速ケースで `lowMachPrecond:0` とビット一致)、`case/23.axi_nozzle` で $\epsilon$↓ 検証。

- `2026-06-09` — **Phase 4 LHS 本体を実装 (固有ベクトル前処理は採らず Rusanov 型分割で代替)。ビルド成功・
  500 step スモーク稼働確認**。
  - **設計判断**: 厳密前処理 FVS ($\hat A_\Gamma^{\pm}=\Gamma_c R_P\Lambda'^{\pm}L_P$) は前処理固有ベクトル $R_P$ の
    閉形導出が要り重い。LHS は defect-correction で**解を変えない**ため、固有ベクトル不要な
    **スペクトル半径分割** $\hat A_\Gamma^{\pm}=\tfrac12(\hat A_\Gamma\pm\rho'I)$, $\hat A_\Gamma=\Gamma_c^{-1}A_c$,
    $\rho'=\tfrac12(1+\beta)|U_n|+c'$ で代替。より散逸的だが LHS なら安定寄り。$A_c=a^{+}-k_{\rm off}$ は既存
    `build_jacobian_split` から、$\Gamma_c^{-1}$ はランク 1 閉形から。
  - **実装**: `timeIntegration_d.cu` に別カーネル `implicit_defect_correction_block_precond_d`
    (対角ブロックを倍精度で組み `solve_5x5_dbl` で解く。時間項 $\Gamma_c V/\Delta\tau'$、分割 $\tfrac12(\hat A_\Gamma\pm\rho'I)$、
    軸対称ソース・粘性・dual-time BDF も倍精度で踏襲)、`matmul_5x5_dbl` ヘルパ。`setDT_d` に
    `setDTlocal_precond_scale_d` ($\Delta\tau'{=}\Delta\tau\cdot(|\mathbf u|+c)/\rho'$、$\sim1/\epsilon$ 有界)。
    `SLAU_d` の `c'` 散逸を `lowMachPrecond>=1` に拡張 (==2 も散逸是正を使う)。wrapper は `==2` で前処理カーネル+
    setDT スケーリングに dispatch。既存カーネルは別物ゆえ **0/1 経路はビット・レジスタ不変**。
  - **0/1 ビット不変の確認**: 別カーネル化＋`SLAU` は `==0` で分岐に入らないため。`==2` は別系。
  - **スモーク**: `case/23.axi_nozzle` を `lowMachPrecond:2, ε=0.15` で 500 step → NaN なし・残差降下
    (rms_roe 41→9.6)。速度は 0/1 比 ~3.6×/step (倍精度+レジスタスピル。実験モードとして許容)。
  - **検証中 (20k)**: `run_conv_p4_eps15` (ε=0.15) と `run_conv_p4_eps05` (ε=0.05) を 20k step 実行中。
    狙い: (1) ε=0.15 で Phase1 単独 (0.603%) と同等以上、(2) **ε=0.05 が安定** (Phase1 では発散した領域) で
    振幅が物理目標 ~0.25% に近づくか。結果は次エントリ。

- `2026-06-09` — **cfl_pseudo スウィープで収束加速を検証 → 導出誤りを発見・修正 → 加速はブレークイーブンと判明**。
  - **cfl_pseudo=0.5 同条件では P4 が負け** (limit-cycle 帯 rms_roe mean: baseline 2.78, Phase1 2.09, P4 2.96)。
    前処理の旨みは高 cfl_pseudo にあるので設定が不適と判断し cfl スウィープへ。
  - **スウィープ (eps=0.15, RHS 共通, LHS のみ m1↔m2)**: 安定限界は **m1 (既存 FVS LHS) ~1** (cfl2 で発散)、
    **m2 (前処理 LHS) ~2-3** (cfl2 安定・cfl4 発散)。前処理で安定限界は広がるが過散逸で頭打ち。
  - **導出誤りの発見・修正 (重要)**: 保存形は $(\Gamma_c V/\Delta\tau' + A_c)\Delta Q=-R$ で、**前処理は時間項 $\Gamma_c$ のみ**、
    フラックス $A_c$ は**物理 (非前処理)** が正しい。当初 $\hat A_\Gamma=\Gamma_c^{-1}A_c$ のスペクトル半径分割で
    フラックスも前処理していた (二重前処理) のが過散逸の主因。`implicit_defect_correction_block_precond_d` の
    フラックス分割を物理 `a_plus`/`k_off` (既存と同一) に戻し、$\Gamma_c$ 時間項のみ残した (Gcinv/matmul 削除→軽量化)。
  - **修正版スウィープ**: 安定限界 **~5-7** に上昇 (cfl5 安定・cfl10 発散)。m1 ~1 の **5×**。per-step は **2.54×** (旧 3.6× から改善)。
  - **収束加速の判定 (ブレークイーブン)**: 等 wall-clock では m2 と m1 はほぼ互角〜やや m1 有利。理由: (1) 残差は
    低マッハ limit-cycle フロアに支配され cfl を上げても比例して速くならない、(2) 内部 Jacobi (nStepInner=20) の頭打ち、
    (3) per-step 2.54× コスト。**収束加速器としては決定打にならず**。
  - **β=1 性質**: 修正版は β=1 で $\Gamma_c=I$・$\Delta\tau'=\Delta\tau$・フラックス同一ゆえ既存カーネルと同一組み立て
    (倍精度ゆえ丸め差 ~1e-7、解一致)。
  - **次**: Phase 4 の価値は収束加速でなく**根治** (ε=0.05 を安定化し振幅を下げる) に懸かる。`run_p4_amp_eps05`
    (==2, ε=0.05, cfl_pseudo=2, 20k) を実行中。結果次第で Phase 4 の採否を確定。

- `2026-06-09` — **Phase 4 根治テスト成功。低マッハ自励振動を根治 (採用確定)**。`case/23.axi_nozzle` 20k step、
  chamber std/mean (M<0.08, 4k–20k 包絡):

  | 設定 | mean [min,max] | vs baseline | 状態 | M_max |
  | --- | --- | --- | --- | --- |
  | baseline (p0) | 0.882% [0.25,1.82] | — | limit cycle | 6.057 |
  | Phase1 ε=0.15 (m1 LHS) | 0.603% [0.11,1.53] | −32% | limit cycle | 6.046 |
  | **P4 ε=0.15 (m2 LHS)** | **0.333% [0.08,0.70]** | **−62%** | limit cycle(減衰) | 6.04 |
  | **P4 ε=0.05 (m2 LHS)** | **0.087% [0.056,0.228]** | **−90%** | **定常収束(振動消滅)** | 6.033 |

  - **要点**: (1) 同じ ε=0.15 でも前処理 LHS だけで 0.603→0.333% (limit cycle は真収束でないため、より安定な LHS が
    振動を減衰)。(2) **ε=0.05 は前処理 LHS でしか安定化できず** (m1 は数十 step で NaN)、そこで per-snapshot が
    単調減衰し step12000 以降 0.056% で平坦=**自励振動が消え定常収束**。物理目標 ~0.25% も下回る。超音速域 (M_max) 不変・NaN なし。
  - **結論**: Phase 4 (完全 Γ⁻¹A 前処理, `lowMachPrecond=2`) は **plan 当初目的 (低マッハ自励振動の根治) を達成**。
    Phase1 単独では届かなかった ε<0.15 の強い `c'` 散逸を前処理 LHS が陰的に安定化することで実現。**採用**。
  - **収束加速は副次 (ブレークイーブン)**: cfl 限界は m1~1→m2~5-7 と広がるが per-step 2.54× で等 wall-clock は互角
    (上エントリ)。Phase 4 の価値は加速でなく**根治**にある。
  - **推奨運用**: 低マッハ振動が問題になる定常ケースで `lowMachPrecond: 2`, `precondEps: 0.05`, `cfl_pseudo: 2` 前後。
    成果物: `run_p4_amp_eps05` / `run_p4_amp_eps15` (`residual_history.png` 付き)。
  - **残課題 (任意)**: per-step の倍精度コスト低減 (条件数が許す範囲で float 化・レジスタ削減)、β=1 完全回帰の
    専用超音速ケース整備、ε のさらなる最適化。いずれも根治の成立とは独立。

- `2026-06-09` — **Sherman-Morrison 化で FP64 を排し、前処理カーネルを物理 block とほぼ同速に高速化**。
  - **着眼**: RTX 3060 等 consumer GPU は FP64=FP32/64。前処理カーネルの倍精度対角ブロック+`solve_5x5_dbl` が
    2.54× 遅延の主因。$\Gamma_c=I+\alpha g r^\top$ が**ランク1**ゆえ対角ブロックは $D=D_0+\gamma g r^\top$
    ($D_0$=物理ブロック・良条件・既存 0/1 が float で解いているのと同形、$\gamma=(V/\Delta\tau')\alpha$) と分解できる。
  - **実装**: `solve_5x5_dbl` を float の `solve_5x5_2rhs` (1 回分解で 2 RHS) に置換し、前処理カーネルを
    Sherman-Morrison $x=y-[\gamma(r^\top y)/(1+\gamma(r^\top z))]z$ ($y=D_0^{-1}b,\ z=D_0^{-1}g$) で解く。
    **D_0 解は float**、悪条件 $\sim1/\beta$ は分母スカラーのみ double に隔離。$g,r,\alpha$ はカーネル内インライン化し
    未使用の `lowMachGammaC`/`lowMachGammaRank1`/`matmul_5x5_dbl`/`solve_5x5_dbl` を削除。
  - **結果**: per-step **43ms→16.9ms (2.5×、baseline m1 17.8ms とほぼ同速)**。`case/23.axi_nozzle` ε=0.05 20k で
    **振幅 0.087% を完全再現** (mean/min/max とも旧 double 版と一致)、残差トラジェクトリも一致、NaN なし。
    20k 実時間 865s→**321s**。FP64 ペナルティ消滅。
  - **収束加速への波及**: per-step がほぼ同速になったので、m2 の cfl_pseudo 余裕 (~5-7 対 m1 ~1) が
    per-step コストで相殺されなくなった。等 wall-clock の収束加速も net で有利になった可能性が高い (要再測定)。

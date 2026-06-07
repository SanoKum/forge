# 低マッハ前処理 (Weiss–Smith) — 密度ベース経路の散逸スケール是正と陰解法固有値前処理

## メタ

- **area**: `time_integration` (convection 横断)
- **status**: `in_progress`  <!-- draft / in_progress / done / superseded -->
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
- **やらない**:
  - 完全 $\Gamma^{-1}A$ 前処理 (固有ベクトル $R,L$ の前処理化)。Phase 2 は固有値差し替えの最小差分を第一候補とし、必要性が出た場合のみ別途評価。
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
**今後の方向は未確定** (full preconditioned-flux Jacobian / 別系の低マッハ補正 / Phase 1 で確定、の三択。§9 末尾)。

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
- docs: `docs/convection/{theory,implementation}.md`, `docs/time_integration/{theory,implementation}.md`
  (更新済)。`docs/index.md` は構成不変につき確認のみ。
- 既存ケース: `lowMachPrecond` 既定 0 のため挙動不変。

## 8. 完了条件

- [x] 関連 `docs/convection/{theory,implementation}.md` 更新済み (収束ベース所見・Phase 2 不採用を反映)
- [x] 関連 `docs/time_integration/{theory,implementation}.md` 更新済み (固有値/setDT 前処理は不採用と明記)
- [x] Phase 1 実装・検証完了 (収束ベース: 低マッハ自励振動を ε=0.15 で −32% 減衰。§9)
- [x] Phase 2 (LHS 固有値 ＋ setDT) 実装・検証 → **棄却** (block DPLUR で有害、revert。§9)
- [x] `.github/plans/README.md` の状態を更新
- [ ] 根治 (完全 preconditioned-flux Jacobian / 別系補正) は未着手 (§9 の選択肢 a/b)。本 plan は Phase 1 を成果として確定

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

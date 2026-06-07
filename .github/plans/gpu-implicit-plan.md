# GPU 陰解法化計画

## メタ

- **area**: `time_integration`
- **status**: `in_progress`
- **related_docs**:
  - [`docs/time_integration/theory.md`](../../docs/time_integration/theory.md)
  - [`docs/time_integration/implementation.md`](../../docs/time_integration/implementation.md)
- **related_plans**:
  - [`architecture-rans-sst.md`](architecture-rans-sst.md) (SST の陰解法連成が本計画を参照)
- **created**: 不明 (既存計画につきメタを後付け)
- **owner**: `Copilot`

## 0. 改訂メモ (2026-06) — 制御フロー整理・古典 DPLUR 化・dual-time 受け入れ

Phase 1–6 の中核（`residual_history.csv` ログ基盤、block DPLUR 5×5 kernel、`dq_block_*`/`diag_block_*`/`rhs_block_*` バッファ、setDT 局所擬似時間）は**実装済み**。本改訂では「実装済みだが未検証・未実用」だった陰解法を、**整理して正式運用に乗せる**ことに焦点を移す。

本改訂のスコープ（詳細は別途プランニング文書に基づく）:

1. **制御フロー整理**: `main.cpp` `advanceOneStep` の `[&]` 巨大 lambda（`assembleCurrentState`/`refreshCommittedState`/`logResidualSnapshot`）と `frozen_scalar_implicit` 早期 return を廃し、`StepContext` 構造体＋自由関数（dispatcher / `assembleResidual` / `advanceExplicitRK` / `implicitNonlinearUpdate` / `blockDPLURSolve` / `advanceImplicitSteady`）に分解。
2. **古典 DPLUR 化**: 1 非線形更新あたり残差/Jacobian を 1 回構築し、`Q` 固定で `dq` を `nStepInner` 回 Jacobi sweep、最後に 1 回だけ commit。kernel の inline `Q+=dq` を撤去し `applyBlockImplicitCorrection`（commit）と `swapBlockImplicitCorrectionBuffers`（swap）をドライバ側で明示化。
3. **config 意味確定**: `nStepInner`=DPLUR sweep 回数、`nStepOuter`=擬似時間（定常）ステップ。
4. **frozen scalar 隔離**: `blockDPLUR==0`（scalar 対角版）は config で reject（kernel/buffer は温存、frozen scalar フェーズで再有効化）。
5. **dual-time 受け入れ**: 非定常 dual-time 陰解法は**必ず実装する**。本改訂では dispatcher 分岐（throw スタブ）・共有核 `implicitNonlinearUpdate`・残差の `addUnsteadyTimeTerm` フック・config 命名分担（`nSubIterDualTime` を後続で追加）だけ用意し、本体は後続フェーズ。
6. **検証**: 層流・定常ケースで陽解法収束解を真値に block DPLUR の収束解一致を確認（古典 DPLUR 化は意図的挙動変更のため現 block 実装とは比較しない）。

block DPLUR（簡易 Jacobian＝面 Roe 絶対値 Jacobian＋スペクトル半径粘性対角、matrix-free）が**実現済みの主モード**。§2 以降の「初回スコープ」記述はこの実現状況に読み替える。

## 1. 目的

`solver_density_cuda` に対して、GPU 上で実用的に回る陰解法基盤を段階導入する。

初回導入では、定常問題向けの pseudo-time implicit を対象とし、既存の陽的 Runge-Kutta 更新を置き換える形で実装する。

この計画は今後の拡張でも参照する基準文書とし、スカラー輸送、凝縮、化学反応などの追加物理は、この方針と整合する形で個別計画を作成する。

## 2. 現時点の実装方針

初回の本命は、次の 3 点を同時に満たす構成とする。

1. `matrix-free`
2. `近似 Jacobian`
3. `block-Jacobi / defect-correction`

これは、現在の `solver_density_cuda` が face-parallel の residual 計算を GPU 上で持っており、coloring なしの forward/backward sweep や full sparse Jacobian assembly よりも整合的だからである。

## 3. 採用方針

### 3.1 初回スコープ

- 対象は定常問題向けの pseudo-time implicit とする。
- 既存の移流・粘性 residual 計算経路は極力そのまま再利用する。
- Jacobian は明示構築せず、face 幾何と現在状態から implicit 用係数を別 kernel で再構成する。
- 線形更新は、cell ごとの対角 block を解く block-Jacobi / defect-correction 反復とする。
- 粘性項は residual には含める。implicit 係数には初回は対角強化として入れるか、必要に応じて lagged 扱いとする。
- 各保存方程式について、RMS Residual をログファイルへ出力する機能を初回スコープに含める。

### 3.2 初回で採らないもの

- full sparse Jacobian の明示構築
- coloring なしの LU-SGS / forward-backward sweep
- unsteady dual-time を含む完全連成 implicit
- Poisson 系の既存 matrix 基盤の再利用前提設計
- chemistry / scalar / condensation まで含む full coupled implicit

## 4. 設計原則

### 4.1 residual と implicit operator を分離する

既存の residual 計算は face-loop を中心に成立しているため、以下を分離する。

- residual の評価
- implicit 用係数の再構成
- `dQ` の反復解法
- residual の集計とログ出力

これにより、既存 flux kernel の責務を肥大化させずに陰解法を追加する。

あわせて、GPU と CPU 間のデータ転送は最小限に抑えることを明示方針とする。特に residual 監視や implicit 反復の観測量は、可能な限り GPU 上で集計し、host 側へ戻すのはログ出力や制御分岐に必要な少数の集計値に限定する。

### 4.4 Residual の定義とログ方針

ログ出力の対象は、各保存方程式

- `ro`
- `roUx`
- `roUy`
- `roUz`
- `roe`

に対する residual とする。

初回実装では、solver 本体が出力するのは各式の RMS Residual のみとする。

相対 Residual は Python などの後処理で算出する前提とし、solver 側 CSV には含めない。これにより出力ファイルの容量増加を抑え、GPU から CPU へ戻す値とディスク書き込み量を最小限に保つ。

Residual の代表値は各保存方程式の cell residual に対する root mean square とし、評価処理そのものは GPU 上で行う。host 側では GPU reduction 後の少数の集計値のみを受け取り、ログへ書き出す。

ログは人が読める text 形式を基本とし、step、inner iteration、各方程式の 3 種類の residual を時系列で追えるようにする。必要に応じて post process しやすい列形式に揃える。

### 4.5 Residual ログ出力仕様案

初回実装では、Residual ログは `CSV 1 本` を基本とする。

推奨ファイル名:

- `residual_history.csv`

出力場所:

- 各 run directory 直下

1 行は 1 回の residual 計測イベントに対応させ、少なくとも次のタイミングで出力する。

1. outer step 開始直後の residual
2. inner loop 開始時の residual
3. 各 inner iteration 後の residual
4. outer step 更新後の residual

初回は解析しやすさを優先し、long format ではなく wide format の CSV を採用する。

推奨列:

- `step`
- `inner`
- `phase`
- `rms_ro`
- `rms_roUx`
- `rms_roUy`
- `rms_roUz`
- `rms_roe`

`phase` は少なくとも次の値を持てるようにする。

- `outer_begin`
- `inner_begin`
- `inner_iter`
- `outer_end`

`inner` は outer begin / outer end では `-1` を許容し、inner loop 中のみ `0, 1, 2, ...` を使う。

Residual の集計値は、初回から各保存方程式のセル residual に対する RMS を基準とする。集計は GPU 上で reduction し、CPU 側へは 5 変数ぶんの RMS 値のみを戻す。

相対 Residual を後処理で算出する場合、`rel_init3_*` の基準値は計算開始時の最初の 3 step における `outer_begin` residual の平均値とする。3 step 未満では、取得済み step の平均を暫定基準とする。

また、`rel_inner_*` の基準値は、その outer step における `inner_begin` residual とする。

後処理での relative residual は

$$
\mathrm{rel} = \frac{R}{\max(R_{\mathrm{ref}}, \varepsilon)}
$$

で計算し、`\varepsilon` は 0 除算防止用の小さい定数とする。

初回は CSV 1 本に限定し、追加の整形済みサマリーファイルは作らない。必要になった場合のみ別計画で summary log を追加する。

CSV を採用する理由は次の通りとする。

- solver 側の実装が単純で壊れにくい
- Python による後処理が容易
- Git 差分や簡易確認がしやすい
- 追加依存なしで他ツールへ流し込みやすい

初回はこの CSV を読む Python 可視化スクリプトも合わせて提供し、必要な相対 Residual を後処理で再構成できる状態にする。

### 4.2 full Jacobian ではなく近似作用素を使う

解く対象は

$$
\left(\frac{V}{\Delta \tau}I + \tilde{J}\right) \delta Q = -R
$$

とする。

ここで `\tilde{J}` は厳密 Jacobian ではなく、移流の upwind spectral radius と粘性の対角強化から構成する近似作用素とする。

### 4.3 初回は全セル並列を維持する

coloring なしで逐次 sweep を入れると GPU 並列性が崩れるため、近傍 coupling は Jacobi 型の lagged neighbor で入れる。inner linear iteration は次の形で回す。

$$
D_i \delta Q_i^{m+1} = -R_i - \sum_{j \in N(i)} K_{ij} \delta Q_j^m
$$

ここで `D_i` はセル対角 block、`K_{ij}` は lagged neighbor coupling である。neighbor の `\delta Q_j` は常に 1 つ前の inner iteration の値を参照し、同一 kernel 内での latest neighbor / Gauss-Seidel 更新は採らない。

## 5. 実装フェーズ

### Phase 1: 設定面の追加

- `solverConfig` に implicit mode id、pseudo CFL、`nStepOuter`、`nStepInner`、update damping、viscous diagonal weight を追加する。
- 旧キー `nInnerLoop` と `dualTime_InnerLoop` は受理せず、`nStepInner` へ名前変更するよう設定エラーで停止する。
- 旧キー `nStep` は受理せず、`nStepOuter` へ名前変更するよう設定エラーで停止する。

### Phase 2: 作業バッファの追加

- `variables` に `dQ_old`、`dQ_new`、`rhs`、対角係数を追加する。
- `variables` に residual 集計用の作業領域と、各式の基準 residual を保持するための変数を追加する。
- 保存形は scalar または block-diagonal を優先し、face ごとの full 5x5 block 保存は避ける。

初回の近傍 coupling 実装では、`dQ_old` と `dQ_new` を 2 面化し、inner iteration ごとに pointer swap する。これにより coloring なしで Jacobi 型反復を維持する。

Residual ログ機能のために、少なくとも現在 step / inner の RMS residual を取得できるようにする。

relative residual 用の基準値は solver 本体には保持せず、初回は CSV を読む後処理側で再構成する方針とする。

### Phase 3: CUDA kernel の追加

- work buffer 初期化
- implicit 対角係数の再構成
- `rhs = -R - K dQ_old` の形成
- 1 回分の Jacobi / defect-correction 更新
- residual の RMS 集計

を担当する CUDA 実装単位を追加する。

Residual の集計は GPU 上で reduction し、ログ出力に必要な最小量だけを host 側へ戻す。

### Phase 4: 時間積分分岐の追加

- `timeIntegration_d` から新しい implicit mode を dispatch する。
- 既存 explicit RK と共存させる。

初回の `timeIntegration = 11` は、`D = V / \Delta\tau + \tilde{J}_{diag}` を主対角とし、近傍項は scalar 近似の `K dQ_old` として lagged に入れる Jacobi / defect-correction とする。これは `matrix-free` かつ GPU 全セル並列を維持しつつ、対角のみ baseline よりも近傍 coupling を 1 段追加する着地点である。

境界 ghost cell については、独立 unknown として off-diagonal へ入れない方針を維持する。一方で、wall / slip などで ghost state が interior state に従属している影響を 0 扱いしないため、暫定措置として non-periodic boundary face は追加の自己対角強化として数える近似を許容する。これは厳密な BC Jacobian ではなく、boundary Jacobian 導入までの暫定安定化策として扱う。

### Phase 5: メインループ統合

- 既存の residual assembly は維持する。
- implicit mode では stage loop を pseudo iteration として再解釈するか、専用 inner loop に分岐する。
- `roN` / `roM` の意味が崩れないように state 保存ルールを整理する。

`nStepOuter` は物理 step / 外側反復回数、`nStepInner` は pseudo-time / implicit inner 反復回数として分離する。

- 明示法では `nStepOuter` で外側ループを回し、`nStage` は `timeIntegration` から決める。
- 陰解法・定常では generic な `iteration_count` へ押し込めず、`nStepOuter` 回の非線形 outer loop と、その内側の `nStepInner` 回の DP-LUR / defect-correction sweep loop を専用に書く。
- 陰解法・定常の `timeIntegration_d_wrapper(..., loop, ...)` 1 回は 1 sweep とみなし、host 側の `iInner` loop が全セル同期付きの inner iteration を表す。
- 陰解法・非定常では `nStepOuter` で外側ループ、`nStepInner` で内部反復を回す。

初回は `\delta Q` の累積ではなく、各 inner iteration で現在反復値へ加える更新量として既存の作業配列を流用する。`dQ_old` は 1 つ前の inner iteration の更新量、`dQ_new` は現在の inner iteration の更新量とする。

implicit 更新の基本形は `Q \leftarrow Q + \delta Q` とし、`implicitRelax` は現在の `\delta Q` に直接掛ける damping とする。`roN` は将来の dual-time で物理時間基準状態を保持する用途に残す。

### Phase 6: 安定化と観測量

- residual RMS norm
- `dQ` norm
- pseudo CFL
- 対角係数の最小値

を確認できるようにする。必要なら damping を導入する。

あわせて、各保存方程式の residual をログファイルへ継続出力できるようにする。初回実装では、まずこのログ機能を優先的に導入し、陰解法更新の健全性確認に使う。

初回の具体仕様は `residual_history.csv` の出力を最小実装とし、GPU 内 RMS reduction ベースで安定して記録できることを優先する。

この CSV を入力として、RMS Residual と後処理で再構成した相対 Residual を対数軸で描画できる Python スクリプトを `solver_density_cuda/tools/` に追加する。

### Phase 7: 回帰確認

- 実装の健全性確認は `case/20.naca_ml` を基準ケースとして行う。
- 既存実装の baseline と、陰解法化後の結果について、壁面の静圧分布を比較する。
- 判定量は壁面静圧変化の比率とし、`0.1%` 未満の変化量であれば問題なしと判断する。
- `0.1%` を超えた場合は自動で許容扱いにせず、結果を整理した上でユーザーに判断を促す。
- あわせて residual ログを確認し、CSV から各保存方程式の RMS Residual が追跡でき、後処理で計算開始基準および inner loop 基準の相対 Residual を再構成できることを確認する。
- 上記の健全性確認を通過した後に、必要に応じて pseudo CFL や inner iteration 数を段階的に上げて収束性を確認する。

## 6. 将来拡張の指針

### 6.1 スカラー輸送 (SST k-ω) — 実装済 (2026-06)

SST (k/ω) を **segregated point-implicit** で実装済み（7×7 連成ではなく平均流 5×5 と分離）。消散項のみ陰化
（$\partial D_k/\partial(\rho k)=\beta^\*\omega$, $\partial D_\omega/\partial(\rho\omega)=2\beta\omega$）。移流・拡散・生産は lagged。
`ransSource_d` が `src_jac_k/ω` を出力、`applySSTPointImplicit_d` が平均流 commit 後に k/ω を更新（旧「凍結」解除）。
詳細は [`docs/turbulence/`](../../docs/turbulence/) と [`docs/time_integration/`](../../docs/time_integration/)。

凝縮・化学など他のスカラー/ソースは引き続き本方針（輸送 implicit ＋ source local implicit の分離）で拡張する。

### 6.2 凝縮計算

凝縮や source term が強く stiff になる場合は、輸送 implicit と source local implicit を分ける IMEX / operator splitting も選択肢に入れる。最初から full coupled Jacobian を組む前提にしない。

### 6.3 化学反応

反応は cell-local stiff ODE に近いため、輸送側 implicit と source 側 local implicit を分けた方が GPU と相性が良い。これも別計画として管理する。

## 7. 変更管理ルール

- 陰解法化に関する実装を行う前に、この文書を参照すること。
- 実装方針を変更する場合は、先にこの文書を更新し、その差分を明示してから実装に移ること。
- この文書の方針と異なる実装を暫定で入れる場合は、理由と暫定期間を記載すること。
- スカラー輸送、凝縮、化学反応など別テーマの計画は `.github/plans/` 配下に個別の `.md` として追加し、それぞれに対象・非対象・方針変更ルールを書くこと。

## 8. 現時点の判断

- first step は `matrix-free + approximate Jacobian + block-Jacobi / defect-correction`
- first implementation target は residual ログ機能の追加と、そのための `variables.cpp` / `variables.hpp` 側の保持量整理とする
- second step で `colored sweep` や `multicolor SGS` を検討する
- full sparse Jacobian assembly は、その後も必要性が明確な場合に限って評価する

## 9. 変更ログ

- 2026-06: 制御フロー整理（`StepContext`＋自由関数分解、frozen 分岐削除）、block DPLUR の古典 DPLUR 化（残差 1 回構築＋`nStepInner` sweep＋単一 commit、`applyBlockImplicitCorrection`/`swapBlockImplicitCorrectionBuffers` をドライバ側に明示化）、`blockDPLUR==0` の config 隔離、dual-time 受け入れ構造（dispatcher 分岐スタブ・共有核・`addUnsteadyTimeTerm` フック・config 命名分担）を実装。docs/time_integration の theory/implementation を同期更新。

- 2026-06: **block DPLUR 演算子の修正（収束化）**。上記リファクタ後の検証（`case/20.naca_ml/001.test`, 既定 SLAU/MUSCL/venkata, cfl 1.5）で block DPLUR が収束せず発散することが判明。切り分けの結果、新旧コードとも**対角・近傍に絶対値 Jacobian $|\widetilde A|$ を符号付き分割 $A^\pm$ の代わりに誤用**しており、対角が upwind 自己 Jacobian と不一致（対流中心項欠落）・近傍結合が逆符号だった（対角のみ＝頭打ち、近傍 sweep 有効＝NaN。旧 2026-05 実装も同様に発散）。正しい LU-SGS 分割 $D_i=V/\Delta\tau I+\sum_f A^{+}_f S_f$、RHS 近傍項 $-\sum_f A^{-}_f S_f \Delta\mathbf Q_{\text{nbr}}$ に置換（`build_jacobian_split` で $A^{+}$ と $-A^{-}$ を同時構築）。局所反復行列のスペクトル半径 $\rho\approx (V/\Delta\tau)/(V/\Delta\tau+\lambda^{+})<1$ を数値確認。theory/implementation を同期更新。

  検証結果 (2026-06-07, native build, `case/20.naca_ml/001.test`): 修正後の block DPLUR は発散せず収束（roe 残差 18→0.5）。pseudo CFL を上げると加速（cfl_pseudo=50 で 500 step に roe≈12、explicit は同 step で≈86）。両者を収束させた状態（explicit 40000 step roe≈0.25、implicit cfl_pseudo=20 4000 step roe≈0.5）で**壁面静圧の最大相対差 0.0552% < 0.1%**（平均 0.0053%）で一致＝陽解法と同一定常解に収束することを確認（plan Phase 7 合格）。陰解法の 4000 step 壁時計時間 25s は explicit 同 step 41s より短い。

- 2026-06: **軸対称ケースの陰解法対応**。幾何 ($r$ 重み付き `volume`/`ss`) は既に整合しており、平均流 `A⁺/A⁻` 修正だけで軸対称 block DPLUR は収束する（ソース lagged でも可）ことを確認。さらに半径運動量ソース $S_{\rho u_r}=(p-\tau_{\theta\theta})A_{\text{planar}}$ の局所ヤコビアン（圧力 $\partial p/\partial Q$ ＋ 粘性フープ $2\mu/(\rho r_{\text{eff}})$）を `implicit_defect_correction_block_d` の roUy 行対角に追加（`isAxisymmetric`/`A_planar` 引数追加）。検証 (`case/23.axi_nozzle` M4 ノズル): 陽解法収束解と壁面静圧が平均 0.02% で一致（最大 0.32% は未収束プルーム最低圧部、explicit/lagged/jac で同一＝case 残差フロア由来）、ソースヤコビアンで過渡収束 ~2 倍速・回帰なし。超音速始動の擬似 CFL 上限は case 律速（planar でも発散）で軸対称ソース律速ではない。

- 2026-06: **乱流 SST (k-ω) の陰解法化（segregated point-implicit）**。旧実装は block 陰解法中 k/ω を凍結していた。平均流 block DPLUR の commit 後、同一擬似時間ステップで k/ω を point-implicit 更新する経路を追加（`main.cpp` `implicitNonlinearUpdate`）。消散項の stiff 性のみ対角に陰化（`ransSource_d` が `src_jac_k=β*ω`・`src_jac_omega=2βω` を出力、`applySSTPointImplicit_d` が $D_\phi=V/\Delta\tau+V\partial D/\partial(\rho\phi)$ で更新、realizability floor 付き）。移流・拡散・生産は lagged。検証 (`case/23.axi_nozzle` M4 軸対称 SST ノズル, `run_0093` 複製を res_1000 から陰解法化): k/ω が凍結せず収束（roK ↓8倍, roOmega ↓14倍）、壁面静圧が陽解法収束解と **0.003%** で一致（PASS）。なお超音速プルームせん断層の影響で大域 roe 残差ノルムは高止まりするが、壁面・積分量は収束。強い圧縮性乱流では `cfl_pseudo ≲ 0.3` 推奨。

- 2026-06: **scalar 対角版 (`blockDPLUR==0`) を有効化**。古典 DPLUR 制御フローに対応（`blockDPLURSolve` を blockDPLUR で分岐、`swapScalarImplicitCorrectionBuffers` 追加、commit は `applyScalarImplicitCorrection`）、config の reject を撤廃（0/1 を受理）。対角はスカラー スペクトル半径 $D=V/\Delta\tau+\sum_f(|U_n|+c+\rho^\nu)S_f$。検証 (`case/20.naca_ml`): 収束先は block と同一（収束場保持・壁面静圧平均 0.02% 一致）だが近似が粗く**安定 cfl_pseudo が大幅に低い**（scalar ≲ 1〜2 vs block 20〜50）。supercritical 始動では block が cfl_pseudo=20/4000step/25s で roe→0.5 収束する一方、scalar は cfl_pseudo=1/12000step でも未収束・過渡オーバーシュート大（roe ピーク ~186 vs ~85）。**収束は block より遅い**。位置づけ: block DPLUR を既定、scalar 対角は 5×5 を避けたい軽量・低レジスタ用途のフォールバック。

- 2026-06: **非定常 dual-time 陰解法を実装**。`advanceImplicitDualTime` が 1 物理ステップ = 時間レベルシフト (`shiftDualTimeLevels_d`) → 擬似時間サブ反復 `nSubIterDualTime` 回（各回: `assembleResidual` + BDF 物理時間項 `addUnsteadyTimeTerm_d` (`res*=res-(V/Δt)(aQ-bQ^n+cQ^{n-1})`) + block DPLUR (対角に `aV/Δt` を `cfg.unsteadyDiagCoef` 経由で加算) + in-place commit `applyBlockImplicitCorrectionInPlace_d`) → 物理時間前進。BDF1(初回)/BDF2(以降)、SST k/ω も BDF＋point-implicit 対応。config: `nSubIterDualTime`(既定20), `bdfOrder`(既定2)。使用条件 `unsteady=1,dualTime=1,blockDPLUR=1,control=0`。検証 (`case/20.naca_ml`): 各物理ステップで擬似サブ反復が R* を ~4 桁低減（rms_roe 22.9→0.0035, ~10 反復で収束）、陽解法 CFL 上限超の Δt でも安定（Δt=2e-6 vs 陽解法限界 ~9e-7）、同一 Δt=5e-7 で陽解法時間精度解と壁面静圧一致（平均 0.006%, 最大 0.2%=RK3/BDF2 の O(Δt²) スキーム差）。

- 残: dual-time の scalar 対角版 (`blockDPLUR==0`) 対応、より厳密な非定常検証ケース（渦放出等）は後続。
# per-phase 詳細プロファイリングに基づく host 律速の局所最適化

## メタ

- **area**: `architecture`
- **status**: `in_progress`
- **related_docs**:
  - `docs/architecture/overview.md`
- **related_plans**:
  - `architecture-residual-monitor-async.md` (残差モニタの device 常駐化; 本 plan の前段)
- **created**: `2026-06-21`
- **owner**: `self`

## 1. 目的

残差モニタ device 常駐化後も、定常 implicit は wall ~4.79ms/step に対し GPU 実働 ~2.87ms/step で
**~40% が host/launch オーバーヘッド**として残る。「全 wrapper の `cudaDeviceSynchronize()` 一括撤廃」は
過去に効果が薄かった (= 大半が memory-bound の GPU 実働か、特定箇所の同期に局在しているため)。
そこで **盲目的な sync 撤廃をやめ、phase/カーネル単位で "詳細に" 計測してボトルネックを局在化し、
箇所ごとに最小リスクの打ち手を当てる**。各打ち手は A/B 実測で効果を数値化してから採用する。

## 2. スコープ

- **やる**:
  - (A) **計測の精緻化**: 既存 `RuntimeProfiler` (cudaEvent, host-stall gap も含む) の粒度を上げ、
    各 phase の「cudaEvent 時間」と「ncu 純カーネル時間」を突合して **同期/launch オーバーヘッド = event − kernel**
    を箇所別に定量化する。
  - (B) **データ駆動の個別最適化** (効果順, A/B で検証):
    1. `setDT`: `thrust::max_element` の per-step host 同期と 1 step 2 回呼び出しを削減
       (max を device 常駐化 / `max cfl` 表示・dt 適応を N step ごと / 冗長呼び出し統合)。
    2. 小カーネル群 (BC・fill_limiter・setCFL 系) の融合や launch 削減。
    3. 必要なら **CUDA Graph** で 1 step のカーネル列を capture/replay し launch+同期を一括削減。
  - (C) 「過去に sync 撤廃が効かなかった理由」を A/B でデータとして示す (仮説の検証)。
- **やらない**:
  - 全 wrapper の `cudaDeviceSynchronize()` 一括撤廃 (本 plan は局所・計測駆動。全廃は CUDA Graph 化に包含)。
  - 数値スキーム自体の変更 (SLAU/limiter/SST のアルゴリズム)。
  - 案3 (残差二乗和のカーネル相乗り; 前 plan の残課題として別扱い)。

## 3. 関連 docs と前提

- 計測値の初期スナップショット (case/36 run_0050, solid 79.4k cell, 500 step, FORGE_PROFILE=1):
  time_integr 15.8% / **set_dt 13.6% (avg 0.397ms×2、ncu 純カーネル ~0.11ms = 同期過大)** / limiter 11.2% /
  convective 8.6% / turbulence 8.2% / viscous 5.5% / gradient 5.1% / bconds・gas・dep 他 ~12%。
- `RuntimeProfiler` (`main.cpp`) は `FORGE_PROFILE=1` で有効、`measureCuda` の cudaEvent はカーネル間
  GPU アイドル (host 同期待ち) も含むため、event≫純カーネルの phase が同期律速の局在点。
- 各 GPU wrapper は末尾で `gpuErrchk(cudaDeviceSynchronize())` を呼ぶ (例: `setDT_d.cu`)。
- `setDT_d_wrapper` は `thrust::max_element` で `max cfl` を host へ取り (同期)、`printf` し、
  `dtControl==1` 時に dt を適応。定常経路では 1 step に 2 回呼ばれる
  (`implicitNonlinearUpdate` と `advanceImplicitSteady`)。

## 4. 設計方針 (計測 → 局在 → 打ち手)

### Phase A: 計測の精緻化
- `ProfileSection` を細分化: `SetDt` を「setCFL 系カーネル」と「max_element+host 同期」に分離計測。
  `blockDPLURSolve` / `applyBlockImplicitCorrection` / `applySSTPointImplicit` を個別 section 化。
- 各 hot section について **(event 時間) と (ncu 純カーネル時間)** を表で突合し、
  overhead = event − kernel を算出。overhead/step の降順でボトルネックを確定。
- 補助 A/B: `FORGE_PROFILE` off の素の wall も併記 (profiling 自体の同期で wall が増える分を分離)。

### Phase B: 打ち手 (確定した overhead の大きい順に、各々 A/B)
- **setDT**: ① `dtControl==1` の dt 適応と `max cfl` 表示を毎 step ではなく N step ごとにする
  (`monitorInterval` 等, 既定で従来挙動)。② max を device に置き必要時のみ host 読み。
  ③ 定常経路の 2 回呼び出しを 1 回へ統合できるか確認。
- **小カーネル融合 / launch 削減**: ncu の per-step 列で 1-3µs の極小カーネル (fill_limiter×5, BC 群) を
  融合候補として評価。
- **CUDA Graph (任意)**: 1 step のカーネル列が静的 (定常 implicit) なら `cudaGraph` capture/replay で
  launch+同期を一括削減。効果が大きければ採用、複雑さに見合わなければ見送り (A/B で判断)。

各打ち手は「残差 CSV 回帰 (前 plan の手順) + `check_convergence` + wall/step」で採否を決める。

## 5. 実装ステップ

1. (A) `RuntimeProfiler` の section 細分化 + setDT/DPLUR の内訳計測を追加 (`main.cpp`)。
2. (A) FORGE_PROFILE 計測 + ncu 突合表を作り overhead を局在化 (本 plan の §3 を更新)。
3. (B-1) setDT の host 同期削減を実装し A/B 実測。
4. (B-2/B-3) 残る overhead に応じ 小カーネル融合 / CUDA Graph を評価・実装。
5. 各段で残差 CSV 回帰 + 収束 + 速度を確認。

## 6. 検証

- **ビルド**: build-native。
- **計測**: case/36 run_0050 (solid 79.4k, 500-2000 step, GPU 単独)。phase 別 overhead 表 + wall/step。
- **回帰**: 残差 CSV が従来と整合 (列・行・solver 非決定性ノイズ床内)、`check_convergence` VERDICT 不変。
- **判定**: 各打ち手の wall/step 短縮を数値で提示。効果が誤差内なら不採用としその旨記録。

## 7. 影響範囲

- `main.cpp` (RuntimeProfiler 細分化), `cuda_forge/setDT_d.cu` (+`.cuh`), 必要なら `solverConfig` (報告間隔),
  CUDA Graph 採用時は step driver。
- docs: `docs/architecture/overview.md` の §8.4 近傍に計測・最適化方針を追記。

## 8. 完了条件

- phase 別 overhead 表で律速が局在し、採用した打ち手の効果が数値で示されている。
- 残差 CSV 回帰 OK・収束挙動不変。plan status=done + 変更ログ、`.github/plans/README.md` 同期。

## 変更ログ

- 2026-06-21: plan 作成 (draft)。初期計測 (run_0050, 500 step) で set_dt の同期過大 (event 0.40ms ≫ ncu 0.11ms) を検出。
- 2026-06-21: **Phase A 完了** (新バイナリ, case/36 run_0050 相当 solid 79.4k cell, GPU 単独)。
  `RuntimeProfiler` (FORGE_PROFILE=1, 500 step; cudaEvent は同期由来 GPU idle も含む) と ncu per-kernel (純カーネル時間)
  を突合し **overhead = prof − ncu = host 同期/launch オーバーヘッド** を phase 別に定量化。コード変更なし (既存計測で取得)。

  | phase | prof(ms/step) | ncu(ms/step) | overhead | overhead share |
  |---|--:|--:|--:|--:|
  | **SetDt** | 0.716 | 0.206 | **0.510** | **28%** |
  | TurbulenceModel | 0.417 | 0.162 | 0.255 | 14% |
  | TimeIntegration (DPLUR×5) | 0.882 | 0.661 | 0.221 | 12% |
  | CalcGradient | 0.280 | 0.123 | 0.157 | 9% |
  | ApplyBconds | 0.155 | 0.026 | 0.129 | 7% |
  | UpdateInner | 0.153 | 0.035 | 0.117 | 6% |
  | ConvectiveFlux | 0.475 | 0.378 | 0.097 | 5% |
  | Limiter | 0.626 | 0.529 | 0.097 | 5% |
  | ViscousFlux | 0.308 | 0.230 | 0.078 | 4% |
  | GasProperties/DependentVars/Ducros/UpdateOuter | — | — | 0.16 | 8% |
  | **合計** | 4.35 | 2.53(+残差0.07) | **~1.82** | 100% |

  **所見**:
  - **#1 SetDt が突出 (0.51ms/step, overhead の 28%)**。prof 0.716 ≫ ncu 0.206 で、差は `thrust::max_element` の
    D2H 同期 + `printf("max cfl")` の host 読み + 末尾 `cudaDeviceSynchronize`×2 を **1 step 2 回**やっているため。
    これは「ただの sync」ではなく **必須の host 読み出し** が毎 step 入っているのが本質 → device 常駐 max + 報告/dt適応の間引きで除ける。低リスク。
  - **残り (~1.3ms) は ~12 phase に薄く分散**した per-kernel `cudaDeviceSynchronize` 由来 (Turb 0.26 / DPLUR 0.22 /
    Grad 0.16 / BC 0.13 / …)。各 phase 〜0.04–0.26ms で、**個別 sync 削除では小さい**。→ 過去「sync 撤廃が効かなかった」
    のはこの分散構造のため (数本消しても全体 1.8ms のごく一部しか取れない)。**まとめて取るには CUDA Graph で
    1 step のカーネル列 (静的) を capture/replay** し launch+同期を一括除去するのが筋。
  - GPU 実働 (ncu) ~2.6ms/step に対し overhead ~1.8ms/step (wall ~4.4–4.8)。残差モニタは 0.074ms/step まで縮小済み (前 plan)。
  - **Phase B の優先度**: ① SetDt の host 読み出し除去 (低リスク, ~0.4ms/step 見込み) → ② CUDA Graph (高リスク高リターン,
    分散 overhead ~1.3ms/step の大半) → ③ 小カーネル融合 (BC/grad)。各々 A/B 実測で採否。
- 2026-06-21: **Phase B-1 完了 (SetDt host 読み出し除去)**。
  - **実装 (最終形)**: `setDT_d_wrapper` の制御を **`adaptDt` と `printCfl` の 2 引数に分離** (「dt を変える」と「出力する」を
    束ねない)。host 読み出し (`thrust::max_element` の D2H+同期) は `adaptDt&&dtControl==1 || printCfl` のときだけ発生。
    中間 `cudaDeviceSynchronize`×2 + 末尾の冗長 sync を撤去 (同一 default stream で順序保証)。
  - **出力間隔の統一**: `residualFlushInterval`/`cflReportInterval` を廃し **`monitorInterval` に統合** (残差 CSV flush と
    `max cfl`/`dt` console 出力の共通頻度。既定 1=従来)。dt 適応はこれと**独立**。
  - **dt 適応の方針 (経路別)**:
    - 定常 implicit: `dt_local=cfl_pseudo·dx/λ` で cfg.dt が打ち消され dt 適応は診断のみ (解に不影響) → adaptDt も printCfl も
      `monitorInterval` ごと。
    - explicit: cfg.dt が時間前進に効くため **dt 適応は毎ステップ** (`adaptDt=true`)、表示のみ間引く。
    - dual-time: pseudo/physical とも CFL に基づき時間を変えうる設計なので **dt 適応は毎サブ反復** (`adaptDt=true`、
      dtControl==1 で cfg.dt 適応)、表示のみ間引く。**※当初 dual-time を「dtControl=0 固定で適応不要」と誤って adaptDt=false に
      していたのを修正** (pseudo の CFL ベース dt_local は setCFL カーネルが無条件に毎回計算するため元から有効)。
  - **重要バグの修正**: 初版で explicit と steady の関数末尾が酷似し throttle が誤って `advanceExplicitRK` 側へ入っていた
    (explicit は unsteady で cfg.dt を物理時間前進に使うため致命的) → 正しい steady 末尾へ移設。
- 2026-06-21: **Phase B-2 ステップ① 完了 (per-kernel sync の peek 化)**。
  - **監査**: cuda_forge 全カーネルが default (legacy) stream、`cudaMemcpyAsync`/明示 stream 無し (residualMonitor の
    memsetAsync は同 stream で順序付き＋flush で同期)。→ per-kernel `cudaDeviceSynchronize` は正しさに不要
    (same-stream 順序保証、host 読みは cudaMemcpy/thrust が intrinsic に同期)。撤去で失うのは実行エラーの即時・行単位検知のみ。
  - **実装**: `cudaWrapper.cuh` に env ゲート `forgeKernelSyncEnabled()` (既定 false=peek化) と
    `gpuErrchkKernelSync()` を追加。従来の `gpuErrchk(cudaDeviceSynchronize())` 計 54 箇所 (cuda_forge/*.cu + boundaryCond.cpp)
    を一括置換。`gpuErrchk(cudaPeekAtLastError())` は残置 (起動エラー検知は維持)。**`FORGE_KERNEL_SYNC=1` で従来同期に復帰**
    (デバッグ用)。
  - **検証** (case/36 run_0050 solid 79.4k, 2000 step, GPU 単独):
    - 速度: sync ON 8.05–8.27s → **sync OFF (peek化) 6.30s (−22%)**。per-step 4.05→3.15ms。
    - 正しさ: 残差 CSV sync-OFF vs ON = 0.128 (非決定性ノイズ床 0.13–0.17 内)。run 完走・GPUassert/NaN 無し。
    - **累計**: 開始時 OLD (毎step残差+detectNaN per-var, sync毎カーネル) 13.38s → **6.30s (−53%)**。
  - **Graph (B-2 ②) の伸びしろ確定**: peek化後 3.15ms/step に対し GPU 実働床 ~2.6ms/step → Graph の残り上限は
    ~0.55ms/step (per-launch CPU コスト) = 最大 −17%・現実 −10〜15%。**peek化が回収可能分の大半を取ったため、
    Graph は高コスト高リスクの割に増分小 → 現時点では見送り推奨**。
  - 注: FORGE_PROFILE の `measureWall` セクションは sync OFF だと kernel 実行を待たず過少計測になる (profiling 時のみの
    注意点・正しさには無関係)。`measureCuda` は cudaEvent 同期で正確。
  - 残: Graph (②) / 小カーネル融合 (③) / device-resident dt (explicit・dual-time 向け) は未着手。
  - **検証** (case/36 run_0050 solid 79.4k, 2000 step, GPU 単独):
    - 速度: setDt 前 9.37s → interval=1 (sync 撤去のみ) 9.09s → **interval=50 8.14s (setDt 前比 −13%)**。
      Phase A 予測 (setDt overhead 0.51ms/step) どおり ~0.5ms/step 回収。開始時 OLD 13.38s からは **累計 −39%**。
    - 正しさ: residual CSV の差は 2000 step の solver 非決定性ノイズ床 (rms_roOmega で同一バイナリ 2 回 0.13–0.17)
      と同等 (interval=50 vs interval=1 = 0.18) → 間引きは非決定性以上の差を生まない。
    - **非定常安全性 (実証)**: explicit unsteady (timeIntegration=1, unsteady=1) で `monitorInterval` 1 vs 50 を比較 →
      **残差 CSV 最大相対差 4e-5 ≈ 0** (dt は毎ステップ適応され結果不変)、`max cfl` print は 11 vs 2 (表示のみ間引き)。
      dt 適応と出力の分離が正しく効くことを確認。
  - **残課題**: Phase B-2 (CUDA Graph で分散 overhead ~1.3ms/step を一括回収) / B-3 (小カーネル融合) は未着手。
    定常末尾 setDt の冗長な setCFL カーネル再計算 (~80µs/step) も削減余地 (別途)。

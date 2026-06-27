# 残差モニタの device 常駐化 (per-step host 同期の除去)

## メタ

- **area**: `architecture`
- **status**: `done`
- **related_docs**:
  - `docs/architecture/overview.md`
- **related_plans**:
- **created**: `2026-06-21`
- **owner**: `self`

## 1. 目的

定常 implicit (timeIntegration=11) の 1 step は GPU 実働が ~2.9ms/step なのに wall ~5.8ms/step
で、約半分が host/launch オーバーヘッドだった (case/36 run_0050 の ncu + 実測)。主因は
**毎 step の残差 RMS 計算 (`thrust::transform_reduce` が変数ごとに host スカラを返す = 変数ごとに
`cudaStreamSynchronize`)** と **detectNaN (`thrust::any_of` が変数ごとに sync)**。
残差情報の毎 step 粒度を保ったまま、これらの per-step 同期を除去して host 律速をほぐす。

## 2. スコープ

- **やる**:
  - 残差 RMS を **fused 1 カーネルで全保存量一括 sum-of-squares → device バッファに格納** (host に返さない=同期しない)。
  - device バッファを `residualFlushInterval` step ごとに 1 回だけ D2H 転送して CSV へ一括書き出し。CSV の列・行構成・値は不変 (毎 step・3 行/step・`rms_*`/`rms_dq_*`)。
  - detectNaN を **fused 1 カーネルで device int フラグ化**し、`detectNaNInterval` step ごとに 1 回だけフラグを host 読み出し。検知時のみ従来の per-var 特定→ダンプ→停止。
  - `solverConfig` に `residualFlushInterval` (既定 1=従来挙動), `detectNaNInterval` (既定 1) を追加。
- **やらない**:
  - 案3 (残差二乗和を assembleResidual / timeIntegration カーネルへ相乗りさせ launch ゼロ化) — 案2 実測後に効果が残れば別途。
  - setDT の `max cfl` reduction (`thrust::max_element`) の同期除去 — 別 launch 系統、本 plan 対象外。
  - 全 wrapper の `cudaDeviceSynchronize()` 撤廃 (大規模・別テーマ)。
  - CPU 経路 (`gpu=0`) は従来通り即時計算・即時書き出しのまま (低リスク維持)。

## 3. 関連 docs と前提

- `docs/architecture/overview.md` に host-device 同期と残差モニタの節を追加。
- 現状の残差経路: `main.cpp` の `gatherResidualSnapshot` → `ResidualCsvLogger` (即時 writeRow + per-row `flush`)、
  reduction 本体は `cuda_forge/residualMonitor_d.cu` の `gatherVariableRms_d` (変数ごと `thrust::transform_reduce`)。
- 確認済み事実: 定常経路は gather 1 回/step・slot 共有 3 行、`rms_dq_*` は常に 0 (CorrectionSnapshot 既定)。

## 4. 設計方針

### device 側 reducer (`residualMonitor_d.cu`)
- `DeviceResidualReducer { int nVar; geom_int nCells; const flow_float** d_ptrs; double* d_sumsq; }`。
  `d_ptrs` は対象変数 (`res_ro`..`res_roOmega`) の device ポインタ配列 (init 時に 1 度だけ host→device コピー)。
- fused kernel `residualSumSq_d(nVar, nCells, d_ptrs, d_sumsq)`: grid-stride + shared-mem block 縮約 →
  ブロックごと 1 回だけ `atomicAdd(&d_sumsq[v], partial)`。累算は **double** (float thrust 版より安定;
  値は ~1e-4 相対で一致)。事前に `cudaMemsetAsync(d_sumsq,0)`。
- `residualFinalize_d`: `d_buf[slot*nVar+v] = sqrt(d_sumsq[v]/nCells)`。
- `reduceResidualToSlot(reducer, d_buf, slot)` = memset+2 kernel、**`cudaDeviceSynchronize` を呼ばない** (peek のみ)。
- detectNaN 用: 同型 reducer を NaN 名 (`ro,roUx,..,P[,roK,roOmega]`) で作り、
  `scanNonFiniteToFlag(reducer, d_flag)` (fused any-nonfinite → atomicOr) と
  `downloadFlag` (sync; `detectNaNInterval` 毎のみ呼ぶ)。

### logger (`main.cpp ResidualCsvLogger`)
- GPU 経路のみ: `ResidualSnapshot` を slot 参照に。logger が device buffer
  (`capacity=residualFlushInterval`, 実確保は `(capacity+margin)*nVar`) と
  行記述子 `rows_{step,inner,phase,slot}` を保持。
- gather (`logGpu`) は slot へ async reduce のみ。`logOuterEnd` で `device_count_>=capacity_` なら
  **step 境界で flush** (D2H 一括 → 記述子を CSV 出力 → reset)。`flush` の `stream_.flush()` は batch 末 1 回。
- destructor / main ループ後 / detectNaN 停止前に最終 flush。
- CPU 経路 (`gpu=0`) は現行コードのまま。

## 5. 実装ステップ

1. `input/solverConfig.hpp` / `.cpp`: `residualFlushInterval` (既定1), `detectNaNInterval` (既定1) 追加。
2. `cuda_forge/residualMonitor_d.cuh` / `.cu`: `DeviceResidualReducer`, fused reduce/finalize kernel,
   `reduceResidualToSlot`, device buffer alloc/download, detectNaN fused scan を追加。既存 `gatherVariableRms_d`
   (CPU 経路と互換のため) は温存。
3. `main.cpp`: `ResidualCsvLogger` を GPU バッファリング対応に改修。`logResidualSnapshot` を GPU/CPU で分岐。
   `checkNonFiniteAndHalt` を fused フラグ + `detectNaNInterval` 間引きに改修。
4. build-native 再ビルド。

## 6. 検証

- **ビルド**: build-native フルリビルド成功。
- **数値一致**: 同一ケース (case/36 run_0050 入力, 2000 step) を旧/新バイナリで実行し
  `residual_history.csv` の `rms_*` を全列比較 (相対差 ≤ ~1e-3、トレンド一致)。`rms_dq_*`=0 維持、行数・phase 不変。
- **性能**: `residualFlushInterval` 大 + `detectNaNInterval` 大で wall/step 短縮を実測 (GPU 単独)。
  detectNaN=1 の従来比短縮を再確認。
- **収束判定**: `tools/check_convergence.py` の VERDICT を従来 run と突き合わせ (挙動不変)。

## 7. 影響範囲

- `input/solverConfig.{hpp,cpp}`, `cuda_forge/residualMonitor_d.{cu,cuh}`, `main.cpp`。
- `solverConfig.yaml` に任意キー追加 (未指定なら既定で従来挙動・ビット非互換は reduction 桁落ちのみ)。
- docs: `docs/architecture/overview.md`, `docs/index.md` (必要なら)。

## 8. 完了条件

- 新バイナリで残差 CSV が従来と整合 (列・行・値) し、`residualFlushInterval`/`detectNaNInterval` 増で wall/step が短縮。
- plan status を done 化し変更ログを追記、`.github/plans/README.md` に追加。

## 変更ログ

- 2026-06-21: plan 作成 (draft→in_progress)。case/36 run_0050 の ncu/実測 (GPU busy 2.87ms/step vs wall 5.77ms/step、detectNaN +1.45ms/step) を根拠に着手。
- 2026-06-21: 実装完了 (status→done)。
  - **実装**: `cuda_forge/residualMonitor_d.{cu,cuh}` に `DeviceResidualReducer` + fused `residualSumSq_kernel`
    (全保存量を 1 カーネルで block 縮約・double 累算) + `residualFinalize_kernel` + 非有限 fused scan
    (`nonFinite_kernel`→device int フラグ) を追加。`main.cpp` の `ResidualCsvLogger` を GPU バッファリング化
    (slot 参照行を buffer→`residualFlushInterval` step ごとに 1 回 D2H+一括書き出し)。`checkNonFiniteAndHalt` を
    fused フラグ + `detectNaNInterval` 間引きに改修。`solverConfig` に `residualFlushInterval`/`detectNaNInterval`
    (既定 1) 追加。CPU 経路は従来通り。
  - **ビルド**: build-native クリーンリビルド成功 (exit 0)。
  - **数値検証**: 旧バイナリ vs 新バイナリ (case/36 run_0050 入力, 79.4k cell, 2000 step)。CSV 列・行構成 (6000 行・
    3 行/step・`rms_dq_*`=0) 不変。rms 値の旧↔新差は **solver 固有の非決定性ノイズ床と同等** (200 step:
    old↔old 4.05% vs old↔new 4.36%; forge solver は flux/勾配 gather の atomicAdd で run-to-run 非決定的で、
    場も old↔old がビット不一致)。残差ノルムは read-only 診断のため解に影響なし。新は double 縮約で旧 float より高精度。
  - **性能 (2000 step, GPU 単独)**: OLD(detectNaN=1) **13.38s** → NEW(flush=1,NaNint=1) **9.92s (−26%)** →
    NEW(flush=200,NaNint=200) **9.58s (−28%)**。**新コードは残差監視+detectNaN を有効にしたまま、旧コードで
    detectNaN を完全 OFF にした場合 (10.42s) より速い**。主因は fused 縮約で per-step の launch+同期が 14→2 に減少。
  - 記録 run: `case/36.passive_pseudoshock_control/run_0050_solid_prof/` (flush=200/NaNint=200, `residual_history.png`,
    `CONVERGENCE_VERDICT.txt`)。
  - **残課題 (別 plan 候補)**: ① 案3 (残差二乗和を assembleResidual/timeIntegration へ相乗りし launch ゼロ化)、
    ② setDT の `max cfl` reduction (`thrust::max_element`) の per-step 同期、③ 各 wrapper の `cudaDeviceSynchronize()`
    全廃 (最大の host 律速要因)。本 plan ではいずれも対象外。
- 2026-06-21 (後続): 本 plan で導入した `residualFlushInterval` は
  [architecture-perphase-profiling-hotspot.md](../active/architecture-perphase-profiling-hotspot.md) (Phase B-1) で
  **`monitorInterval` に改名・統合**した (残差 CSV flush と `max cfl`/`dt` の console 出力頻度を 1 つの間隔に統一)。
  挙動・既定 (1) は不変。`detectNaNInterval` はそのまま。

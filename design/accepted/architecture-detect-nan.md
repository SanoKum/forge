# detectNaN — NaN/Inf 検知・自動ダンプ診断オプション

## メタ

- **area**: `architecture`
- **status**: `done`
- **related_docs**:
  - `.github/forge-solver-settings.md`（`detectNaN` 設定リファレンス）
- **created**: `2026-06-11`
- **owner**: `CFD Dev`

## 1. 目的

発散の発生ステップと場所を特定しやすくするための診断オプション。毎ステップ終端で保存量
（と圧力）の NaN/Inf を検査し、検出したらその時点の解を h5 にダンプして即停止する。
理論変更を伴わない純粋な診断機能のため `docs/<area>` の theory/implementation は新設しない。

## 2. スコープ

- **やる**: `time.deltaT.detectNaN`（既定 0）を追加。`1` のとき毎ステップ終端で
  `ro,roUx,roUy,roUz,roe,P`（RANS 時は `roK,roOmega`）を内部セルにわたり検査し、
  非有限値があれば `res_nan_<step>.h5`/`.xmf` を強制ダンプして `exit(1)`。
- **やらない**: 検出後の自動リカバリ（CFL 縮小・ロールバック等）。原因セルの逐次特定ログ。
  ghost/境界セルまで含めた検査（内部セルのみ）。

## 3. 関連 docs と前提

- `.github/forge-solver-settings.md` に `detectNaN` 節を追加済み。
- AGENTS.md の「NaN / 発散チェック (必須)」運用を実行時に支援する位置づけ。

## 4. 設計方針

- 既定 0 で検査を一切行わない（通常実行はビット不変・性能影響なし）。
- 検査は GPU リダクション `thrust::any_of` + `IsNonFinite` 関手で実装し、最初に該当した
  変数名を返す（`hasNonFiniteCellValue_d`）。CPU 経路（`gpu==0`）はホスト配列を走査。
- ダンプは既存の h5/XDMF 書き出し本体を `writeSolutionH5_XDMF(..., prefix)` に切り出し、
  出力間隔ガードを無視する `dumpSolutionH5_force(..., "res_nan_")` を追加して再利用。
- 呼び出し点は `advanceOneStep` の終端 1 箇所に集約し、explicit / implicit steady /
  dual-time すべてを一様にカバーする。

## 5. 実装ステップ

1. `input/solverConfig.hpp` / `input/solverConfig.cpp`: `int detectNaN=0` を追加し
   `time.deltaT.detectNaN` を optional パース。
2. `cuda_forge/residualMonitor_d.cu` / `.cuh`: `IsNonFinite` 関手と
   `hasNonFiniteCellValue_d(msh,var,names,offending)` を追加。
3. `output/output.cpp` / `output.hpp`: 書き出し本体を `writeSolutionH5_XDMF(...,prefix)` に
   切り出し、`outputH5_XDMF`（従来ガード付き）と `dumpSolutionH5_force`（ガード無し）を提供。
4. `main.cpp`: `checkNonFiniteAndHalt(StepContext&)` を追加し、`advanceOneStep` 終端で
   `cfg.detectNaN==1` のとき呼ぶ。検出時は stderr 出力 → `dumpSolutionH5_force` → `exit(1)`。

## 6. 検証

- **ビルド**: `tools/build_native_wsl.sh`（native, sm_86）で警告のみ・エラー無しでリンク成功。
- **検証ケース**: `case/05.sod_shock_tube`（CPG）を複製した使い捨て run 2 本（検証後削除）。
  - 発散ケース（`dt` 過大で CFL≫1）: step 2 終端で `ro` の NaN を検出し
    `res_nan_2.h5`/`.xmf` をダンプして `exit(1)`。ダンプ h5 に NaN が含まれること、
    XDMF が `res_nan_2.h5` を正しく参照することを確認。
  - 安定ケース（通常 dt + `detectNaN:1`）: 300 step 完走・`exit(0)`・`res_nan_*` 無し・
    残差に NaN 無し（誤検出しないこと）。
- **判定基準**: 発散時に停止+ダンプ、安定時に無影響。既定 0 の挙動はビット不変。

## 7. 影響範囲

- 触るファイル: `input/solverConfig.{hpp,cpp}`, `cuda_forge/residualMonitor_d.{cu,cuh}`,
  `output/output.{cpp,hpp}`, `main.cpp`。
- 既存ケースは `detectNaN` 未指定で既定 0 → 挙動・性能とも不変。
- ドキュメント: `.github/forge-solver-settings.md` に節追加。

## 8. 完了条件

- [x] 理論変更なし（theory/implementation 新設不要）
- [x] 実装・検証完了（§6）
- [x] `.github/plans/README.md` に追記
- [x] `status: done`

## 9. 変更ログ

- `2026-06-11` — 初稿。実装・native ビルド・sod_shock_tube での発散/安定 2 ケース検証まで完了し `done`。

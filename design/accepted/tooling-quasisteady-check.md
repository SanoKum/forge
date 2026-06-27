# 準定常確認ツール `check_quasisteady.py` と運用ルール

## メタ

- **area**: `architecture` (検証ツール / 運用フロー)
- **status**: `done`
- **related_docs**:
  - `AGENTS.md` 「準定常確認 (必須)」節 (本ルールの正本)
  - `solver_density_cuda/tools/check_quasisteady.py` (実体)
- **related_plans**:
  - 既存 `check_convergence.py` (残差収束) と対をなす
- **created**: `2026-06-27`
- **owner**: `CFD Dev`

## 1. 目的

「残差が下がった (またはプラトー)」=「報告してよい」ではない。**報告する派生量そのもの**
(衝撃位置・上下非対称・CL/CD・massflux・推力・peak μt 等) が**定常化 (頭打ち) したか**は別問題で、
**過渡ピーク ≠ 飽和値**である。これを怠り過渡ピークの量を定常値として誤報告した事例 (case/36 で
visc=0 Euler の非対称が過渡 0.25→減衰 0.05 なのに ~12k step の 0.25 を定常偏りと誤判定) を受け、
**量の定常化をツールで判定し VERDICT を必須化**する。`check_convergence.py` (残差) と対をなす。

## 2. スコープ

- **やる**: (a) 全 `res_*.h5` スナップショット時系列から対象量を計算し末尾の頭打ちを判定する
  `check_quasisteady.py`、(b) AGENTS.md「準定常確認 (必須)」ルール (量報告に VERDICT 必須)。
- **やらない**: probe 時系列ベース判定 (将来)、量の物理妥当性チェック (別)。

## 3. 関連 docs と前提

- `AGENTS.md`「収束確認 (必須)」(`check_convergence.py`) の姉妹ルール。残差収束 ≠ 量の定常化。
- 量の抽出はメッシュ `/CELLS/centCoords` と `res_*.h5/VALUE` に依存。

## 4. 設計方針

- VERDICT 4 値: `STEADY` / `DRIFTING` (単調トレンド継続) / `OSCILLATING` (リミットサイクル→平均±振幅で報告) /
  `TRANSIENT-UNSETTLED` (スナップショット不足 or 全系列極値が末尾=過渡が減衰しきっていない)。
- 末尾 `--tail` (既定 0.4) について、トレンド `|傾き×tail幅|/平均 > --drift` (既定 0.05) で `DRIFTING`、
  振れ `(max-min)/平均 > --osc` (既定 0.10) で `OSCILLATING`、それ以外 `STEADY`。
- 組み込み量: `shock` (中心線 P 上昇位置) / `asym` (鏡像 Mach 非対称、非対称メッシュは自動 skip) /
  `machmax` / `pmax`。`--quantity` で選択。終了コード STEADY=0/それ以外=1 で gate 可能。
- メッシュは入力 h5 (`/CELLS/centCoords` を持つ非 res ファイル) を自動検出 (`--mesh` で上書き)。

## 5. 実装ステップ

1. `solver_density_cuda/tools/check_quasisteady.py` 新規 (numpy/h5py/scipy)。
2. `AGENTS.md`「準定常確認 (必須)」節を「収束確認」直後に追加。
3. メモリ `convergence-check-discipline` に「残差収束 ≠ 量の定常化」を追記。

## 6. 検証

- 本セッションの run で判定確認: Euler bp1.70 (過渡・移動中) → `DRIFTING`、SST bp1.70 (収束) → `STEADY`、
  非対称メッシュ (flat plate) で `asym` 自動 skip、終了コード STEADY=0 / DRIFTING=1。

## 7. 影響範囲

- 追加: `solver_density_cuda/tools/check_quasisteady.py`、`AGENTS.md` ルール節。既存挙動への影響なし (純追加)。

## 8. 完了条件

- [x] AGENTS.md ルール追加
- [x] ツール実装・判定検証
- [x] `status: done` / accepted 配置 / `design/README.md` 同期

## 9. 変更ログ

- `2026-06-27` — 実装・検証完了 (commit 9bc65a3)。AGENTS.md「準定常確認 (必須)」+ `check_quasisteady.py`。
  以後、派生量を報告する応答は本ツールの VERDICT を引用する。

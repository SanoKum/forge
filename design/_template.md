# <タスク名>

<!--
ファイル名規約: `<area>-<short-slug>.md` (例: `gradient-least-squares.md`)
長期テーマは `<area>-<theme>/README.md` をインデックスにしてサブディレクトリ化する。
-->

## メタ

- **area**: `<gradient | limiter | convection | diffusion | time_integration | boundary | poisson | architecture | その他>`
- **status**: `draft`  <!-- draft / in_progress / done / superseded -->
- **related_docs**:
  - `docs/<area>/theory.md`
  - `docs/<area>/implementation.md`
- **related_plans**: <!-- 親 / 子 / 関連 plan があれば列挙、無ければ削除 -->
- **created**: `YYYY-MM-DD`
- **owner**: `<担当>`

## 1. 目的

このタスクが解決する課題と、完了時に得られる状態を 2〜4 行で記述する。

## 2. スコープ

- **やる**: 本計画で実装する範囲を箇条書き。
- **やらない**: 将来計画 / 別 plan に切り出す範囲を明示。

## 3. 関連 docs と前提

理論的背景・既存実装の参照先 (`docs/<area>/theory.md` 等) と、依存する別 plan を列挙する。
理論や実装方針に変更が及ぶ場合は、本 plan の実装に入る前に対応する `docs/<area>/*.md`
を先に更新する。

## 4. 設計方針

数式・データ構造・呼び出し関係など、実装に直結する設計判断を記述する。
スキーム名や物理プロセスはここでは theory.md / implementation.md へリンクし、
本 plan には差分のみ書く。

## 5. 実装ステップ

1. <ステップ 1: 触るファイルと変更概要>
2. <ステップ 2>
3. <ステップ 3>

各ステップで触る主要ファイルを併記すること。

## 6. 検証

- **単体 / ビルド**: コンパイル・ユニット検証手順。
- **検証ケース**: `case/<case_name>/` のどれを使うか、`.github/forge-verification-cases.md` 参照。
- **判定基準**: 残差・物理量・性能の合否ライン。

## 7. 影響範囲

- 触るモジュール / ファイル一覧
- 既存ケース・実行手順への影響
- ドキュメント (`docs/<area>/*`, `docs/index.md`) の更新箇所

## 8. 完了条件

- [ ] 関連 `docs/<area>/theory.md` 更新済み
- [ ] 関連 `docs/<area>/implementation.md` 更新済み
- [ ] 実装・検証完了 (本 plan の §6 を満たす)
- [ ] `.github/plans/README.md` の状態を `done` に更新
- [ ] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `YYYY-MM-DD` — 初稿。

# Forge Agent Rules

このファイルは forge リポジトリの開発・運用ルールの正本であり、Claude (CLAUDE.md 経由) と GitHub Copilot (copilot-instructions.md 経由) の双方から参照される。

## 文書構成

各文書の役割は次の通り。詳細手順は本ファイルに重複記載せず、該当文書を参照する。

| 文書 | 役割 |
| --- | --- |
| `AGENTS.md` (本ファイル) | 全ルール・方針の正本。常時参照される基準 |
| [`.github/forge-calculation-workflow.md`](.github/forge-calculation-workflow.md) | 計算ケース準備・メッシュ生成/変換・`forge` 実行の標準手順 |
| [`.github/forge-solver-settings.md`](.github/forge-solver-settings.md) | `convMethod` / `limiter` などの数値設定リファレンス |
| [`.github/forge-development-environment.md`](.github/forge-development-environment.md) | 開発環境とビルド (Docker / WSL native) の方針 |
| [`.github/forge-coding-conventions.md`](.github/forge-coding-conventions.md) | ソース構成・C++/CUDA 命名規約・ビルド/テスト実行手順 |
| [`.github/forge-verification-cases.md`](.github/forge-verification-cases.md) | コード変更時の検証ケース選定と確認手順 |
| [`.github/plans/`](.github/plans/README.md) | 機能ごとの実装計画 (着手前に参照する基準文書) |
| [`.github/verification-cases/`](.github/verification-cases/) | 各標準検証ケースの個別手順 |
| [`docs/`](docs/index.md) | 理論 (`<area>/theory.md`) と実装解説 (`<area>/implementation.md`) |

## 言語方針

- 文書本文・コード内コメントは日本語で記述する。
- 識別子・関数名・スキーム名・ファイルパスなどコード由来の語は原語 (英語) のまま `` `code` `` 表記とする。
- git commit メッセージは英語の命令形で記述する (識別子・コード由来語は原語のまま)。

## 計算・実行ルール

計算の実行方法、ケース準備、メッシュ生成、メッシュ変換、`forge` の起動、Docker 経由の Gmsh/ParaView 利用について回答するときは、まず `.github/forge-calculation-workflow.md` を参照し、その手順に合わせて案内すること。計算手順の本文はこのファイルに重複記載しない。

`solverConfig.yaml` の `convMethod`・`limiter` などの数値設定を変更・確認するときは、必ず `.github/forge-solver-settings.md` を参照すること。設定値の意味を記憶や推測で判断しないこと。

エージェント自身が計算検証を実行する場合も、既存の `run_*` ディレクトリをそのまま使い回さず、必ず複製した新しい `run_*` ディレクトリで実行すること。
また、計算を実行した場合は `residual_history.csv` から `residual_history.png` も生成して残すこと。

開発環境に関する既定ルールは `.github/forge-development-environment.md` を参照すること。通常の開発は Docker コンテナを基本とし、NVIDIA 提供の GPU 速度計測・プロファイリングツールを使う場合は、Docker 内でうまくいかないケースを想定して WSL native または Linux native の手順を優先する。さらに、**計算速度の評価・比較・プロファイル作業は原則 native で行う** (Docker は計測値が不安定でツールも揃わないため、速度評価の基準環境としない)。詳細は `.github/forge-development-environment.md` の「速度評価・ベンチマークは native を既定とする」を参照。

コード変更時の検証ケースの選び方と個別の確認手順は `.github/forge-verification-cases.md` を参照すること。

`.github/plans/` 配下の `.md` は、実装計画ごとの基準文書として扱うこと。対象タスクに対応する計画が存在する場合、実装や調査を進める前にまずその `.md` を参照すること。計画一覧は [`.github/plans/README.md`](.github/plans/README.md) を使う。

実装方針を変更する場合は、先に対応する `.github/plans/*.md` を更新してから実装に移ること。計画未更新のまま実装だけを先行させないこと。

ユーザーが solver のコード構造、アーキテクチャ、モジュール責務について質問した場合は、`docs/architecture/overview.md` を既定の参照先とすること。

## ドキュメント整備方針

forge の理論的背景と実装解説は `docs/` 配下に機能単位 (物理プロセス / 数値処理) のサブディレクトリで蓄積する。運用方針の本文は `docs/README.md` を参照し、本ファイルには重複記載しない。

- 新規に理論背景や実装解説を書くときは `docs/<area>/` (例: `gradient/`, `limiter/`, `convection/`, `diffusion/`, `time_integration/`, `boundary/`, `poisson/`, `architecture/`) の配下に追加する。スキーム名 (Roe / SLAU / KEEP など) は対応する機能ディレクトリ内のセクションまたはファイルとして扱い、新たなトピック単位ディレクトリは作らない。
- 既存ドキュメントがある話題は新規ファイルを作らず、該当ファイルを更新する。
- 理論・実装に関する質問に答えるときは、該当する `docs/<area>/*.md` を既定の参照先とする。
- ファイルの追加・改名・削除を行った場合は `docs/index.md` の目次を必ず同期させる。
- 数式は KaTeX 記法 (インライン `$...$`、ブロック `$$...$$`)、図は mermaid フェンスまたは画像で統一する。
- `docs/` は恒常的な理論・実装解説、`.github/forge-*.md` / `.github/plans/` / `.github/verification-cases/` はそれぞれ運用手順・実装計画・検証ケース手順という棲み分けを守る。

## 開発フロー (新規機能・設計変更)

バグ修正と軽微な変更を除く「新規機能追加」「数値スキームや設計方針の変更」を行うときは、必ず次の 4 ステップを順に踏むこと。

1. 該当する [`docs/<area>/theory.md`](docs/) を更新する (理論的な変更がある場合)。
2. 該当する [`docs/<area>/implementation.md`](docs/) を更新する (実装方針・ソース対応を反映)。
3. [`.github/plans/_template.md`](.github/plans/_template.md) を雛型に `.github/plans/<area>-<short-slug>.md` を作成し、`related_docs` で 1./2. のファイルをリンクする。[`.github/plans/README.md`](.github/plans/README.md) の一覧にも追記する。
4. 以上 3 つが揃ってから実装に着手する。

完了時には次を行うこと。

- plan の `status` を `done` に更新し、`## 変更ログ` に実装・検証結果を追記する。
- [`.github/plans/README.md`](.github/plans/README.md) の一覧表も同期させる。
- 1./2. で触った docs と [`docs/index.md`](docs/index.md) の整合性を確認する。

例外 (本フローを要しないもの):

- typo やコメント修正のみ
- 1 ファイル内で閉じるバグ修正 (振る舞い変更なし)
- リネーム/象徴のない単純リファクタ (振る舞い同一)
- docs のみの記述補正

判断がつかないときはフローを踏む側を選ぶこと。

## コミット・push 運用

実装または修正が一段落し、検証 (ビルド成功 + 該当検証ケースの結果が妥当) まで確認できたら、その機能・計画単位で git commit し、現在の feature ブランチへ push すること。本節は検証済みマイルストーンでの commit/push を事前承認するものとして扱い、都度の確認は不要とする。

- **commit タイミング**: 機能・計画単位でまとめて commit する。段階検証ごと・小さな編集ごとには commit しない。plan を伴う作業では plan の `status: done` 化と整合させる。
- **commit 対象**: `solver_density_cuda/` などのソース、`docs/`、`.github/plans/` といった意味あるソース・文書変更を対象とする。`case/` 配下の検証 run 成果物 (`res_*.h5`, `*.xmf`, `residual_history.csv` / `.png`, `forge_run.log` など) は commit に含めない。検証 run の入力 config (`solverConfig.yaml`, `bcondConfig.yaml` 等) は、リファレンスとして意図的に残したい場合のみ明示的に `git add` する。
- **push**: commit 後は現在の feature ブランチへ push する。`main` へ直接 commit / push しない。`main` 上にいる場合は先に feature ブランチを切ってから作業する。
- **commit メッセージ**: 英語の命令形で記述する (例: `Add Menter SST turbulence model with dilatation correction`)。
- 無関係な既存の作業ツリー変更を巻き込まないよう、commit 前に `git status` / `git diff` で対象を確認すること。

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
| [`.github/forge-divergence-and-startup.md`](.github/forge-divergence-and-startup.md) | 発散の主因と安定起動の手順 (新規計算 / 発散時に必ず参照) |
| [`.github/forge-verification-cases.md`](.github/forge-verification-cases.md) | コード変更時の検証ケース選定と確認手順 |
| [`.github/forge-su2-cross-check.md`](.github/forge-su2-cross-check.md) | 同一メッシュ・同一 BC で SU2 と比較し forge 固有の問題を切り分ける手順 |
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

**新規計算の投入時、および計算が発散 (NaN / 残差爆発) したときは、まず [`.github/forge-divergence-and-startup.md`](.github/forge-divergence-and-startup.md) を参照し、そのチェックリストと段階起動手順に従うこと。** 発散の大半は物理やメッシュでなく投入設定 (初期値が入口流れと不整合 / 超音速向け BC を亜音速に使用・出口逆流 Pt/Tt 未設定 / 初手から 2 次移流・乱流・no-slip / 非直交での free-stream 桁落ち) が原因であり、易しい条件で収束させてから引き継ぎ計算で段階的に上げるのが基本。

エージェント自身が計算検証を実行する場合も、既存の `run_*` ディレクトリをそのまま使い回さず、必ず複製した新しい `run_*` ディレクトリで実行すること。
また、計算を実行した場合は `residual_history.csv` から `residual_history.png` も生成して残すこと。

**計算結果ディレクトリの明示 (必須)**: 後からどの検証がどこにあるか追えるよう、計算の「投入時」と「結果まとめ時」の両方で、結果が出力される `run_*` ディレクトリを**リポジトリルートからの相対パスで明示**すること。

- **投入時**: 計算を起動する応答内で、各 run の出力先パスを列挙する (例: `case/05.sod_shock_tube/run_0004_slau_tracer_n2n2/`)。複数 run を投入する場合は run ごとに何を変えたか (スキーム・設定差分) を 1 行で添える。
- **結果まとめ時**: 結論や比較表とあわせて、その結論の根拠となった `run_*` パスと主要成果物 (`res_*.h5` / `residual_history.png` / 後処理図 `*.png` 等) を明示する。複数 run を比較したときは「どの図/数値がどの run のものか」を取り違えないよう run パスを併記する。
- 一時 run や破棄予定の run を作った場合も、その旨とパスを明記する (放置して所在不明にしない)。

**case の run 索引 (必須・恒久)**: 応答での明示はスクロールで流れて後から辿れないため、**各 case ディレクトリの `README.md` に「## 計算 run 一覧」表を維持する**。これが「どの検証がどの run か」の恒久的な一次情報であり、結果をまとめる応答ではこの表へのポインタも示す。

- **対象 case に `README.md` が無ければ新規作成**し、ケース概要に続けて run 一覧表を置く。
- **run_* を作成・破棄したら同表を必ず同期**する。列は最低限: `run_*` ディレクトリ名 / 目的・主要設定差分 / 主要結果・成果物 / 状態 (`active` / `破棄予定` / `ref`=入力リファレンスとして保持)。
- **命名規約**: `run_NNNN_<slug>` とし、`<slug>` は目的が一目で分かる語にする (例 `run_0015_prod_double`)。**既存 run と番号が衝突しない連番**を使い、衝突や系列違いが起きたら表の備考でグルーピングを明記する。
- 表の本文は簡潔に (1 run = 1 行)。詳細な考察は `.github/plans/` の該当 plan に書き、表からはその plan へリンクする。

**NaN / 発散チェック (必須)**: 計算が「完走 (exit 0)」しても妥当とは限らない。**初期ステップから発散している**ことが多いため、結果を使う前に必ず NaN/Inf の有無を確認すること。

- **確認タイミング**: ① 投入直後 (序盤数十ステップで NaN/Inf が出ていないか早期確認)、② 結果まとめ前 (最終 `res_*.h5` と全 `residual_history.csv` を確認)。長時間 run では途中でも定期的に見る。
- **確認内容**:
  - `residual_history.csv` の `rms_*` 列に NaN/Inf が無いか、**最初に NaN になったステップ** (`step 0`/`step 1` からなら初手発散) を特定する。
  - 最終および中間 `res_*.h5` の `VALUE/*` (特に `ro`,`P`,`T`,`roe`,勾配 `dUxdx` 等) に NaN/Inf が無いか、値が物理的か (例: 静圧 ≤ 全圧、`ro>0`、`T>0`) を確認する。
- **発散していた場合**: 「収束した」と報告しない。最初に NaN が出たステップ・場所 (どの境界・どのセル) を切り分け、原因 (BC のゼロ割・IC とBC の不整合・CFL 過大・メッシュ不良 等) を特定してから対処する。`max cfl` のログだけ見て安定と判断しないこと (CFL が有限でも残差が NaN のことがある)。
- 早期確認の例: `python3` で `residual_history.csv` の先頭数十行の `rms_ro` と `res_0/res_1` 相当出力の `isnan` を見る。専用ツールがあれば優先する。

**収束確認 (必須・NaN チェックとは別)**: 「NaN が無い」=「収束した」ではない。結果を使う/報告する前に、各保存量の残差が**実際に下がりきっているか**を必ず確認すること。

- **判定は手作業の目視ではなく `solver_density_cuda/tools/check_convergence.py <run_dir> ...` を実行して行う** (本ルールの実体化ツール)。全残差列の低下桁数・トレンド・NaN を判定し `PASS / NOT CONVERGED / DIVERGED` を返す。**「収束した」「一致した」と報告する応答には、このツールの VERDICT を根拠として必ず貼ること** (これを怠り `rms_ro` だけで「収束」と誤報告した事例があるため、ツール経由を必須とする)。未収束なら「未収束」と明記し、未収束のトランジェント同士の場の比較を「一致」と表現しない。
- **`rms_ro` だけで判断しない**。`residual_history.csv` の **全列** (`rms_ro`,`rms_roUx`,`rms_roUy`,`rms_roUz`,`rms_roe`、RANS 時は `rms_roK`,`rms_roOmega`) のトレンドを見る。`rms_ro` が低くても、運動量 (特に軸対称の `rms_roUy`) や乱流 (`rms_roK`/`rms_roOmega`) が**下げ止まり・横ばい・上昇**していれば未収束。実例: 軸対称 SST で `rms_ro`≈3e-5 でも `rms_roUy`≈1e-2 停滞・`rms_roK` 増大=近軸が未収束だった ([architecture-axisym-axis-singularity.md](.github/plans/architecture-axisym-axis-singularity.md))。
- **残差プラトーは「収束」ではない**。下げ止まる場合は、積分量 (massflux/推力/出口諸量) が**定常化**しているか、場が**発達しきっている**かを併せて確認する (リミットサイクルの可能性)。
- **場の発達も確認する** ([develop-flow-before-reporting] と同趣旨): 残差が下がっていても、境界層・乱流・衝撃などが発達途中なら結果は使えない。中間 `res_*.h5` を時系列で見て、注目量が定常化したことを確認する。
- 外部ソルバ (SU2 等) でクロスチェックする場合の収束確認も同様。手順は [`.github/forge-su2-cross-check.md`](.github/forge-su2-cross-check.md) を参照。

forge の結果が「軸対称・乱流・近軸で forge だけ妙な値になる」ようなときは、推測で結論づけず [`.github/forge-su2-cross-check.md`](.github/forge-su2-cross-check.md) の手順で **同一メッシュ・同一 BC の SU2 と比較**して切り分けること。

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

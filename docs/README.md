# docs/

このディレクトリは forge の**現在の仕様と解説** (理論的背景 + 実装解説) を蓄積する場所である。
「読めば現在の挙動が分かる」ことを目的とし、設計判断の経緯 (なぜそうしたか・これから何を変えるか)
は `design/`、調査メモは `notes/` に置く。個別の検証手順・計算ワークフローは `.github/` 配下の
メタ文書を参照すること (役割分担は後述)。

## 構成方針

機能単位 (物理プロセス / 数値処理) でサブディレクトリを切る。スキーム名 (Roe / SLAU /
KEEP など) は機能ディレクトリ内のセクションまたはファイルとして扱い、トピック単位の
ディレクトリは増やさない。

| ディレクトリ | 対象 |
| --- | --- |
| `architecture/` | ソルバ全体構成・モジュール責務・データフロー |
| `gradient/` | 勾配再構成 (Green-Gauss) |
| `limiter/` | スカラ・リミッタ (Barth-Jespersen, Venkatakrishnan, Nishikawa) |
| `convection/` | 対流フラックス (Roe, SLAU, KEEP, AUSM, AUSM+UP, FVS, HLLE など) |
| `diffusion/` | 粘性・熱伝導フラックス |
| `time_integration/` | 陽解法・陰解法 (RK, defect correction, block DPLUR, CFL 制御) |
| `boundary/` | 境界条件 (壁面, 流入流出, 対称, 周期) |
| `poisson/` | Poisson 求解 (壁面距離, AMGCL/AMGX) |

新しい領域 (例: 乱流モデル) が必要になった時点でディレクトリを新設し、
[`index.md`](index.md) と [`.github/copilot-instructions.md`](../.github/copilot-instructions.md)
の例示にも追記する。

## ファイル命名

各機能ディレクトリ配下は、その領域の**現在仕様ドキュメント**を 1 つ以上置く。
理論と実装の**強制的な 2 ファイル分割 (`theory.md` / `implementation.md`) は廃止**し、
次の指針で構成する。

- 小〜中規模の領域: 1 ファイル内に節で分ける (例: `## 目的` / `## 理論` / `## forge での実装`
  / `## 設定項目` / `## 既知の落とし穴` / `## 検証ケース` / `## 関連 design`)。
- 大規模な領域: 理論と実装を分けたい場合のみ `theory.md` / `implementation.md` に分割してよい
  (既存の分割ファイルはそのまま維持し、肥大化したものから統合していく)。
- スキーム比較が肥大化した場合は `convection/schemes/<name>.md` のようにファイル単位で分割する。

## 記法

- 数式は KaTeX 記法 (インライン `$...$`、ブロック `$$...$$`)。
- 図は mermaid フェンス (` ```mermaid `) または画像ファイル。
- ソース参照は相対リンクで記載する。

## 既存メタ文書との棲み分け

| 場所 | 役割 |
| --- | --- |
| `docs/` (本ディレクトリ) | 現在の仕様と解説 (恒常的) |
| `design/` | 変更単位の設計判断 (`active/` 進行中 / `accepted/` 現役 / `archived/` 終了) |
| `notes/` | 調査メモ・作業ログ (`investigations/` / `sessions/`) |
| `.github/forge-calculation-workflow.md` | 計算実行・メッシュ生成手順 |
| `.github/forge-development-environment.md` | 開発環境・ビルド・プロファイル |
| `.github/forge-verification-cases.md`, `.github/verification-cases/` | 検証ケース選定と確認手順 |

## 目次の維持

ファイルを追加・改名・削除した場合は [`index.md`](index.md) を必ず同期する。

## `design/` との関係

`docs/` は**現在の仕様**を蓄積する場所で、[`design/`](../design/README.md) は
それを変えるための**設計判断 (変更単位)** を置く場所である。
新規機能や設計変更を行うときは、先に該当する `docs/<area>/` の現在仕様を更新し、
その上で [`design/_template.md`](../design/_template.md) を雛型に `design/active/` へ計画を
作って実装に着手する。完了したら計画を `design/accepted/` (または `archived/`) へ移す。
フローの詳細は [`../AGENTS.md`](../AGENTS.md) の「## 開発フロー」節を参照。

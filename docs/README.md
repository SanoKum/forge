# docs/

このディレクトリは forge の理論的背景と実装解説を蓄積する場所である。
個別の検証手順・計算ワークフロー・実装計画は `.github/` 配下のメタ文書を参照すること
(役割分担は後述)。

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

各機能ディレクトリ配下では原則として次の 2 ファイルを使う。両方必須ではない。

- `theory.md` — 支配方程式・離散化・スキームの理論的背景。
- `implementation.md` — 該当ソースファイル・関数とデータフロー、設定項目。

スキーム比較が肥大化した場合は `convection/schemes/<name>.md` のように
ファイル単位で分割する余地を残す。

## 記法

- 数式は KaTeX 記法 (インライン `$...$`、ブロック `$$...$$`)。
- 図は mermaid フェンス (` ```mermaid `) または画像ファイル。
- ソース参照は相対リンクで記載する。

## 既存メタ文書との棲み分け

| 場所 | 役割 |
| --- | --- |
| `docs/` (本ディレクトリ) | 理論背景・実装解説 (恒常的) |
| `.github/forge-calculation-workflow.md` | 計算実行・メッシュ生成手順 |
| `.github/forge-development-environment.md` | 開発環境・ビルド・プロファイル |
| `.github/forge-verification-cases.md`, `.github/verification-cases/` | 検証ケース選定と確認手順 |
| `.github/plans/` | 実装計画 (タスク単位) |

## 目次の維持

ファイルを追加・改名・削除した場合は [`index.md`](index.md) を必ず同期する。

## `.github/plans/` との関係

`docs/` は恒常的な理論・実装解説を蓄積する場所で、`.github/plans/` は
それを変えるための一時的な実装計画 (タスク単位) を置く場所である。
新規機能や設計変更を行うときは、先に該当する `docs/<area>/theory.md` /
`implementation.md` を更新し、その上で `.github/plans/_template.md` を雛型に
plan を作って実装に着手する。フローの詳細は
[`../.github/copilot-instructions.md`](../.github/copilot-instructions.md) の
「開発フロー」節を参照。

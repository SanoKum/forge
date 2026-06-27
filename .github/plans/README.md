# （移転）`.github/plans/` → `design/` + `notes/`

旧 `.github/plans/` 配下の実装計画・調査メモは、文書の役割ごとに次へ移転した。

| 旧 | 新 | 役割 |
| --- | --- | --- |
| 実装計画 (検討中・進行中) | [`design/active/`](../../design/active/) | いま手を動かしている計画 |
| 実装計画 (実装・検証済み) | [`design/accepted/`](../../design/accepted/) | 現役の設計判断 (`done`) |
| 実装計画 (superseded / 終了) | [`design/archived/`](../../design/archived/) | 役目を終えた計画 |
| 技術調査・サーベイ | [`notes/investigations/`](../../notes/investigations/) | コード変更を伴わない研究 |
| 作業セッション用プロンプト | [`notes/sessions/`](../../notes/sessions/) | 使い捨て・恒久参照しない |

- 計画インデックス: [`design/README.md`](../../design/README.md)
- 計画テンプレート: [`design/_template.md`](../../design/_template.md)
- 運用フロー: [`../../AGENTS.md`](../../AGENTS.md) の「## ドキュメント体系」「## 開発フロー」節

**新規の計画・調査は本ディレクトリに追加しない。** 上記の移転先に置くこと。

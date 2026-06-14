# forge ドキュメント目次

forge の理論的背景と実装解説の索引。運用方針は [README.md](README.md) を参照。

新規ファイル追加・改名・削除時は本目次を必ず更新すること。

## 機能ディレクトリ

| 領域 | theory | implementation | 状態 |
| --- | --- | --- | --- |
| アーキテクチャ全体 | — | [architecture/overview.md](architecture/overview.md) | 整備済み |
| 離散化レイアウト (cell/node) | [discretization/theory.md](discretization/theory.md) | [discretization/implementation.md](discretization/implementation.md) | 整備中 (median-dual 両対応) |
| 勾配再構成 | [gradient/theory.md](gradient/theory.md) | [gradient/implementation.md](gradient/implementation.md) | 整備済み |
| リミッタ | [limiter/theory.md](limiter/theory.md) | [limiter/implementation.md](limiter/implementation.md) | 整備済み |
| 対流フラックス | [convection/theory.md](convection/theory.md) | [convection/implementation.md](convection/implementation.md) | 整備済み |
| 粘性・熱伝導 | [diffusion/theory.md](diffusion/theory.md) | [diffusion/implementation.md](diffusion/implementation.md) | 整備済み |
| 時間積分 | [time_integration/theory.md](time_integration/theory.md) | [time_integration/implementation.md](time_integration/implementation.md) | 整備済み |
| 境界条件 | [boundary/theory.md](boundary/theory.md) | [boundary/implementation.md](boundary/implementation.md) | 整備済み |
| Poisson 求解 | [poisson/theory.md](poisson/theory.md) | [poisson/implementation.md](poisson/implementation.md) | 整備済み |
| 軸対称 | [axisymmetric/theory.md](axisymmetric/theory.md) | [axisymmetric/implementation.md](axisymmetric/implementation.md) | 整備済み |
| 乱流モデル | [turbulence/theory.md](turbulence/theory.md) | [turbulence/implementation.md](turbulence/implementation.md) | 整備済み |
| 多成分熱物性 (TP gas) | [thermophysics/theory.md](thermophysics/theory.md) | [thermophysics/implementation.md](thermophysics/implementation.md) | M1-M4 + TP境界条件 完了 |
| 非平衡凝縮 (4 モーメント) | [condensation/theory.md](condensation/theory.md) | [condensation/implementation.md](condensation/implementation.md) | 整備中 (Phase 1 受動スカラー骨格) |

新規領域を追加した場合は本表に行を追加すること。

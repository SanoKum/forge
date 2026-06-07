# forge 実装計画インデックス

`.github/plans/` 配下の実装計画 (`.md`) の一覧と状態。
新規 plan を追加・状態を変えた場合は、本表を必ず同期させる。

計画の作り方は [`_template.md`](_template.md) を起点に、命名規約 `<area>-<short-slug>.md`
(長期テーマは `<area>-<theme>/README.md`) で配置する。
運用フロー全体は [`../../AGENTS.md`](../../AGENTS.md) の
「## 開発フロー」節を参照。

## 一覧

| Plan | area | status | related_docs | 概要 |
| --- | --- | --- | --- | --- |
| [diffusion-viscous-shear-flux.md](diffusion-viscous-shear-flux.md) | diffusion | done | [`docs/diffusion/`](../../docs/diffusion/) | 粘性せん断フラックスの修正 (内部面の法線項 `delta_x→delta` + 転置項 + 壁面 no-slip 抗力、軸平行格子で横方向粘性が落ちる不具合、SU2 検証済) |
| [gpu-implicit-plan.md](gpu-implicit-plan.md) | time_integration | in_progress | [`docs/time_integration/`](../../docs/time_integration/) | GPU 上での pseudo-time implicit (matrix-free + 近似 Jacobian + block-Jacobi/defect-correction) 段階導入 |
| [architecture-axisymmetric.md](architecture-axisymmetric.md) | architecture | done | [`docs/axisymmetric/`](../../docs/axisymmetric/) | 軸対称モード Phase 1 (B 流儀 r 重み付け + 圧力ソース + 軸 BC、非粘性) |
| [architecture-axisym-nozzle-geometry.md](architecture-axisym-nozzle-geometry.md) | architecture | done | [`docs/axisymmetric/`](../../docs/axisymmetric/) | 軸対称検証ノズルの喉部直後を 5 次多項式で滑らかに接続する geometry 改善 |
| [architecture-rans-sst.md](architecture-rans-sst.md) | architecture | done | [`docs/turbulence/`](../../docs/turbulence/) | Menter SST (k-ω) を explicit 軸対称ノズル (run_0087〜0090) で 4 段階検証完了。advection・diffusion・source (F1/F2 ブレンド)・渦粘性すべて実装済。軸対称 geometric source は子 plan へ |
| [architecture-axisym-sst.md](architecture-axisym-sst.md) | architecture | done | [`docs/turbulence/`](../../docs/turbulence/) | 親 SST の子 plan。フープひずみ $2(u_r/r)^2$ (run_0091) + 生産項の圧縮性補正 `dilatationCorrection` (0:off/1:deviatoric/2:+等方項) を実装・段階検証。run_0093(A+B) で k −12%・vis_turb −14% の膨張減衰を確認。既定値 2 (全 SST ケースに適用) |

## 状態の意味

- `draft` — 設計検討中、実装未着手。
- `in_progress` — 実装または検証進行中。
- `done` — 関連 docs 更新・実装・検証がすべて完了。
- `superseded` — 別 plan に置き換えられた。後継 plan を `related_plans` に明記すること。

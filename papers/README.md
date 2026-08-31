# papers

参照文献 (一次資料) の置き場。理論・スキームの裏取りに使う **PDF 本体**をここに置く。

## 運用方針

- **調査内容ごとにサブディレクトリ `papers/<topic>/` へまとめる**。例: `nozzle_design/` (下位に `dual_bell/`), `condensation/`, `pintle_nozzle/`, `les_des/`, `multicomponent/`, `shock_control/`, `precision/`。
- **PDF 本体は `papers/` に置き、git では追跡しない**(大容量バイナリのため。`<topic>/` 配下に置く)。
- **調査の「書きもの」(サーベイ・文献まとめ・論文要約・レビュー成果物 `*.pptx` 等) は `papers/` ではなく [`notes/investigations/`](../notes/README.md) に置く**。そこから本ディレクトリの PDF を相対リンク (`../../papers/<topic>/...`) で参照する。
  - 例: [`notes/investigations/nozzle-pintle-literature-review.md`](../notes/investigations/nozzle-pintle-literature-review.md) ↔ [`pintle_nozzle/`](pintle_nozzle/)、[`notes/investigations/condensation-droplet-hypersonic-survey.md`](../notes/investigations/condensation-droplet-hypersonic-survey.md) / [`condensation-nitrogen-nozzle-flows-summary.md`](../notes/investigations/condensation-nitrogen-nozzle-flows-summary.md) ↔ [`condensation/`](condensation/)、[`notes/investigations/nozzle-throat-curvature-shape-representation-survey.md`](../notes/investigations/nozzle-throat-curvature-shape-representation-survey.md) ↔ [`nozzle_design/`](nozzle_design/)。
- ビルド成果物の生成スクリプトや図切り出し (`build_pptx.py`, `.figs/`, 生ログ `*_raw.json` 等) は、対応する `papers/<topic>/` に置く (PDF と同じ場所で完結させる)。

> ルートに散在している既存 PDF は順次 `<topic>/` へ寄せてよい。新規追加は最初から `papers/<topic>/` に置くこと。

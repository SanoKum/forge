# forge 調査・作業ノート (`notes/`)

`notes/` は恒久仕様 (`methods/`) でも変更単位の設計判断 (`plans/`) でもない、**調査メモ・作業ログ・一時的知見**の逃がし先。
ここに置いたものは「現在地」を示す `plans/` の一覧を汚さない。

- `investigations/` — 技術調査・サーベイ・外部ソルバ棚卸し (コード変更を伴わない研究)。結論は該当 `plans/` 計画や `methods/` に反映する。
- `sessions/` — 使い捨ての作業セッション用プロンプト・引き継ぎメモ。恒久参照を意図しない。

## investigations

| ノート | area | 概要 |
| --- | --- | --- |
| [condensation-droplet-hypersonic-survey.md](investigations/condensation-droplet-hypersonic-survey.md) | `condensation / multiphase` | 極超音速燃焼風洞・飛翔体における凝縮液滴/雨滴の挙動 (風洞凝縮の蒸発vs成長・燃焼器到達 end-to-end・迎角依存) |
| [node-slip-tangential-density-spurious-flow.md](investigations/node-slip-tangential-density-spurious-flow.md) | `boundary / discretization` | node slip 境界 + 接線密度勾配で市松状スプリアス接線流が定在する未修正バグの切り分け記録 (case/24 純伝導検証の副産物) |
| [condensation-nitrogen-nozzle-flows-summary.md](investigations/condensation-nitrogen-nozzle-flows-summary.md) | `condensation` | "On Nitrogen Condensation in Hypersonic Nozzle Flows" 論文精読まとめ + Arthur 1952 形状の forge dry 再現 (PDF: [papers/condensation/](../papers/condensation/)) |
| [nozzle-pintle-literature-review.md](investigations/nozzle-pintle-literature-review.md) | `nozzle` | ピントル/ニードル弁式 超音速ノズル推力調整機構 — 文献・特許まとめ (PDF/成果物: [papers/pintle_nozzle/](../papers/pintle_nozzle/)) |
| [convection-central-scheme-oscillation-control.md](investigations/convection-central-scheme-oscillation-control.md) | `convection` | LES/DES における中心差分 (KEEP) のスプリアス振動抑制 技術調査 (振動源4分類の統括) |
| [convection-pep-scheme-survey.md](investigations/convection-pep-scheme-survey.md) | `convection` | PEP (Pressure-Equilibrium-Preserving) 系スキーム技術調査と forge 実装方針 |
| [nozzle-optimization-tool-survey.md](investigations/nozzle-optimization-tool-survey.md) | `—` | 超音速・極超音速ノズル最適化ツール — 技術動向調査と開発フロー提案 |
| [su2-nemo-contact-thermo-investigation.md](investigations/su2-nemo-contact-thermo-investigation.md) | `convection / thermophysics` | SU2-NEMO contact/interface thermo 取り扱い調査 (forge mixed-order face-state 比較) |
| [turbulence-des-flux-survey.md](investigations/turbulence-des-flux-survey.md) | `turbulence` | DES/DDES/IDDES 用 低散逸対流 flux 設計 技術調査（圧縮性 LES/DES） |
| [turbulence-des-wmles-survey.md](investigations/turbulence-des-wmles-survey.md) | `turbulence` | DES / WMLES 最新動向サーベイ（超音速ノズル・ピントルバルブ適用向け） |

## sessions

| ノート | 概要 |
| --- | --- |
| [condensation-nonequilibrium-session-prompt.md](sessions/condensation-nonequilibrium-session-prompt.md) | 次セッション用プロンプト: 非平衡凝縮 (4 モーメント方程式) の forge 実装 |
| [cutler-tp-multispecies-cfl-analysis-prompt.md](sessions/cutler-tp-multispecies-cfl-analysis-prompt.md) | 分析依頼プロンプト: 多成分 TP の陰解法 (block-DPLUR) で安定 CFL が ~1 に頭打ちする原因 |
| [implicit-acceleration-session-prompt.md](sessions/implicit-acceleration-session-prompt.md) | 別セッション用プロンプト: 陰解法の安定 CFL 引き上げ |
| [precision-mixed-axisym-session-prompt.md](sessions/precision-mixed-axisym-session-prompt.md) | 新セッション引き継ぎプロンプト — 軸対称 近軸の陰解法を混合精度で root-fix |
| [species-in-dplur-session-prompt.md](sessions/species-in-dplur-session-prompt.md) | 引継ぎプロンプト: 化学種 `roY_s` を block-DPLUR に結合する (roe↔roY coupling) |

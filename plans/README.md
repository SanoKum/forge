# forge 計画・設計判断インデックス (`plans/`)

`plans/` は forge の**変更単位 (変更計画・設計判断・検証ログ)** の文書を置く。`methods/` が「現在の仕様・解説」なのに対し、
`plans/` は「**なぜその設計にしたか / これから何を変えるか**」を担う。

計画の作り方は [`_template.md`](_template.md) を起点に、命名規約 `<area>-<short-slug>.md` で配置する。
運用フロー全体は [`../AGENTS.md`](../AGENTS.md) の「## 開発フロー」節を参照。

> **凡例 (旧パス)**: これらの計画は以前 GitHub 配下の `.github/plans/` に置かれていた (現在の root `plans/` とは別物)。
> 各計画の完了ログや起票履歴に出てくる `.github/plans/<name>.md` という記述は、**現在の root
> `plans/{active,accepted,archived}/<name>.md`** を指す (移行済み)。履歴の文言はそのまま残してあるので、
> 当時の記録として読むこと。

ライフサイクルは**フォルダ移動**で表す (ファイル名は固定し、移動でリンクを壊さない)。

- `active/` — 検討中・実装/検証中 (`draft` / `in_progress`)。**いま手を動かしている計画はここだけ**。
- `accepted/` — 実装・検証済みで、**いまも現在挙動を支配している設計判断** (`done`)。仕様は `methods/` 側、ここには「なぜ」を残す。
- `archived/` — superseded / 棄却 / 役目を終えた計画 (`superseded`)。後継を本文に明記する。

各計画の詳細 (背景・検証 run・結論) は当該 `.md` 本文を参照。本インデックスは航行用の一覧に徹する。

## active (検討中・進行中)

| Plan | area | 概要 |
| --- | --- | --- |
| [architecture-node-centroid-value-position.md](active/architecture-node-centroid-value-position.md) | `architecture` | node-centered の centCoords を「値の位置 (ノード座標)」に統一し、双対重心/軸半径を分離する |
| [architecture-perphase-profiling-hotspot.md](active/architecture-perphase-profiling-hotspot.md) | `architecture` | per-phase 詳細プロファイリングに基づく host 律速の局所最適化 |
| [condensation-nonequilibrium.md](active/condensation-nonequilibrium.md) | `condensation` | 非平衡凝縮 (4 モーメント方程式) の forge 実装 |
| [convection-freestream-preserving-flux.md](active/convection-freestream-preserving-flux.md) | `convection` | 対流流束の free-stream 保存 (基準静圧差分) |
| [convection-multispecies-contact-pressure.md](active/convection-multispecies-contact-pressure.md) | `convection` | 多成分 TP 接触面の圧力振動 (face-state 整合化 → PEP 切り分け) |
| [convection-keep-revive-node.md](active/convection-keep-revive-node.md) | `convection / turbulence` | KEEP スキーム復活 (modern API・cell/node) + node WALE で LES/ILES |
| [discretization-lsq-gradient.md](active/discretization-lsq-gradient.md) | `gradient / architecture` | node-centered 最小二乗 (LSQ) 勾配 |
| [discretization-median-dual.md](active/discretization-median-dual.md) | `architecture` | cell-centered / node-centered (median-dual) 両対応化 |
| [discretization-median-dual-3d.md](active/discretization-median-dual-3d.md) | `architecture` | 3D median-dual (M4): 3D 双対生成 + periodic 双対面 (3D node DDES/LES の起点) |
| [discretization-node-boundary-ghostless.md](active/discretization-node-boundary-ghostless.md) | `boundary` | node-centered 境界のゴースト撤廃 (段階導入: まず壁, 次に流入出/slip) |
| [gpu-implicit-plan.md](active/gpu-implicit-plan.md) | `time_integration` | GPU 陰解法化計画 |
| [turbulence-iddes-sst.md](active/turbulence-iddes-sst.md) | `turbulence` | SST-DDES / SST-IDDES 実装計画 |
| [verification-passive-pseudoshock-control.md](active/verification-passive-pseudoshock-control.md) | `verification` | 多孔壁パッシブコントロールによる擬似衝撃波抑制の逆解析 (case/36) |

## accepted (現役の設計判断)

| Plan | area | 概要 |
| --- | --- | --- |
| [boundary-cell-periodic-conservation.md](accepted/boundary-cell-periodic-conservation.md) | `boundary` | cell 全周期境界の非保存バグ修正: device `bint_d[partnerCellID]` 未転送で ghost が誤値→運動量注入。setPeriodicPartner 直後に H2D コピーで根治。TGV で cell/node 一致を確認 |
| [thermophysics-eos-positivity-floor-config.md](accepted/thermophysics-eos-positivity-floor-config.md) | `thermophysics` | EOS 正値化フロア (pMin/roMin/tMin) の config 化。無次元・低圧ケース (Taylor-Green) で既定 1.0 Pa フロアが場を破壊する問題を解消、既定値据え置きでビット不変 |
| [boundary-inlet-profile.md](accepted/boundary-inlet-profile.md) | `boundary` | 入口分布プロファイル: CSV テーブルで inlet bvar を非一様化 (x/y/z 1D 線形 or xyz 最近傍)、壁法則 helper |
| [turbulence-node-wall-function-coverage.md](accepted/turbulence-node-wall-function-coverage.md) | `boundary` | node SST 壁関数の生産置換を第一内層ノードにも適用 (近壁 k 暴走修正、cell 不変・x_R が SU2 整合) |
| [architecture-axisym-axis-singularity.md](accepted/architecture-axisym-axis-singularity.md) | `architecture` | 軸対称 近軸の数値問題 (軸中心 k スパイク) の根本原因特定 |
| [architecture-axisym-nozzle-geometry.md](accepted/architecture-axisym-nozzle-geometry.md) | `architecture` | architecture-axisym-nozzle-geometry |
| [architecture-axisym-sst.md](accepted/architecture-axisym-sst.md) | `architecture` | architecture-axisym-sst — 軸対称 SST 幾何項 子計画 |
| [architecture-axisymmetric.md](accepted/architecture-axisymmetric.md) | `architecture` | architecture-axisymmetric — 軸対称ソルバ Phase 1 |
| [architecture-detect-nan.md](accepted/architecture-detect-nan.md) | `architecture` | detectNaN — NaN/Inf 検知・自動ダンプ診断オプション |
| [architecture-rans-sst.md](accepted/architecture-rans-sst.md) | `architecture` | architecture-rans-sst — RANS SST 親計画 |
| [architecture-residual-monitor-async.md](accepted/architecture-residual-monitor-async.md) | `architecture` | 残差モニタの device 常駐化 (per-step host 同期の除去) |
| [boundary-outlet-characteristic-backflow.md](accepted/boundary-outlet-characteristic-backflow.md) | `boundary` | 静圧固定流出 BC の特性ベース化と逆流統一 |
| [convection-slau2-lowmach.md](accepted/convection-slau2-lowmach.md) | `convection` | SLAU2 圧力流束 (低マッハ圧力–速度カップリング是正) |
| [diffusion-node-boundary-real-distance.md](accepted/diffusion-node-boundary-real-distance.md) | `diffusion` | node 境界半割面 拡散の実距離 over-relax 化 (∇·S 弱形式の置換) |
| [diffusion-node-wall-viscous-distance.md](accepted/diffusion-node-wall-viscous-distance.md) | `diffusion` | node-centered 粘性壁 flux の法線距離修正 (mirror-ghost 退化対策) |
| [diffusion-turbulent-thermal-conductivity.md](accepted/diffusion-turbulent-thermal-conductivity.md) | `diffusion` | 乱流熱伝導 (turbulent thermal conductivity) の追加 |
| [diffusion-viscous-shear-flux.md](accepted/diffusion-viscous-shear-flux.md) | `diffusion` | 粘性せん断フラックスの修正 (内部面の法線項 + 壁面 no-slip 抗力) |
| [discretization-node-wall-implicit-dirichlet.md](accepted/discretization-node-wall-implicit-dirichlet.md) | `boundary / time_integration` | node-centered 壁 no-slip の implicit Jacobian Dirichlet 化 (SU2 DeleteValsRowi 相当) |
| [output-node-wall-surface-viz.md](accepted/output-node-wall-surface-viz.md) | `architecture` | node 壁サーフェス出力の修正 (primal 面トポロジ退避 + h5 スキーマ + 壁関数 twall の ρu_τ² 化) |
| [thermophysics-multicomponent-tpgas.md](accepted/thermophysics-multicomponent-tpgas.md) | `thermophysics` | 多成分 thermally-perfect gas ソルバ機能 |
| [tooling-quasisteady-check.md](accepted/tooling-quasisteady-check.md) | `architecture` | 準定常確認ツール check_quasisteady.py と運用ルール (残差収束 ≠ 量の定常化) |
| [thermophysics-species-implicit-coupling.md](accepted/thermophysics-species-implicit-coupling.md) | `time_integration` | 化学種の陰解法緩和整合 (matched-relaxation scalar-DPLUR) |
| [time_integration-explicit-pointimplicit-sst.md](accepted/time_integration-explicit-pointimplicit-sst.md) | `time_integration` | 陽解法 RK のスカラー (k/ω) 源項 point-implicit 化 |
| [time_integration-general-eos-jacobian.md](accepted/time_integration-general-eos-jacobian.md) | `time_integration` | block-DPLUR の一般EOSヤコビアン化 (TP gas の陰解法律速の根治) |
| [time_integration-implicit-stable-cfl.md](accepted/time_integration-implicit-stable-cfl.md) | `time_integration` | 陰解法 (block DPLUR) の安定 cfl_pseudo 引き上げ |
| [time_integration-lowmach-preconditioning.md](accepted/time_integration-lowmach-preconditioning.md) | `time_integration` | 低マッハ前処理 (Weiss–Smith) — 密度ベース経路の散逸スケール是正と陰解法固有値前処理 |
| [time_integration-scalar-dplur-axisym-source.md](accepted/time_integration-scalar-dplur-axisym-source.md) | `time_integration` | scalar DPLUR の軸対称ソース Jacobian 整合 (TP 対応) |
| [turbulence-enhanced-wall-treatment.md](accepted/turbulence-enhanced-wall-treatment.md) | `turbulence` | SST Enhanced (Automatic / y⁺ 非依存) Wall Treatment |
| [turbulence-kato-launder.md](accepted/turbulence-kato-launder.md) | `その他 (turbulence)` | SST 生産項 Kato–Launder 補正 (`katoLaunder`) |
| [turbulence-node-sst-wallfunction.md](accepted/turbulence-node-sst-wallfunction.md) | `turbulence` | node-centered SST 壁関数の代表点修正 + τ_w 付与 (ic=壁ノード u=0 トラップ) |
| [turbulence-node-inlet-dirichlet-conserved.md](accepted/turbulence-node-inlet-dirichlet-conserved.md) | `turbulence` | node 入口 Dirichlet の保存量整合+残差除外 (rms_roOmega プラトー根治、cell 不変、x_R が実験に接近) |
| [turbulence-sst-freestream-init.md](accepted/turbulence-sst-freestream-init.md) | `turbulence` | SST 初期乱流の freestream 初期化 (kInit / omegaInit) |

## archived (superseded / 終了)

| Plan | area | 概要 |
| --- | --- | --- |
| [diffusion-node-scalar-nonortho-limit.md](archived/diffusion-node-scalar-nonortho-limit.md) | `diffusion` | node-centered flux の dcc に node 座標を使う (双対 CV 重心由来の見かけ非直交を排し SST omega 爆発を根治) |
| [precision-mixed-axisym.md](archived/precision-mixed-axisym.md) | `precision` | 混合精度 (iterative refinement) で軸対称 近軸の陰解法を root-fix する |

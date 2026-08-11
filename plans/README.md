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
| [convection-freestream-preserving-flux.md](active/convection-freestream-preserving-flux.md) | `convection` | 対流流束の free-stream 保存 (基準静圧差分) |
| [convection-keep-cb-pressure-correction.md](active/convection-keep-cb-pressure-correction.md) | `convection` | 高周波圧力欠陥 δp^HF 駆動の mass-flux 補正 (Rhie–Chow 型市松キラー, σ_min 独立で丘頂 spanwise 鋸歯を減衰) |
| [convection-keep-diss-lowmach-precond.md](active/convection-keep-diss-lowmach-precond.md) | `convection` | KEEP ES 散逸レイヤの低マッハ前処理化 (圧力ジャンプ増強で市松と解像物理を分離) |
| [convection-keep-revive-node.md](active/convection-keep-revive-node.md) | `convection / turbulence` | KEEP スキーム復活 (modern API・cell/node) + node WALE で LES/ILES |
| [discretization-lsq-gradient.md](active/discretization-lsq-gradient.md) | `gradient / architecture` | node-centered 最小二乗 (LSQ) 勾配 |
| [discretization-median-dual-3d.md](active/discretization-median-dual-3d.md) | `architecture` | 3D median-dual (M4): 3D 双対生成 + periodic 双対面 (3D node DDES/LES の起点) |
| [axisymmetric-su2-source-formulation.md](active/axisymmetric-su2-source-formulation.md) | `axisymmetric / discretization` | 軸対称を SU2 流 (planar 幾何+1/y ソース) に切替える `axisymMethod: 1`。node 軸行真空化の根治 (nodeAxisDirichlet 撤去) |
| [boundary-node-nozzle-wall-outlet-stability.md](active/boundary-node-nozzle-wall-outlet-stability.md) | `boundary / discretization` | node ノズル安定性。課題1/2/3 解決済 (§2.6-2.9)、**node 出口列不整合も根治 (§2.11, outlet Psb 動的化 — η 出口積分が正値に)**。残: 壁 T エントロピー市松 (§2.8) |
| [discretization-node-boundary-ghostless.md](active/discretization-node-boundary-ghostless.md) | `boundary` | node-centered 境界のゴースト撤廃 (段階導入: まず壁, 次に流入出/slip) |
| [tooling-cloud-gpu-env.md](active/tooling-cloud-gpu-env.md) | `tooling / infra` | クラウド GPU 環境 (AWS EC2 spot): Docker 構築・計算投入・速度計測基盤。MPI 化の実行基盤を先行整備 |
| [tooling-nozzle-design-tool.md](active/tooling-nozzle-design-tool.md) | `tooling / optimization` | 超音速ノズル設計ツール親計画 (5 機種・forge 評価器・サロゲート MOO・帰還エンジン・確認 CFD メニュー・AI 対話問題定義)。フェーズごとに子 plan を起票 |
| [tooling-nozzle-phase0-foundation.md](active/tooling-nozzle-phase0-foundation.md) | `tooling / optimization` | ↑の Phase 0 子 plan: 問題定義 YAML・区分構成ジオメトリ・TFI→forge h5 メッシュ・バッチ評価 CLI・目的関数ライブラリ |
| [turbulence-iddes-sst.md](active/turbulence-iddes-sst.md) | `turbulence` | SST-DDES / SST-IDDES 実装計画 |
| [turbulence-sst-thermal-flux-model.md](active/turbulence-sst-thermal-flux-model.md) | `turbulence / boundary` | SST 壁関数のエネルギー流束モデル置換 (Kader q_w)。等温壁×粗メッシュの熱負荷予測と T_aw 強閉包の前提 (in_progress: 平板合格 ±7%・Kader T⁺ 原式修正済。残 = T⁺ 圧縮性補正 [ベル +87% 実測]) |
| [turbulence-wmles-wall-stress.md](active/turbulence-wmles-wall-stress.md) | `turbulence / boundary` | WMLES 用代数壁応力モデル (Reichardt + Kader)。既存 SST 壁関数資産 (Normal_Neighbor / AddTauWall) を流用し τ_w/q_w で壁粘性流束を置換 |

## accepted (現役の設計判断)

| Plan | area | 概要 |
| --- | --- | --- |
| [turbulence-sst-adiabatic-taw-fluxmodel.md](accepted/turbulence-sst-adiabatic-taw-fluxmodel.md) | `turbulence / boundary` | SST 断熱壁 T_aw の W-I 流束注入 (全辺/代表辺とも) は**実測発散で棄却** (壁半CVの復元項喪失→EOS床まで異常冷却)。最終形 = Taw は境界出力 (Tsb/Taw_diag) 専用・W-I 拡散は常に DOF 状態・壁熱は境界半割面 q_w (断熱=厳密0)・res_roe[W] は生かす |
| [architecture-node-boundary-gradient-dof-only.md](accepted/architecture-node-boundary-gradient-dof-only.md) | `architecture / discretization` | node 境界勾配 (GG/LSQ) から bvar を排除し owner-state のみに統一 (node outlet P/T interior 化の一般化)。outlet 非退行・cell 不変を検証済。turbulence-sst-adiabatic-taw-fluxmodel (未完了) の前提 |
| [architecture-median-dual-3d-double-geometry.md](accepted/architecture-median-dual-3d-double-geometry.md) | `architecture / discretization` | 3D median-dual 幾何の Newell ローカル原点化+境界蓄積 double 化 (堅牢化)。**監査結論: 3D は元から double 演算で 2D のような実害なし** (露出見積もりを訂正)。wall_dist 定義は 2D と一貫 (双対重心間距離) |
| [architecture-limiter-negative-skip-fix.md](accepted/architecture-limiter-negative-skip-fix.md) | `architecture` | `limiter: -1` (off) が早期 return を素通りし Venkatakrishnan フル計算 (KEEP では未使用) に落ちるバグ修正。KEEP 系 run 全体で ~20% 高速化・挙動不変 |
| [turbulence-sst-omega-crossdiff-jacobian.md](accepted/turbulence-sst-omega-crossdiff-jacobian.md) | `turbulence / time_integration` | SST ω 交差拡散の point-implicit Jacobian (sstCrossDiffJac): dual-time サブ反復の ω 収束改善 (1.5x→2.7x)、収束解は不変 |
| [timeint-bodyforce-massflow-control.md](accepted/timeint-bodyforce-massflow-control.md) | `time_integration` | 体積力の質量流量一定制御 (Benocci & Pinelli 圧縮性版, γ=0.02): 周期丘 DDES 駆動。case/39 で検証済 |
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
| [architecture-perphase-profiling-hotspot.md](accepted/architecture-perphase-profiling-hotspot.md) | `architecture` | per-phase 計測駆動の host 律速最適化 (13.4→5.9s/2000step, −56%)。残レバー (CUDA Graph/renumbering) は不採用理由ごと記録 |
| [boundary-outlet-characteristic-backflow.md](accepted/boundary-outlet-characteristic-backflow.md) | `boundary` | 静圧固定流出 BC の特性ベース化と逆流統一 |
| [convection-keep-es-dissipation.md](accepted/convection-keep-es-dissipation.md) | `convection` | KEEP 用 entropy-stable 散逸レイヤ (scalar/matrix, CPG/TP/多成分, σ較正済) |
| [convection-keep-diss-recon-jump.md](accepted/convection-keep-diss-recon-jump.md) | `convection` | KEEP ES 散逸の再構成ジャンプ化 (keepDissJump: 市松無傷で解像スケール保護, node LES 推奨構成) |
| [turbulence-wale-fix.md](accepted/turbulence-wale-fix.md) | `turbulence` | WALE 不活性バグ (壁なしメッシュ wall_dist≡0) + Sd テンソル誤式の修正; 64³ TGV では WALE off の ILES+ES が最良と判明 |
| [turbulence-sigma-model.md](accepted/turbulence-sigma-model.md) | `turbulence` | σ-model SGS (Nicoud 2011, LESmodel:2) 追加; 64³ TGV では WALE 同様過散逸で ILES+ES が最良のまま |
| [convection-multispecies-contact-pressure.md](accepted/convection-multispecies-contact-pressure.md) | `convection` | 多成分 TP contact 圧力振動: S2 (`speciesFaceReconstruction:1`, ρ-limiter face 組成) を production 採用。S3/中心補間/PEP は棄却 |
| [convection-slau2-lowmach.md](accepted/convection-slau2-lowmach.md) | `convection` | SLAU2 圧力流束 (低マッハ圧力–速度カップリング是正) |
| [diffusion-node-boundary-real-distance.md](accepted/diffusion-node-boundary-real-distance.md) | `diffusion` | node 境界半割面 拡散の実距離 over-relax 化 (∇·S 弱形式の置換) |
| [diffusion-node-wall-viscous-distance.md](accepted/diffusion-node-wall-viscous-distance.md) | `diffusion` | node-centered 粘性壁 flux の法線距離修正 (mirror-ghost 退化対策) |
| [diffusion-turbulent-thermal-conductivity.md](accepted/diffusion-turbulent-thermal-conductivity.md) | `diffusion` | 乱流熱伝導 (turbulent thermal conductivity) の追加 |
| [diffusion-viscous-shear-flux.md](accepted/diffusion-viscous-shear-flux.md) | `diffusion` | 粘性せん断フラックスの修正 (内部面の法線項 + 壁面 no-slip 抗力) |
| [discretization-median-dual.md](accepted/discretization-median-dual.md) | `architecture` | cell/node (median-dual) 両対応化の親計画 (M1–M3 完了)。現役の入口は子 plan (median-dual-3d / ghostless / centroid-value-position / lsq-gradient) |
| [discretization-node-wall-implicit-dirichlet.md](accepted/discretization-node-wall-implicit-dirichlet.md) | `boundary / time_integration` | node-centered 壁 no-slip の implicit Jacobian Dirichlet 化 (SU2 DeleteValsRowi 相当) |
| [output-node-wall-surface-viz.md](accepted/output-node-wall-surface-viz.md) | `architecture` | node 壁サーフェス出力の修正 (primal 面トポロジ退避 + h5 スキーマ + 壁関数 twall の ρu_τ² 化) |
| [condensation-nonequilibrium.md](accepted/condensation-nonequilibrium.md) | `condensation` | 非平衡凝縮 (4 モーメント)。N2=Arthur Fig.2 1–2% / H2O=Wyslouzil Fig.3 ~5% 一致。TP<200K は棚上げ (CPG が正) |
| [thermophysics-multicomponent-tpgas.md](accepted/thermophysics-multicomponent-tpgas.md) | `thermophysics` | 多成分 thermally-perfect gas ソルバ機能 |
| [tooling-quasisteady-check.md](accepted/tooling-quasisteady-check.md) | `architecture` | 準定常確認ツール check_quasisteady.py と運用ルール (残差収束 ≠ 量の定常化) |
| [thermophysics-species-implicit-coupling.md](accepted/thermophysics-species-implicit-coupling.md) | `time_integration` | 化学種の陰解法緩和整合 (matched-relaxation scalar-DPLUR) |
| [time_integration-explicit-pointimplicit-sst.md](accepted/time_integration-explicit-pointimplicit-sst.md) | `time_integration` | 陽解法 RK のスカラー (k/ω) 源項 point-implicit 化 |
| [time_integration-general-eos-jacobian.md](accepted/time_integration-general-eos-jacobian.md) | `time_integration` | block-DPLUR の一般EOSヤコビアン化 (TP gas の陰解法律速の根治) |
| [gpu-implicit-plan.md](accepted/gpu-implicit-plan.md) | `time_integration` | GPU 陰解法基盤 (block-DPLUR / 軸対称 / SST point-implicit / scalar 版 / dual-time BDF2)。親計画 |
| [time_integration-implicit-stable-cfl.md](accepted/time_integration-implicit-stable-cfl.md) | `time_integration` | 陰解法 (block DPLUR) の安定 cfl_pseudo 引き上げ |
| [time_integration-lowmach-preconditioning.md](accepted/time_integration-lowmach-preconditioning.md) | `time_integration` | 低マッハ前処理 (Weiss–Smith) — 密度ベース経路の散逸スケール是正と陰解法固有値前処理 |
| [time_integration-scalar-dplur-axisym-source.md](accepted/time_integration-scalar-dplur-axisym-source.md) | `time_integration` | scalar DPLUR の軸対称ソース Jacobian 整合 (TP 対応) |
| [turbulence-enhanced-wall-treatment.md](accepted/turbulence-enhanced-wall-treatment.md) | `turbulence` | SST Enhanced (Automatic / y⁺ 非依存) Wall Treatment |
| [turbulence-kato-launder.md](accepted/turbulence-kato-launder.md) | `その他 (turbulence)` | SST 生産項 Kato–Launder 補正 (`katoLaunder`) |
| [turbulence-node-sst-wallfunction.md](accepted/turbulence-node-sst-wallfunction.md) | `turbulence` | node-centered SST 壁関数の代表点修正 + τ_w 付与 (ic=壁ノード u=0 トラップ) |
| [turbulence-node-inlet-dirichlet-conserved.md](accepted/turbulence-node-inlet-dirichlet-conserved.md) | `turbulence` | node 入口 Dirichlet の保存量整合+残差除外 (rms_roOmega プラトー根治、cell 不変、x_R が実験に接近) |
| [turbulence-sst-freestream-init.md](accepted/turbulence-sst-freestream-init.md) | `turbulence` | SST 初期乱流の freestream 初期化 (kInit / omegaInit) |
| [verification-passive-pseudoshock-control.md](accepted/verification-passive-pseudoshock-control.md) | `verification` | 多孔壁パッシブコントロールの逆解析 (case/36)。2D で中心軸波状 ~30% 低減=ショックレス化を定性再現 (Phase 1 完了) |

## archived (superseded / 終了)

| Plan | area | 概要 |
| --- | --- | --- |
| [diffusion-node-scalar-nonortho-limit.md](archived/diffusion-node-scalar-nonortho-limit.md) | `diffusion` | node-centered flux の dcc に node 座標を使う (双対 CV 重心由来の見かけ非直交を排し SST omega 爆発を根治) |
| [precision-mixed-axisym.md](archived/precision-mixed-axisym.md) | `precision` | 混合精度 (iterative refinement) で軸対称 近軸の陰解法を root-fix する |
| [turbulence-sst-thermal-wall-function.md](archived/turbulence-sst-thermal-wall-function.md) | `turbulence / boundary` | SST 壁関数の熱的閉包 (Crocco 型 T_aw「弱閉包」旧版)。**superseded**: 境界勾配が当時 bvar を読んでいたため場に触れる経路が実在し、かつ実データでは壁ノード実状態がほぼ動いていなかった (「4K 一致」は出力 Taw 診断値同士)。後継: `accepted/turbulence-sst-adiabatic-taw-fluxmodel.md` |

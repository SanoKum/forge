# plans/active 監査: 「弱いモデルでも実行できるか」の観点 (2026-07-19)

## 目的と経緯

Fable (最上位モデル) が使えなくなる前提で、`plans/active/` 全 14 件を「Opus/Sonnet 級のモデルが
plan だけを頼りに安全に実行できる状態か」の観点で監査した。判定軸は次の 4 つ:

- **(a) 次の一手**が具体的か (触るファイル・関数・手順の順序)
- **(b) 検証コマンドと数値合格基準**が書かれているか (モデルの主観判断に頼らない)
- **(c) 却下済み代替案と理由**が書かれているか (弱いモデルの再提案・再発明を防ぐ)
- **(d) ユーザー判断が必要な分岐**が明示されているか (勝手に進んではいけない箇所)

**注意**: 監査時点で `plans/active/*.md` に別セッションの未コミット編集があるため、本監査は
plan 本体を編集せず所見のみ記す。反映は当該セッション終了後に行うこと。

## 総括

リポジトリ全体としては AGENTS.md の VERDICT 必須ルール + check_* ツール群が強力な防波堤であり、
plan の多くは (b) を満たす。最大のリスクは次の 3 点:

1. **実質完了 plan が active に滞留**している (7 件)。弱いモデルは「active = やるべき仕事」と
   解釈して蒸し返す危険がある。→ **クローズ作業 (accepted へ移動 + README 同期) を最優先で実施**。
   これは機械的作業で、弱いモデルでも安全にできる。
2. **ユーザー設計判断待ちの分岐が plan 内に埋もれている** (centroid-value-position の
   omega BC 再定式化、lsq-gradient の動機再検証)。Fable がいるうちに判断を確定して plan に
   書き込むのが最も利回りが高い。
3. **plan 間の陳腐化した相互参照** (iddes-sst §4.8 の KEEP_SLAU 設計 vs その後確定した
   KEEP+ES 散逸系、keep-revive-node の WALE 前提 vs wale-fix の「WALE off 推奨」)。
   弱いモデルは古い節をそのまま実装しかねない。→ 整合メモの追記が必要。

## 分類サマリ

| Plan | 実質状態 | 弱モデル実行可否 | 主な欠落 |
| --- | --- | --- | --- |
| architecture-perphase-profiling-hotspot | **完了** (「一区切り推奨」明記) | クローズのみ | status/移動/README 同期 |
| gpu-implicit-plan | **完了** (dual-time まで検証済) | クローズのみ | 残 2 件を別 plan 化 or 明示棚上げ |
| convection-multispecies-contact-pressure | **結論確定** (S2 採用/S3 棄却) | クローズ作業可 | docs 反映・回帰登録・チェックボックス消化 |
| convection-freestream-preserving-flux | **主要部完了** (pRef+advGauge) | クローズ作業可 | 残課題 (TP/node/SLAU/種) の別 plan 化、変更ログ重複整理 |
| condensation-nonequilibrium | **Phase 1–3 完了** | クローズ作業可 | TP<200K 不安定の扱い決定 (別 plan 化 or 棚上げ明記) |
| verification-passive-pseudoshock-control | **主目的達成** (~30% 低減再現) | 要ユーザー判断 | 続行 (3D 丸穴等) か Phase 1 クローズかの決定 |
| discretization-median-dual (親) | **M1–M3 完了・M4 は子 plan** | クローズ作業可 | 未確認 1 点: case/29 node viscous が ghostless 修正後に回るか |
| discretization-median-dual-3d | 進行中 (現ブランチ) | **可** (良い見本) | 残タスク 3 件の優先順のみ |
| discretization-node-boundary-ghostless | 進行中 | 条件付き可 | 5c の要否判定が宙吊り、species+node 未検証の検証ケース未指定 |
| architecture-node-centroid-value-position | 進行中 (Stage1+2 済) | **不可 (判断待ち)** | node 壁 omega BC の再定式化 = ユーザー判断。残ステップの分割も粗い |
| discretization-lsq-gradient | 負結果で停止 | **不可 (判断待ち)** | 動機 (checkerboard 実害) の再検証が先、というゲートを plan 冒頭に明示すべき |
| convection-keep-revive-node | ほぼ完了・一部陳腐化 | 条件付き可 | WALE 前提が wale-fix/ES 散逸の結論と不整合、backstep LES 検証の残否を明確化 |
| turbulence-iddes-sst | Phase 1 完了・1.5 待ち | 条件付き可 | §4.8 と KEEP ES 散逸系 (keepDissType/Jump) の整合メモが必要 |
| turbulence-sigma-model | 起票直後 | ほぼ可 | 実装ステップのファイル粒度と合格基準数値の明記 |

## 個別所見

### 1. クローズすべき 7 件 (機械的作業・弱モデルで安全)

**共通手順** (AGENTS.md の完了時フローそのまま): ① `status: done` 化、② 変更ログに最終状態を
1 段落追記、③ `plans/active/` → `plans/accepted/` へ `git mv`、④ `plans/README.md` の両セクション
同期、⑤ related_docs の methods/ 整合確認。

- **architecture-perphase-profiling-hotspot**: 変更ログ末尾に「性能フェーズはここで一区切り推奨」
  「残レバー (DPLUR 帯域/SLAU/CUDA Graph) はコスパ低」と結論済み。累計 13.38→5.9s (−56%)。
  残レバーは「やらない理由」ごと accepted に残せば再提案を防げる。
- **gpu-implicit-plan**: block-DPLUR・軸対称・SST point-implicit・scalar 版・dual-time まで全て
  実装・検証済み。残 (dual-time の scalar 対角版、渦放出系の非定常検証) は必要になったときの
  別 plan とし、その旨を明記してクローズ。
- **convection-multispecies-contact-pressure**: §13–15 で S2 production / S3 棄却 / S2c 劣位まで
  決着し、原理 (「ρ と Y は同一再構成で揃える」) も 3 実験で確立。未消化は §9 の docs 反映
  (`methods/convection/{theory,implementation}.md`)、回帰登録、commit 分離のみ。
  **production 設定 (`speciesFaceReconstruction: 1`) を `procedures/solver-settings.md` にも載せる**こと
  (弱いモデルが設定値を推測しないための一次情報)。
- **convection-freestream-preserving-flux**: pRef (SLAU/KEEP) + advGauge (KEEP CPG×cell) が実装・
  機械精度検証済み。残課題 (ROE/AUSM への pRef、SLAU の移流ゲージ、TP 枝、node 境界、種輸送) は
  優先度が低い旨を明記して followup 別 plan に切り出し。**ファイル末尾に「## 変更ログ」が 2 つ
  重複**しているので統合する。`unsteady:0` の音響過渡不安定は memory 化済み
  ([[steady-localdt-acoustic-transient-instability]]) だが `procedures/divergence-and-startup.md`
  への転記を確認。
- **condensation-nonequilibrium**: N2 (Fig.2 1–2% 一致)・H2O (Fig.3 ~5% 一致・分圧汎化) とも検証済。
  未決は「TP 気相の冷・高マッハ域数値不安定」(真因未特定、フラックス計装が要る) — これは
  凝縮ではなく convection/thermophysics の問題なので、**別 plan として切り出すか「>200K で
  実用可のため棚上げ」と明記**してクローズ。
- **verification-passive-pseudoshock-control**: 完了条件の主項目「ショックレス再現可否」は
  **達成済** (Ps=1.82–1.84 窓で ~30% 低減)。ただし変更ログに新旧 2 セッションの結論が併存し
  「未達」(2026-06-21 自走) と「再現」(2026-06-21 続き) が読み手を混乱させる。**冒頭に最終結論の
  サマリ 3 行を追記**した上で、続行項目 (3D 丸穴・10°ディフューザ・matched-x_f) を続けるかは
  ユーザー判断。判断が出るまでは触らないこと。
- **discretization-median-dual (親)**: M1–M3 の本文は歴史的価値が高いが、M4 は
  `discretization-median-dual-3d` に、境界は `node-boundary-ghostless` に、centCoords は
  `node-centroid-value-position` に分離済み。親としてクローズし「現役の入口は子 plan」と明記。
  **クローズ前の唯一の確認事項**: case/29 node viscous (軸対称) の exit-lip 発散が、その後の
  ghostless 5a/5b/5e + 弱形式化で解消したかの再検証 run が見当たらない。1 run 流して VERDICT を
  取ってから閉じるのが望ましい (未解決なら ghostless plan 側に残課題として移す)。

### 2. ユーザー判断を先に確定すべき 2 件 (Fable/ユーザーの残り時間をここに使う)

- **architecture-node-centroid-value-position**: 変更ログ 2026-06-22 の未解決 3 項目のうち
  「**SST 壁 omega を壁ノード y=0 でなく第一内部ノード距離で定式化する**」は明示的に
  「ユーザーと方針決定」待ち。ここが決まらないと Step 1 (ghostless 完遂) 以降が着手できない。
  **推奨**: 方針決定を先に行い、決定内容 (式・参照実装・却下案) を plan §4 に書き足してから
  実装セッションに渡す。なお ghostless plan 側 (2026-06-22) で「node は omega[ic] を直接ピン +
  wall_y_eff MIN + point-implicit decouple」が既に実証済みなので、centCoords=node 化後も
  同方式 (壁 omega は第一オフ壁距離基準) を踏襲するのが自然、という整合メモも足すとよい。
  残ステップ (ghost 生成停止 `nCells_all=nCells` / convert centCoords=node / r_eff 分離 /
  axisCentroidShift 撤去) は §5 にあるが、**各ステップの検証 run とビット不変確認の対応表**
  (どの run を流し何を比較するか) を §6 の粒度まで割り付けてから着手させること。
- **discretization-lsq-gradient**: 負結果 (近壁 M 退化で発散) で停止中。修正候補 ①double 蓄積
  ②重み/QR ③GG フォールバック、が並ぶが、**その前段の「動機の再検証」(GG が本ケースで収束して
  おり、そもそも近壁 checkerboard の実害がどのケースで出ているのか) が本文中で「再確認すべき」の
  ままゲート化されていない**。弱いモデルは修正候補①〜③に直行しがち。**推奨**: plan 冒頭に
  「次の一手 = 固定 GG 場で checkerboard 指標を測る診断 (安価)。実害が確認できなければ
  archived (動機消滅)」というゲートを明記する。なお `nodeMidpointFx` (fx=0.5) が checkerboard を
  −36% 低減した結果 (median-dual plan 2026-06-14) との関係も整理しておくと二重投資を防げる。

### 3. 進行中 3 件の整合性メモ (陳腐化した節の明示)

- **convection-keep-revive-node**: §1 の目的「KEEP + WALE で LES」は、その後の
  `turbulence-wale-fix` の結論 (**64³ TGV では WALE off の ILES + ES 散逸が最良**) および
  `convection-keep-es-dissipation`/`keep-diss-recon-jump` (σ=0.02 + keepDissJump:2 推奨) と
  部分的に矛盾する。§7 (陰解法化) は §7.6 で決着済み (新規実装不要)。**推奨**: 「現在の推奨
  LES 構成は KEEP + ES 散逸 (jump2, σ0.02) + WALE off。WALE 併用は σ-model (turbulence-sigma-model)
  の結果待ち」という現在地サマリを冒頭に追記し、残作業 (backstep 実 LES の定量検証があるなら
  それ) を列挙してから渡す。実質クローズ可能な可能性が高い。
- **turbulence-iddes-sst**: Phase 1 完了・検証ログ充実で品質は高い。ただし §4.8 (Phase 1.5
  低散逸 flux) の本命設計「`KEEP_SLAU` の `duc` に f_d 注入」は 2026-06-21 時点の設計で、
  その後 KEEP 系は純化 (`keepDissipation` 廃止) + ES 散逸レイヤ (keepDissType/keepDissJump) に
  再編された。`KEEP_SLAU_d` 経路が現存するか・ES 散逸との住み分け (σ を f_d で駆動する案も
  あり得る) を**実装前に再調査して §4.8 に追記**しないと、弱いモデルが存在しないコードパスを
  前提に実装する危険がある。§5.7 の失敗予測と T1-B 判定基準 (RMS・−5/3・スパン相関) は
  そのまま強力なガードなので維持。
- **discretization-node-boundary-ghostless**: mean-flow (5a/5b/5e)・SST (omega pin + wall_dist +
  massflux) まで決着し node SST が cell 基準一致。宙吊りは ① **5c (壁粘性 CSR カーネル) の要否**
  — u_τ 問題は 5c なしで解決したため、5c は「精度改善オプション」に格下げか撤回かを明記する。
  ② **species+node の未検証** (2026-06-26 追補で「構造は k/ω パターン準拠だが未検証」) —
  検証に使うケース (多成分 TP × node の run が存在しない) の指定がない。case/28 の node 版など
  検証ケースを指定するか「node×species は使用前に要検証」と solver-settings 側に警告を書く。

### 4. 良い見本と新規 plan

- **discretization-median-dual-3d**: 本監査で最も「弱いモデルが実行できる」形式に近い。
  設計判断の理由 (§4.5.1 なぜ素朴案が駄目か)・却下案 (§4.5.5)・実装状況の逐次追記・検証の
  数値基準が揃う。残: §4.5.8 回転周期 / x_R settle / wall_y_eff 修正の SST 較正影響。
  この 3 件に優先順位と「どこまでやったら done か」を 1 行ずつ足せば完全。
- **turbulence-sigma-model**: スコープ・事前検証 (2 成分流で ν_t≡0) は明確。弱いモデル向けに
  足すべきは ① 実装ステップのファイル粒度 (`turbulent_viscosity_d.cu` に `SIGMA_d` 追加、
  `solverConfig` の `LESmodel` 範囲拡張、WALE_d の起動グリッド踏襲 = normalcell)、
  ② 合格基準の数値 (TGV 64³ node で KE 散逸ピーク t* が WALE の 7.84 より DNS 8.98 側へ、
  層流期 ν_t/μ ≪ WALE の 0.26)、③ 比較対象 run の指定 (case/09 の既存 ILES+ES run を対照に)。

## 横断アクション (優先順)

1. **クローズ 7 件の実施** (§1)。弱いモデルで可。1 件ずつ commit を分ける。
2. **ユーザー判断 2 件の確定** (§2: omega BC 定式化 / LSQ 動機ゲート) → plan へ書き込み。
3. **整合メモ 3 件の追記** (§3)。KEEP_SLAU 現存確認は grep 1 回で済む。
4. `plans/README.md` に `turbulence-sigma-model` を追記 (active 表に未掲載。別セッションが
   対応中の可能性があるため要確認)。
5. 検証系メモ: 監査中に見つけた**未検証のまま残っている主張**は
   (a) case/29 node viscous の ghostless 後再検証、(b) species×node 経路、
   (c) `run_0021` 系 prism/pyramid の 3D 双対 (同一コードパスだが未実行)。
   回帰スイート拡充の候補にすること。

## 追記 (2026-07-19 反映後の訂正)

- **turbulence-sigma-model の所見は撤回**: 監査時に読んだ `plans/active/turbulence-sigma-model.md`
  (41 行, in_progress) は stale な下書きの重複だった。実体は commit 92fd11f で**最初から
  `plans/accepted/` に完了済みとして追加**されており (結論: σ-model も 64³ TGV では WALE 同様に
  過散逸で、ILES+ES 散逸が最良のまま)、README の accepted 表にも掲載済み。active 表への追記は不要
  だった (反映作業中に是正済み)。
- 反映済みアクション: クローズ 6 件 (perphase / gpu-implicit / condensation / multispecies+docs /
  pseudoshock Phase1 / median-dual 親)、整合メモ 3 件、ユーザー判断 2 件の plan 書き込み
  (omega BC 踏襲・LSQ ゲート)、case/29 再検証の ghostless への移管。**未反映**:
  convection-freestream-preserving-flux のクローズ (残課題の別 plan 化含む) は当該セッションの
  依頼により保留 (a648e50 後の内容を前提に実施すること)。

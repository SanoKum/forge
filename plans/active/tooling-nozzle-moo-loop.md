# ノズル設計ツール Phase 2: ③ベル MOO ループ (TOP 幾何 dv + pymoo/SMT + Rao 照合)

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (現在仕様)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) — §4.2 (pymoo+SMT)・§4.6③ (TOP 幾何 dv, 2026-08-13 改訂)・§6 Phase 2 を実装する
  - 前提: [`tooling-nozzle-phase0-foundation.md`](tooling-nozzle-phase0-foundation.md) (評価パイプライン。帰還エンジンには依存しない)
- **created**: `2026-08-13`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

親計画 Phase 2 の実体化。③ベルを **TOP 直接幾何 dv** ($\theta_n, \theta_e, L/r_t$ の 3 変数) で
パラメータ化し、DOE → Kriging サロゲート → NSGA-II → EHVI infill のループを回して
η ($C_F/C_{F,ideal}$) vs $L/r_t$ の 2 目的パレート前線を取得する。
**最重要マイルストーン = Rao TOC/TOP 照合** (前線が Rao 点と整合し、不当に支配しないこと)。
検証対象は評価器 + MOO ループであり、帰還エンジンは関与しない (①の CONTUR 照合が担当 — 親計画 §4.7)。

## 2. スコープ

- **やる**:
  - **TOP 放物線壁**: `geometry/wall.py` の下流壁に `bell_type: top` を追加 ($P_a$ で接線
    $\theta_n$・出口 $(L,\sqrt\varepsilon)$ で接線 $\theta_e$ の 2 次 Bézier。C1 接続、$P_a$ 曲率不連続は
    Rao 慣行通り許容)。既存 3 次 Hermite は `bell_type: hermite` として残す (回帰対照)
  - `probdef`: dv に `theta_e_deg` 追加、TOP の幾何整合チェック (接線交点が $P_a$–出口間に存在)
  - **opt/ モジュール新設**: LHS DOE (`doe.py`) / SMT `KRG` サロゲート (`surrogate.py`) /
    pymoo NSGA-II サロゲート上探索 (`moo.py`) / **EHVI infill 自作 (2 目的解析式)** (`infill.py`)
  - **バッチドライバ**: 評価キュー (逐次, 1 GPU)、warm seed 管理 (幾何最近傍の収束済み run から
    `interp_field`、発散時は段階起動 fallback)、`check_mesh_quality` / `check_convergence` /
    `check_quasisteady` の VERDICT ゲート (PASS 以外はサロゲートに入れず評価失敗扱い)
  - **Rao 照合**: case/29 の Rao/円錐資産 + TOP チャート点 (θn, θe の慣行値) を DOE に混ぜ、
    パレート前線 vs Rao 点を比較・図化
- **やらない** (後続):
  - 帰還エンジン v1–v3 (Phase 3 = ①風洞と共起票)、モード K・Kliegel–Levine (③2 巡目の任意課題)
  - co-kriging 多忠実度 (Euler 低忠実度層は①で導入)、3 目的以上の MC-HVI、並列バッチ評価
  - 剥離制約の自動判定 (Summerfield) は骨格のみ (今回の動作点は設計背圧で剥離しない)

## 3. 関連 docs と前提

- 評価パイプライン: Phase 0 の runner (node 生産レシピ = 全域 2 次 + 陰解法 cfl_pseudo 4 +
  nodeAxisDirichlet + warm start、12000 step ≈ 20 秒)。E2E 回帰 = run_0071 (2026-08-13,
  現バイナリ, η=0.9734, ALL STEADY)
- η の絶対値はメッシュ/壁処理依存が残る (case/40 README: y+1 で 0.978±0.003) が、本 Phase の
  目的は**同一メッシュ規約での形状間ランキング** — DOE 全点で同一メッシュパラメータを使い、
  系統誤差を共通化する
- Python 依存 (pymoo / SMT) は `.venv-opt` を新設して隔離 (システム python3 は Phase 0 通り
  numpy/scipy/h5py/yaml のみ)

## 4. 設計方針 (差分のみ)

- **EHVI は 2 目的解析式で自作** (親計画 §4.2 案 A)。Emmerich 系の区分積分で実装し、
  **ZDT1/ZDT2 ベンチで BoTorch qNEHVI (または文献値) と照合**してから実戦投入。
  照合で挙動が疑わしければ A' (BoTorch 借用) へ昇格 — 親計画の条項通り
- **決定性**: DOE seed・infill 候補生成・サロゲート学習の乱数を YAML 固定。同一 YAML → 同一
  提案点列を単体テストで保証
- **評価失敗の扱い**: 幾何フィルタ NG / メッシュ FAIL / 発散 / NOT STEADY は「実行不能点」として
  optimizer に返す (ペナルティでなく制約違反マーク — pymoo の feasibility 扱い)
- **run 台帳**: 評価 run は `case/40.nozzle_design_tool/` に `run_NNNN_moo_<tag>` で置き、
  README 一覧表を維持。ただし DOE 数十点は 1 行にグルーピングして記載 (点ごとの行は作らない)

## 5. 実装ステップ

1. `wall.py` TOP 放物線 + 単体テスト (接線一致・C1・決定性・整合チェック)
2. `probdef` dv 拡張 + YAML 雛型 (`problem_bell_top.yaml`) + E2E 1 点スモーク (TOP 壁で PASS)
3. `.venv-opt` + `opt/` 骨格 (doe / surrogate / moo / infill) + EHVI ZDT ベンチ照合
4. バッチドライバ (評価キュー + VERDICT ゲート + warm seed 管理)
5. DOE (LHS ~30 点) 実行 → サロゲート → infill (~20 点) → パレート
6. Rao 照合 (case/29 資産) + 図化 + 親計画 §6 Phase 2 判定

## 6. 検証

- **単体**: TOP 壁の端点接線・決定性、EHVI の ZDT1/ZDT2 照合 (hypervolume 進展が文献レンジ)、
  提案点列の再現性
- **E2E**: TOP 壁 1 点で VERDICT PASS (収束 + 準定常)
- **Rao 照合 (合格基準)**: (i) パレート前線が Rao/TOP チャート点を ±(評価器ノイズ ~0.3%) で
  通る、(ii) Rao 点を有意に支配する点が出た場合は評価系のバグとしてまず疑い、原因を特定する
- 全評価 run で 3 ツール VERDICT を自動記録 (§4 ゲート)

## 7. 影響範囲

- 新規: `design/forge_design/opt/`、`design/.venv-opt` (git 管理外)、`problem_bell_top.yaml`
- 変更: `geometry/wall.py` (bell_type 追加)、`probdef.py` (dv 追加)、`evaluate/runner.py` (seed 管理)
- 既存ソルバ無改造。`methods/design/overview.md` を実装と同期

## 8. 完了条件

- [ ] §5 の 1–6 実装・実行
- [ ] Rao 照合の合格判定 (§6) を親計画 §6 Phase 2 に記録
- [ ] `methods/design/overview.md` 同期
- [ ] status を done にして `accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-08-13` — §5 ステップ 4 完了: `opt/driver.py` (バッチドライバ)。スモーク
  (`run_0073_moo_smoke`, DOE4+infill1) で 2 つの実装欠陥を発見し修正した:
  (1) **cross-geometry warm start の陰解法 cfl4 直投入は発散せず「不始動流」という偽アトラクタ
  に収束し得る** (machmax~1.4・CF~7・ṁ −13%) — 残差・quasisteady だけでは弾けないため、
  **物理ゲート (ṁ/ṁ_1D ∈ [0.94,1.005], η ∈ (0.5,1.02), res_step=規定 step) を追加**し、
  **全評価を soft 段 (convMethod 0+cfl0.5, 3000 step) → 本段の 2 段起動に統一** (run_0019/0034
  レシピの常用化)。(2) soft 段の `outStepInterval` が段長を超え最終場が未ダンプ → 本段が seed
  直投入相当になる罠 (+ quasisteady の「NOT ALL STEADY」部分文字列誤パース) を修正。
  修正後スモーク = **5/5 PASS** (~33s/評価, conv PASS・ALL STEADY・ṁ 比 0.9858±0.0003,
  η は L に単調応答 0.952@L5.1→0.978@L10.9)。本番キャンペーン `run_0074_moo_top_r1`
  (DOE24+infill2×13=50 評価, seed=run_0072, ref=(−0.90,13)) を投入。
- `2026-08-13` — §5 ステップ 3 完了: `opt/` モジュール新設 (`ehvi.py` = 2 目的 EHVI 閉形式
  [Hupkens–Yang–Emmerich 区分和] + hypervolume + 非劣解 / `doe.py` = LHS / `surrogate.py` =
  SMT KRG 束 [既定で決定的・補間的を確認] / `moo.py` = NSGA-II サロゲート探索 + EHVI infill
  提案 [候補 = 最終個体群 ∪ LHS 探査点、正規化距離で重複排除])。**検証 ALL PASS**
  (`design/tests/run_opt_tests.py`, .venv-opt): EHVI は (i) 空前線の解析解一致、(ii) σ→0 極限 =
  決定的 HVI 厳密一致、(iii) **Monte-Carlo 照合 4 配置で ≤0.1% 一致**、(iv) **ZDT1 ミニ BO
  (d=4, DOE24+infill2×14=52 評価) で HV 96.6→120.62 = 真値 120.667 の 99.96%**。LHS/infill の
  同一 seed 決定性・KRG 補間性も PASS。§4 の「A' 昇格」条項は不発 (自作 EHVI で足りる)。
- `2026-08-13` — §5 ステップ 1–2 完了: `wall.py` に `bell_type: top` (TOP 放物線 = 2 次 Bézier。
  $t(x)$ 逆変換は Citardauq 型で $A\to0$ 安定、整合条件 $\tan\theta_n >$ 弦勾配 $> \tan\theta_e$ を
  構築時検査)。単体テスト 10 本追加 ALL PASS (端点/接線/C1/凸性/数値微分整合/NG 検出)。runner の
  `theta_e_deg` dv 化 + `problem_bell_top.yaml`。E2E スモーク
  `case/40.nozzle_design_tool/run_0072_bell_top_smoke` = 品質 PASS / 収束 **ALL PASS**
  (3.4–6.0 桁) / quasisteady **ALL STEADY**、η=0.9731・ṁ=1.2928 (同条件 Hermite run_0071 と
  −0.03%/一致)。
- `2026-08-13` — 起票 (親計画 §4.6③ の 2026-08-13 改訂 = TOP 直接幾何 dv を受けて、旧予定の
  「モード K + 帰還前提の MOO」から「TOP 幾何 dv の MOO」にスコープ確定)。

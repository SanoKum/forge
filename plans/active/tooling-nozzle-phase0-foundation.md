# ノズル設計ツール Phase 0: 基盤 (問題定義 YAML・ジオメトリ・TFI メッシュ・評価 CLI・目的関数)

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (現在仕様)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) (§4.6/§4.7 の決定を実装する)
- **created**: `2026-08-03`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

親計画 Phase 0 の実体化。設計変数ベクトルから「壁形状 → forge メッシュ → 計算 → 目的関数値」
まで無人で辿れる評価パイプラインの土台を作る。帰還エンジン (Phase 1)・最適化ループ (Phase 2)
はこの上に載る。

## 2. スコープ

- **やる**:
  - `design/` パッケージ骨格 (`forge_design`)
  - 問題定義 YAML スキーマ + 検証 (spec/derived/dv、過拘束チェック)
  - ジオメトリ生成器: 区分構成壁 (Bell–Mehta 収縮 + 円弧 $R_u/R_d/\theta_a$ + 下流壁)、
    中心線マッハ Bézier (端点拘束・自由度勘定)、Kliegel–Levine 遷音速級数、kernel MOC マーチ
    (アンカー $x_k, M, M', M''$)、事前フィルタ
  - 下流壁の**初期形状**: 擬 1D 面積則 (目標 $M_c$ → 面積比 → 半径) の C1 接続。
    完全な MOC 逆設計は実装しない (帰還 Phase 1 が収束させるため初期値精度は収束速度にのみ影響
    — 親計画 B3 の理屈。将来必要なら別 plan)
  - TFI 構造化メッシュ生成 → forge HDF5 直書き (トポロジ固定・壁クラスタリング・品質ゲート)
  - バッチ評価 CLI (run dir 準備 → forge 起動 → 収束/NaN 判定 → メトリクス出力)
  - 目的関数ライブラリ最小セット: ③ベル用 ($\eta = C_F/C_{F,ideal}$, $L/r_t$) + 軸上/出口
    分布抽出 ($\varepsilon_M$, $\varepsilon_\theta$ は骨格まで)
- **やらない** (後続フェーズ):
  - 帰還エンジン (Phase 1)、pymoo/SMT 最適化ループ (Phase 2)
  - ②側壁 δ\*・④2 作動点・⑤壁圧/3D (Phase 4)
  - $q_{peak}$ の壁関数/壁解像ランキング検証 (Phase 3)

## 3. 関連 docs と前提

- 仕様: [`methods/design/overview.md`](../../methods/design/overview.md)
- 判断の出典: 親計画 §4.1 (YAML)、§4.3 (R1 = TFI トポロジ固定)、§4.6 (dv・目的の初期セット、
  Bézier 自由度勘定、閉ループ派生)、§4.7 (壁所有権・アンカー)
- forge 側は無改造 (評価器としてブラックボックス利用)

## 4. 設計方針 (差分のみ)

- Python は WSL native の `python3` (numpy/scipy/h5py/yaml 確認済み)。追加依存なしで動く
  範囲を Phase 0 とする (pymoo/SMT は Phase 2 で venv 導入)。
- forge h5 スキーマは既存コンバータの出力仕様に厳密一致させる (調査結果を
  `design/meshing/` 実装コメントと本 plan 変更ログに記録)。
- メッシュは当面 cell 中心を既定 (③ベルの Rao 照合資産 [case/29] が cell 前提のため)。
  node (median-dual) は変換オプションとして温存 ([user-prefers-node-base] は認識しつつ、
  Phase 0 の照合は既存資産との一致を優先)。
- 決定性: 同一 YAML → bit 同一メッシュを単体テストで保証する。

## 5. 実装ステップ

1. `design/forge_design/` 骨格 + `probdef.py` (YAML スキーマ・検証)
2. `geometry/`: `bezier.py` → `transonic.py` (Kliegel–Levine) → `moc_kernel.py` →
   `wall.py` (区分構成 + 擬 1D 初期壁)
3. `meshing/`: `tfi.py` → `forge_writer.py` (h5 直書き) → 品質ゲート呼び出し
4. `metrics/`: `extract.py` (res h5 読取・軸/出口分布・$C_F$/$\eta$/$L$)
5. `evaluate/`: `runner.py` (run dir 準備・forge 起動・判定) + CLI エントリ
6. 検証 (§6) → 親計画 §5 表の Phase 0 行を本 plan にリンク

## 6. 検証

- **単体**: Bézier 端点条件 (指定次数で $M,M',M''$ 一致・単調性 repair)、Kliegel–Levine の
  既知値照合 (壁マッハ・流量係数のオーダー)、TFI の決定性 (同一入力 → 同一 sha256)
- **メッシュ**: 生成した③ベルメッシュが `check_mesh_quality.py` で VERDICT PASS
  (AR≤1000 / skew≤0.9)
- **E2E スモーク**: 生成メッシュ + 既存 config 雛型で forge を実行し完走・収束
  (`check_convergence.py` VERDICT を記録)。run は新規 case ディレクトリ
  (`case/40.nozzle_design_tool/`) の `run_*` に置き README run 一覧表を作成
- **判定基準**: 上記 3 点全て PASS。$C_F$ 抽出値が 1D 理論値と整合するオーダー
  (詳細照合は Phase 2 の Rao 照合で行う)

## 7. 影響範囲

- 新規: `design/` 一式、`case/40.nozzle_design_tool/`、`methods/design/overview.md`
- 既存コードへの変更なし。`methods/index.md` に design 行を追加済み

## 8. 完了条件

- [ ] §5 の 1–5 実装
- [ ] §6 の検証 3 点 PASS (VERDICT 記録)
- [ ] `methods/design/overview.md` を実装と同期
- [ ] status を done にして `accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-08-03` — 起票 (親計画 Phase 0 の子 plan)。実装開始。
- `2026-08-03` — Phase 0 実装・検証 (初日, Claude 自走):
  - **メッシュ経路の実装判断**: forge h5 直書きでなく **msh4.1 テキストを Python から直書き →
    既存 `convertGmshToForge`** に変更 (`meshing/mesh2d.py`)。wall_dist・幾何量・境界構造を
    検証済み変換器に任せられ、スキーマ追従リスクを回避。R1 の本質 (トポロジ固定・決定的
    再生成) は不変。前例 = case/23 `generate_3d_ogrid_msh.py`。
  - **初期下流壁**: 擬 1D 面積則でなく **3 次 Hermite ベル** (P_a で C1、出口半径 √ε・出口角
    θ_e) に変更 — 実装が単純で、初期形状の役割 (帰還の出発点) には十分。
  - **単体テスト ALL PASS** (`design/tests/run_tests.py`): Bézier 端点/自由度勘定、Sauer 壁
    接線 (幾何スロートと軸ソニック点のオフセット込み)、壁 C1/C2 接続、メッシュ/msh の
    bit 決定性、PM 往復、**kernel MOC 3 系統検証** (放射源流厳密解 0.1%・2D 平面単純波則
    1%・軸対称格子収束 <0.5%)。
  - **kernel MOC の実装知見**: (i) 軸対称 1/r 源項は「軸から遠い方の点」で評価しないと
    δθ/r 増幅で発散、(ii) 壁ステップは規定ステーション + C+ の足を前フロント補間で求める
    方式が必須 (近壁点間隔任せは破綻)、(iii) 初期値線の楔領域は「初期線各点からの C- を
    順にフロント化」する充填が必要 (鉛直線直結は O(Δ²) バイアス、Sauer 場での C- トレースは
    有効域外挿で破綻)。**既知の制約: Rd ≲ 1 (ロケット慣行 0.4) は Sauer パッチが円弧終端を
    超え明示エラー — Kliegel–Levine 高次級数の実装が Phase 1 以降の宿題** (モード K を
    ③に使う前に必要)。
  - **メトリクスの定義修正**: 出口面運動量法 $F=\int(\rho u^2+(P-P_a))\,dA$ は壁摩擦・壁圧の
    影響を出口状態経由で既に含むため、サーベイ B4.5 の「+壁摩擦抗力」加算は**二重計上として
    不採用** (入力 h5 の twall はソルバ非更新の全ゼロで、実害が出る前に発見)。
  - **dv 感度確認** (run_0003/0004): L/rt = 6/7/9 → η_CF = 0.9708/0.9790/0.9837 の単調応答、
    ṁ は 0.03% で一定 (スロート同一)。評価器としての応答性を確認 — Phase 2 パレート軸の実証。
  - **E2E スモーク** (`case/40.nozzle_design_tool/`, README に run 一覧):
    run_0001 = SST roK/roOmega 未初期化で step 2–6 NaN (既知指紋を再現) → IC パッチに
    roK/roOmega 追加 (`evaluate/ic.py`) で run_0002 = **NaN なし完走**。メッシュ品質
    VERDICT PASS (AR 19/skew 0.42)。残差は 1.1–1.7 桁で プラトー (VERDICT: NOT CONVERGED/
    stalled — cell モードノイズ床水準) だが**推力は 4000→12000 step で 0.002% ドリフト**
    (F=1961 N, η_CF=0.979, ṁ=1.299 kg/s)。収束品質の深掘りは Phase 2 (Rao 照合) で実施。

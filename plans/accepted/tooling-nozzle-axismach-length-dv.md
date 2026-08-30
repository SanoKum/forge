# axis-Mach チェーン: 全長 dv 化 (`Lc_mode: from_length`)

## メタ

- **area**: `tooling / optimization`
- **status**: `done`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (axis-Mach チェーン §点の定義 / §全長を設計変数にする)
- **related_plans**:
  - 親: [`tooling-nozzle-axismach-chain.md`](tooling-nozzle-axismach-chain.md) (A0–A5)
  - 関連: [`tooling-nozzle-axismach-knot-law.md`](tooling-nozzle-axismach-knot-law.md) (knot 則 — 本 plan の逆算は knot にもそのまま適用)
  - 関連: [`tooling-nozzle-shortest-robust-axis.md`](tooling-nozzle-shortest-robust-axis.md) ($x_F$ 最小化探索 — 本 plan で「長さ指定」が直接 dv になった)
- **created**: `2026-08-31`
- **owner**: ユーザ指示 (2026-08-31)

## 1. 目的

axis-Mach チェーンの長さ系設計変数を、加速区間長 $L_c$ (= 軸 M が極大 $M_d$ に達する
$x_E$ までの長さ) から**スロート→物理出口の全長** $L_{\rm total} = x_F$ に替える。
$x_E$ (マッハ数極大点) は設計者が与える量ではなく、**終端特性線の計算から逆算される
従属量**にする。設備設計・機体搭載に効くのは全長であり、$x_F - x_E$ (E→F 一様化区間)
は $\approx r_F\sqrt{M_d^2-1}$ で $L_c$ にほぼ依存しない物理量なので、$L_c$ を dv に
すると「全長がいくつになるかは設計してみないと分からない」という帳簿の向きの悪さがあった。

## 2. スコープ

- **やる**: `design_chain` に `Lc_mode: from_length` を追加 (dv `L_total`、$r_t$ 単位・
  スロート $x=0$ 起点)。全 `axis_law` (quintic / knot / bspline / onepoint) 共通。
  既存 `explicit` / `max` は回帰対照として不変。
- **やらない**: 問題 YAML 資産の一括書き換え (既存 case は `explicit` のまま有効)。
  opt/MOO ドライバの dv 空間定義の変更 (dv 名を `L_total` にした YAML を書けばそのまま
  サロゲート探索に乗る — ドライバは dv dict に対して汎用)。CFD 検証 (壁生成は既存経路と
  bit 同一のため設計側テストで足りる)。

## 3. 関連 docs と前提

- [`methods/design/overview.md`](../../methods/design/overview.md) の「点の定義 (A/E/F)」:
  $x_F$ は終端特性線 (E 発 C⁺ と壁流線の交点, `moc_inverse.terminal_exit`) で決まる。
  直線 Mach line 近似 $x_F = x_E + r_F\sqrt{M_d^2-1}$ と 0.04% 一致 (実測)。
- $L_c$ スイープ ([nozzle-axismach-lc-sweep.md](../../notes/investigations/nozzle-axismach-lc-sweep.md)):
  $x_F(L_c)$ は単調増加・傾き ≈1 ($M_d$=4/R=2 で $L_c$ 6→9.9 に対し $x_F$ 19.1→22.9)。

## 4. 設計方針

- **逆問題の解き方**: $x_F(L_c) = L_{\rm total}$ を $L_c$ について解く。
  1. 解析初期推定 $L_c^{(0)} = L_{\rm total} - x_A - r_F\sqrt{M_d^2-1}$
     (終端 Mach line 近似、0.04% 精度)。
  2. ステップ $=-$残差 の固定点反復 $L_c^{(k+1)} = L_c^{(k)} - (x_F^{(k)} - L_{\rm total})$。
     写像の傾きが ≈1 なので縮小率 ~0.05、実測 1–2 パスで収束。
  3. 最良反復点 (|残差| 最小) を採用。改善が止まったら打ち切り (ノイズ床)。
- **secant を使わない理由 (実装時の発見)**: 離散 MOC の $x_F$ には解像度依存の
  **ノイズ床**がある (n_axis 500 で ~0.05 $r_t$ の階段 — 隣接 $L_c$ で局所傾き 0.2〜1.6
  に暴れる)。secant はこの局所勾配を信じて発振し 12 回で未収束だった (実測)。
  固定点反復は大域傾き ≈1 だけを仮定するのでノイズに頑健。
- **収束判定**: `geometry.L_total_tol` (既定 0.05 $r_t$ = n500 ノイズ床相当)。
  tol を詰めたいときは `n_axis_inv` を上げてから下げる。tol 未達で改善も停止した場合は
  明示エラー (ノイズ床の可能性としてその旨を案内)。
- **窓との関係**: 反復中の $L_c$ は単調ゲートの許容窓に数値マージン付きでクランプ。
  窓端に張り付いたまま残差が残る = 窓内で実現不能な `L_total` であり、窓と窓端での
  到達 $x_F$ を添えて明示エラー。
- **各パスは既存 `design_chain` の設計パスそのもの** (law 構築 → gates → 逆 MOC →
  terminal_exit)。壁 QA / 壁表現は収束後の最終パスにのみ適用 (従来と同一順序)。
  収束した $L_c$ を `Lc_mode: explicit` に渡すと **bit 同一**の壁 / $x_F$ になる (往復一致)。
- **knot 則との合成**: dv は (`L_total`, `M_K`)。$M_K$ は独立のまま (逆算は $L_c$ のみ)。
- **診断**: 返り値 / `prepare_info.json` に `Lc_mode` と `Lc_solve`
  (`L_total_target` / `xF_residual` / `n_design_evals` / `tol`) を出す。

## 5. 実装ステップ

1. `design_chain` の L_c 選択部を `_make_law` / `_checked_law` / `_run_inverse` に
   分解し、`from_length` の固定点反復を追加 (`runner_axismach.py`)。
2. `prepare` / `prepare_ns` の `prepare_info.json` に `Lc_mode` / `Lc_solve` を追記。
3. 統合テスト `design/tests/run_axismach_length_tests.py`。
4. `methods/design/overview.md` の該当節を更新。

## 6. 検証

`design/tests/run_axismach_length_tests.py` (2026-08-31 ALL OK):

- M4/R2 quintic, `L_total`=20: 残差 0.036 ≤ 0.05、設計パス 1 回、$L_c$=6.783 窓内。
- 往復一致: 逆算した $L_c$ を `explicit` で渡すと $x_F$ の差 0 (bit 同一)。
- `L_total`=30 (窓上限超): 「実現不能」ValueError (窓と窓端 $x_F$ を報告)。
- M6/R3 knot (`axis_dx0: 0.03`), `L_total`=100: 残差 0.034、設計パス 2 回。
- 既存回帰 `run_axismach_wall_tests.py` ALL PASS (幾何経路は不変)。

## 変更ログ

- 2026-08-31: 実装・テスト・docs 同期。secant → 固定点反復への置換 (ノイズ床の発見) を
  §4 に記録。status: done。

# ノズル設計 wall-driven チェーン (①風洞): 円弧廃止 + CFD 特性抽出 + wave-cancellation MOC

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`  <!-- Phase W1–W2 実装済み (2026-08-15)。W0 はユーザ判断待ち -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (現在仕様。「wall-driven チェーン」節)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md)
  - 兄弟: [`tooling-nozzle-phase3-windtunnel.md`](tooling-nozzle-phase3-windtunnel.md) (現行 B8 チェーン。**生産構成として維持し壊さない**)
  - 出典調査: [`notes/investigations/nozzle-throat-curvature-shape-representation-survey.md`](../../notes/investigations/nozzle-throat-curvature-shape-representation-survey.md)
  - 引き継ぎ: [`notes/sessions/2026-08-15-handover-wind-tunnel-nozzle-design.md`](../../notes/sessions/2026-08-15-handover-wind-tunnel-nozzle-design.md)
- **created**: `2026-08-15`
- **owner**: `sano`

## 1. 目的

ユーザ持ち込みの外部推奨方針 (ChatGPT 作成の「wall-driven 設計ツール試作仕様」、以下「原案」) を
引き継ぎ §2 の確定事実・§3 の文献結論に照らして評価し、採否と修正を確定して実装する。
完成時には「U→T→D 低自由度多項式壁 → CFD (遷音速+初期超音速) → 特性情報抽出 →
D 下流 wave-cancellation MOC → δ\* 補正 → 全域 B-spline → NS/RANS 検証」の
チェーンが現行 B8 チェーンと**並存**して動き、同一要求 (①風洞 $M_d=4$) で
出口一様性とこぶ振幅を直接比較できる状態を得る。

## 2. 原案の評価 (確定事実 §2 / 文献結論 §3 との照合)

### 2.1 判定サマリ: **骨格は採用** — 3 経路のうち (b)+(c) を同時に実装し、(a) を dv 化する構成

こぶ問題の確定帰結 (引き継ぎ §2) は「効く経路は (a) $R$ を上げる / (b) スロート下流円弧を
設計出力に置き換える (接合廃止) / (c) 遷音速アンカーの高次化、の 3 つだけ。J 下流の壁表現
変更は的外れ」。原案を分解すると:

| 原案の要素 | 3 経路との対応 | 判定 |
| --- | --- | --- |
| 凍結円弧を廃止し U→D を低自由度多項式で設計 ($R_t$ はスロート曲率**条件**として残る) | **(b)** の変種。超音速域内の「κ=const 円弧 → 設計壁」接合 J が消滅。T→D 多項式の κ は $1/R_t$→0 へ**連続に減衰** (κ' 跳びの接合が超音速域から消える) | 採用 |
| 遷音速+初期超音速を CFD で解き特性情報を抽出、設計モデルにしない | **(c) の完成形**。B6 (CFD アンカー化でこぶ 58% 減) の実測と整合。Sauer の「$M''$ が構造的に存在しない」問題 (調査 §9-2) が消滅 | 採用 |
| D 下流は wave cancellation で壁を**従属出力**として生成 | CONTUR の「壁は MOC の従属出力」思想 (調査 §2.1) そのもの。データ線が D からの C⁻ 特性線になるので「縦線 starting line は古典に対応物がない」問題 (調査 §9-5) も解消 | 採用 |
| $M_c(x)$ を設計入力から診断量へ降格、出口 M は $\theta_D$ の root find | B10 が暴いたジレンマ (始端 C2 整合 vs 出口一様性、自由 CP 3 個で両立不能) が**構造的に消滅** — 整合させるべき target 曲線が存在しない | 採用 |
| $R_t$ が明示的な設計変数 | **(a)**。outer loop / MOO で $R$ を直接動かせる | 採用 |

### 2.2 正直な限界 (原案が言っていないこと)

- **軸こぶが「消える」保証はない**。到達不能域 ($x<x_{reach}$) の軸 M は U→D の壁形状と
  $R_t$ だけで決まる物理であり、新チェーンでも下流壁からは届かない (確定事実は不変)。
  新チェーンの利点は (i) スロート域が作る波が**実 CFD 場として正確に**下流設計に入り、
  相殺壁が実場基準になるので試験部への漏れが構造的に減る、(ii) $R_t$ と T→D の κ 分布が
  直接 dv なので、こぶ振幅そのものを outer loop で下げに行ける、の 2 点。
- **仮説 (検証対象)**: T→D 多項式は κ' 連続なので「同じ $R_t$ でも円弧+J 構成より波源が
  弱い」可能性がある。ただし κ' 不連続 (G3 破れ) の効果を G2 破れと分離定量した文献は
  存在しない (調査 §2.4 末尾) ので、これは実測で判定する (Phase W9 の A/B)。
- **Korte の相場観** (調査 §2.5): 遷音速モデル差起因の 0.02 級うねりは古典チェーンの実測
  相場。新チェーンはまさにその誤差源を CFD 実場で置換するので、0.02 を下回れるかが
  成否の定量指標になる。

### 2.3 原案から修正する点

1. **新規ディレクトリ構成 (原案 §26) は不採用**。`design/forge_design` の既存構成に統合する
   (原案自身が「既存 repository に合わせて変更してよい」と許容)。再利用資産:
   `_hermite_quintic` (wall.py) / `moc_kernel` (単位プロセス) / `cminus_cfd` (C⁻ 追跡) /
   `cfd_anchor` (局所多項式フィット抽出) / `mesh2d`+`runner_wt` (メッシュ・node Euler 評価) /
   `deltastar` (δ\* 相関) / `metrics`。
2. **判別実験 W0 (現行チェーン R=5 A/B、調査 §8) を並走枠として残す**。新チェーンの既定
   $R_t$ と「こぶの R 依存 vs 構成依存の分離」という物理情報を最安 (config+1 パス) で得る。
   新チェーン実装と独立に実行可能。**実行するかはユーザ判断** (引き継ぎ §8 の指示通り)。
   run 番号は `run_0039` 以降。
3. **U→T は既存 Bell–Mehta 5 次収縮と同型** (両端 6 条件の quintic Hermite)。差分は
   「円弧に C2 接続」でなく「throat で直接 $r''=1/R_t$」を課すことだけ。原案 §7 の単調条件
   $\mu=L_U^2/(R_t\Delta r)\le20$ は検証済み (§4.1) だが大収縮比では緩い制約で、実用制約は
   $\max|d\kappa/ds|$ と入口境界層 (Bell–Mehta の趣旨) の側にある。
4. **原案 §10 (直管→一定傾斜) は geometry option として実装・保存** ($\lambda\in[5/3,5/2]$、
   $\lambda=2$ で 4 次退化)。ただし①風洞の主経路では使わない。
5. **実在気体 MOC の抽象化 (原案 §14) は初期実装では最小限**。`moc_kernel` は γ=const の
   PM 関数ベースなので、ThermoModel 分離は TP 拡張の必要が生じた時に行う (原案も許容)。
   γ を関数引数に保つ現行流儀は維持する。
6. **D 下流の CFD ドメイン延長壁は θ_D 直円錐**。D からの C⁻ より上流の場は延長壁に
   依存しない (超音速の依存域) ため、データ線抽出には十分。
7. **BL 補正の λ 指定より tolerance 制約優先 (原案 §21) に同意** (MVC 汎関数、調査 §3.2/3.4
   と整合)。physical exit slope=0 を課さない (§22) も同意。ただし Phase W6–W8 は後段。

### 2.4 原案の数式の独立検証 (2026-08-15 実施、実装前に完了)

乱数 200 試行 + 境界 ±0.1% すれすれ試験で数値検証済み (`design/tests/run_walldriven_tests.py`
に恒久化):

- U→T 勾配閉形式 $\frac{dr}{dx}=\frac{\Delta r}{2L_U}\xi^2(\xi-1)[5\mu\xi-3\mu-60\xi+60]$ は
  quintic Hermite の厳密解と一致 (相対 1e-9 以下)。単調収縮 ⇔ $\mu\le20$ は**必要十分**
  (μ=20×0.999 単調 / ×1.001 非単調)。
- T→D の $r_D=r_t+\frac{L_D^2}{12R_t}+\frac{L_D}{2}\tan\theta_D$ は適応求積と機械精度一致
  (3e-16)。$r''$ 閉形式一致。全区間 $r''\ge0$ ⇔ $L_D\le3R_t\tan\theta_D$ は**必要十分**。
- §10 の $\lambda$ 窓 $[5/3,5/2]$ と $\lambda=2$ の 4 次退化も検証済み。

## 3. スコープ

- **やる**: Phase W1–W9 (§6)。①軸対称風洞のみ。ideal gas MOC。既存 B8 チェーンとの比較評価。
- **やらない**: 実在気体 MOC (抽象化のみ準備) / ②–⑤機種への展開 / 現行 B8 チェーンの改廃
  (比較で新チェーンが勝つまで生産構成は B8) / ソルバ側変更。

## 4. 設計方針 (差分のみ — 数式は methods 側)

幾何の定義式・検証条件は [`methods/design/overview.md`](../../methods/design/overview.md) の
「wall-driven チェーン」節に記載。実装上の判断:

- **無次元系は既存踏襲** ($r_t=1$、スロート $x=0$、U は $x=-L_U$)。
- **`geometry/wall_walldriven.py`** に `UpstreamThroatPoly` (U→T) / `ThroatExpansionPoly`
  (T→D) / `WallDrivenThroatRegion` (合成 + validator + 診断) を新設。$r,r',r'',r'''$ は
  全て解析式で持ち、κ と $d\kappa/ds$ も解析評価 (数値微分しない)。
- **validator は「解析条件 + 数値サンプル検査」の二重化**: 解析条件 ($\mu\le20$、
  $L_D\le3R_t\tan\theta_D$、$r>0$) と、密サンプルでの単調性・変曲・$\max|d\kappa/ds|$ を
  両方返す (解析条件は表現を差し替えたら失効するため)。
- **T 点の κ' 跳びは診断値として常時出力する** (両区分とも $\kappa(T)=1/R_t$ で κ 連続だが
  κ' は一般に不連続。跳び量は幾何診断 CSV に載せ、遷音速域内の跳びとして許容可否を
  CFD で確認する — 調査 §2.4 の教訓により「滑らかだから無害」とは推定しない)。
- **特性抽出は node 場 + 局所多項式フィット** (cell の $M''$ 雑音支配・node 124× 改善の
  教訓を踏襲)。データ線上の $P_0$ 変動を等エントロピー仮定の妥当性診断として記録する。
- **cancellation MOC** は `moc_kernel` の単位プロセスを再利用し、データ線 = D からの C⁻、
  壁点は「入射波を反射させない壁角」で逐次決定。回転 MOC は $P_0$ 診断が悪い場合の
  フォールバックとして温存。

## 5. 関連 docs と前提

- 現在仕様: `methods/design/overview.md` (実装と同期して更新)
- 確定事実・撤回履歴: `tooling-nozzle-phase3-windtunnel.md` §9.2 (B5–B10)
- AGENTS 必須ルール (新連番 run / メッシュ VERDICT / NaN・収束・準定常チェック /
  README run 一覧同期 / 検証後 commit) を全 Phase に適用。

## 6. 実装ステップ

| Phase | 内容 | 状態 |
| --- | --- | --- |
| **W0** | 判別実験: 現行 B8 チェーンで R=5 A/B (調査 §8 の計画のまま、run_0039〜)。こぶの R 依存を測る。**ユーザ確認後に実行** | 未着手 (ユーザ判断待ち) |
| **W1** | 幾何生成器: `wall_walldriven.py` (U→T / T→D / 合成 / validator / 解析 κ・dκ/ds / §10 option) | **済 (2026-08-15)** |
| **W2** | 解析テスト + 可視化: 端点条件・境界すれすれ・違反ケースの検出、$r,\theta,\kappa,d\kappa/ds$ 図 | **済 (2026-08-15)** |
| W3 | メッシュ/CFD 接続: `mesh2d` を U→D+θ_D 円錐延長ドメインに対応、`runner_wt` (node Euler) で評価。品質 VERDICT | 未着手 |
| W4 | 特性情報抽出 + cancellation MOC: D からの C⁻ データ線抽出 (`cminus_cfd` 拡張)、`moc_kernel` 単位プロセスで非反射壁マーチ、$P_0$ 診断 | 未着手 |
| W5 | outer loop: $M_e(\theta_D)$ root find (必要なら $R_t, L_D$ も)。$M_c(x)$ ・こぶ振幅を診断出力 | 未着手 |
| W6 | δ\* 補正: `deltastar` 相関 + smoothing spline + 法線オフセット | 未着手 |
| W7 | 全域 B-spline fitting: tolerance 制約 ($\max d\le\varepsilon_{geom}$) 下で $\int(d\kappa/ds)^2ds$ 最小化。throat 等の hard constraint | 未着手 |
| W8 | STEP/CSV 出力 (親計画 5a と合流) | 未着手 |
| W9 | NS/RANS 検証 + **B8 チェーンとの同条件比較** (ε_M / ε_θ / こぶ振幅 / 同一 R での円弧 vs 多項式 A/B) | 未着手 |

## 7. 検証

- W1–W2: 解析 unit test (`design/tests/run_walldriven_tests.py`) — 端点条件 (1e-12)、
  閉形式一致、境界 ±0.1% の判別、違反ケースで非単調/変曲が実際に発生すること、
  解析 κ・dκ/ds と数値微分の一致。
- W3 以降: 各 Phase で AGENTS 必須ゲート (メッシュ VERDICT / `check_convergence.py` /
  `check_quasisteady.py`) を通し、run は `case/41.wind_tunnel_design` の README 表に登録。
- 最終判定 (W9): 同一要求で B8 比較 — 出口 ε_M/ε_θ が要求内、こぶ振幅が B8 (0.019–0.023)
  以下、Korte 相場 0.02 を下回るかを報告。

## 8. 影響範囲

- 追加のみ: `design/forge_design/geometry/wall_walldriven.py`、`design/tests/run_walldriven_tests.py`。
- 既存チェーン (B8 生産構成・`wall_modef`・`moc_inverse`・帰還ループ) は無変更。
- W3 で `mesh2d`/`runner_wt` に option 追加 (既存経路は不変に保つ)。

## 9. 完了条件

1. W9 の同条件比較で新チェーンの採否が定量判定されている (勝てば生産構成の移行 plan を
   別途起票、負ければ本 plan を archived へ移し敗因を記録)。
2. `methods/design/overview.md` が最終仕様と同期している。

## 変更ログ

- 2026-08-15: 起票。原案 (ChatGPT 仕様) を引き継ぎ §2/§3 と照合して採否確定 (§2)。
  原案の設計式を独立数値検証 (§2.4、全て正)。W1–W2 実装・テスト 34 件 PASS・可視化完了。

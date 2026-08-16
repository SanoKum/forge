# 軸 Mach 則 A の最短ロバスト軸長探索

## メタ

- **area**: `tooling / optimization`
- **status**: `done`
- **related_docs**:
  - `methods/design/overview.md` (§ knot 付き区分 C² 軸 Mach 則)
- **related_plans**:
  - [tooling-nozzle-axismach-knot-law.md](../accepted/tooling-nozzle-axismach-knot-law.md)
  - [tooling-nozzle-topology-crossing-speedup.md](../accepted/tooling-nozzle-topology-crossing-speedup.md)
- **created**: `2026-08-17`
- **owner**: `Codex`

## 1. 目的

既存の $M_K$ 感度は各固定 $L_c$ 内で壁角を最小化しており、「可能な限り短い軸長」という目的を
直接解いていなかった。必須品質と工学的余裕を制約にし、閉包後端 $x_F$ が最短となる
$(L_c,M_K)$ を探索・解像度確認する。

## 2. スコープ

- **やる**: M6/R3 で $35<L_c<40$ の境界探索、$M_K$ 同時最適化、600点探索、境界候補のみ
  1200点確認、JSON/CSV成果物、既存軸則・MOCテスト。
- **やらない**: 2400点の一律評価、CFD、余裕閾値そのものの較正、軸則 A の数式変更。

## 3. 関連 docs と前提

現在仕様と評価順序は実装前に `methods/design/overview.md` へ反映した。M6/R3 では
$x_F-x_E$ がほぼ固定で、$x_F$ 最小化は $L_c$ 最小化と同値である。既存の交差検査高速化と
候補process並列を利用する。

## 4. 設計方針

候補は hard gate と $M''$ 単峰条件を満たし、さらに
$g_\mu=\min(\mu-\theta)-1^\circ\ge0$、
$g_t=\min|\sin\angle|-0.02\ge0$ を満たすときロバスト合格とする。目的は辞書順に
(1) $x_F$ 最小、(2) 同一長さで $\theta_{max}$ 最小、(3) 余裕最大とする。
600点で長さ境界を絞り、1200点の再評価で最終採否を決める。1200点で不合格なら長さを増やして
次候補を確認する。

## 5. 実装ステップ

1. `case/42.isobutane_wt/optimize_axislaw_A_shortest.py`: 粗探索と境界細分、候補選択、成果物出力。
2. 600点で $L_c=35$〜40 と $M_K$ を探索。
3. 最短境界候補を1200点で確認し、文書とcase索引へ記録。

## 6. 検証

- **単体**: 既存 `run_axislaw_tests.py`、`run_inverse_tests.py`、`run_axismach_wall_tests.py`、
  `run_moc_diagnostics_tests.py`。
- **検証ケース**: `case/42.isobutane_wt/` の設計のみ。CFD runは作らない。
- **判定基準**: 1200点で全制約合格。600→1200で最短候補の判定が逆転した場合は、隣接する長い
  候補へ進み再確認する。

探索実績: 600点では$L_c=39.03,M_K=2.76$がmargin 1.00027°で最短だった。1200点では
$L_c=39.03/39.04$が0.99748/0.99922°で不合格、39.05が1.00095°で初めて合格した。同長さで
合格する$M_K$は2.76/2.77/2.78、そのうち二次目的の壁角が最小の2.76を採択した。
最小sin角0.04177、内部flip/交差0。全探索87.76 s、最大RSS 308 MiB。

## 7. 影響範囲

- 設計探索スクリプトとJSON/CSV成果物。
- `methods/design/overview.md`、`case/42.isobutane_wt/README.md`、`plans/README.md`。
- 生産既定値やCFD入力は本計画だけでは変更しない。

## 8. 完了条件

- [x] 関連 `methods/design/overview.md` の現在仕様を更新済み
- [x] 実装・検証完了
- [x] 本計画を `done` にし、§9に変更ログを記載
- [x] `plans/accepted/` へ移動
- [x] `plans/README.md` の一覧を同期

## 9. 変更ログ

- `2026-08-17` — 起票。壁角最小から、品質制約下の実全長最小へ評価順序を変更。
- `2026-08-17` — 600点の3段階探索と1200点境界確認を完了。最短合格を
  $L_c=39.05,M_K=2.76$ と確定。関連6 suite ALL PASS、acceptedへ移動。
- `2026-08-17` — 結論・数式・探索ファネル・境界/M_K感度図・トポロジ・不確実性をまとめた
  自己完結型結果ページ `report_axislaw_A_shortest.html` と再生成スクリプトを追加。
- `2026-08-17` — 最良候補を1200点で再構築し、A–K–E–Fの軸中心Mach図、壁内20,710点の逆MOC
  Machコンター、図示用圧縮NPZを追加。結果ページへbase64埋込みし、一次結果パスも冒頭に明記。

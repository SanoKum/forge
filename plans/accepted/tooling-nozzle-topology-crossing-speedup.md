# 逆 MOC 特性線交差検査の候補線分 sweep 高速化

## メタ

- **area**: `tooling / optimization`
- **status**: `done`
- **related_docs**:
  - `methods/design/overview.md` (§ 軸則の滑らかさ比較 → 非隣接特性線の交差検査)
- **related_plans**:
  - [tooling-nozzle-axislaw-smoothness.md](../accepted/tooling-nozzle-axislaw-smoothness.md)
  - [tooling-nozzle-axislaw-onepoint.md](../accepted/tooling-nozzle-axislaw-onepoint.md)
- **created**: `2026-08-16`
- **owner**: `Codex`

## 1. 目的

`characteristic_topology_diagnostics` の非隣接特性線交差検査が、折れ線対ごとに全線分対の
外積行列を生成し、`n_axis=1200` の設計感度掃引を律速している。交差判定の意味を変えずに
候補線分対を事前削減し、高解像度の多点掃引を実用時間で実行可能にする。

## 2. スコープ

- **やる**: x 単調折れ線向け AABB sweep、x 非単調折れ線向け汎用チャンク fallback、旧実装との
  判定一致テスト、M6/R3/$L_c=45$ の n_axis=600/1200/2400 ベンチマーク、独立な $M_K$ 候補の
  解像度別プロセス並列評価。
- **やらない**: サンプルする特性線本数・非隣接の定義・交差の外積判定式・hard gate の変更、
  Numba/Shapely 等の新規依存追加、MOC 本体の変更。

## 3. 関連 docs と前提

現在仕様と計算量は `methods/design/overview.md` に先行反映した。既存の
`_segments_cross_matrix` を参照実装として保持し、ランダム・接触・非単調・実 MOC 網で新旧結果を
照合する。

## 4. 設計方針

折れ線 q が x 単調なら、各線分の x 区間下端・上端も順序付けられる。折れ線 p の各線分について
`searchsorted` で x 区間が重なる q 線分の添字範囲を求め、r 区間の AABB 重複でさらに絞り、残った
線分対だけを従来と同じ4外積の符号判定へ渡す。一方だけがx単調なら引数を交換する。両方が非単調な
場合は、pを固定長チャンクに分けたAABB行列で候補を抽出し、一時メモリを制限する。

## 5. 実装ステップ

1. `design/forge_design/geometry/moc_diagnostics.py`: 線分メタデータ、単調 sweep、汎用 fallback。
2. `design/tests/run_moc_diagnostics_tests.py`: 参照全対行列との一致と実 MOC 網回帰。
3. `case/42.isobutane_wt/sweep_axislaw_A_MK.py`: 600/1200/2400 点をメモリ量に応じて 6/4/2
   process で並列評価し、高速化後に中断した $M_K$ 感度掃引を再開。

## 6. 検証

- **単体**: 手作り交差/非交差/接触、固定seedランダム折れ線、x非単調折れ線で旧全対行列と一致。
- **実網**: M6/R3/$L_c=45$, n_axis=600/1200/2400 で新旧の交差数が一致。
- **性能**: n_axis=1200, n_crossing_sample=40 の交差検査を旧 5.41 s から 1.0 s 未満へ短縮し、
  `characteristic_topology_diagnostics` 全体を 1.5 s 未満とする。判定結果は不変。12 CPU / 7.7 GiB
  環境では、2400 点1候補の最大 RSS 0.83 GiB を根拠に並列数を制限し、swap を使わない。

実績:

- 手作り6ケース + 固定seedランダム1000組で旧全対行列との不一致0。
- M6/R3/$L_c=45$ 実網で crossing は600/1200/2400すべて旧0＝新0。
- 交差検査: 600点 1.1096→0.0149 s (74.6倍)、1200点 3.6879→0.0248 s (148.9倍)、
  2400点 23.7397→0.0523 s (453.9倍)。
- 1200点1候補全体: 6.62→1.15 s、topology診断 5.85→0.45 s。
- $M_K$ 感度156点+代表1200/2400を6/4/2 processで **32.41 s**。swapなし。
- `run_moc_diagnostics_tests.py` と関連5 suite (`axislaw`, `onepoint`, `bspline`, `inverse`,
  `axismach_wall`) は全て ALL PASS。

## 7. 影響範囲

- `moc_diagnostics.py` の診断性能のみ。MOC解・壁・CFD入力は不変。
- `methods/design/overview.md`、本plan、`plans/README.md`。

## 8. 完了条件

- [x] 関連 `methods/design/overview.md` の現在仕様を更新済み
- [x] 実装・検証完了
- [x] 本計画を `done` にし、`plans/accepted/` へ移動
- [x] `plans/README.md` の一覧を同期

## 9. 変更ログ

- `2026-08-16` — 起票。cProfile で n_axis=1200 の全 6.62 s 中、交差検査 5.41 s、MOC 本体
  0.75 s と特定。全線分対行列の候補削減を採用。
- `2026-08-17` — AABB sweep + 非単調chunk fallbackを実装。実網で最大453.9倍、旧判定と完全一致。
  解像度別process並列も追加し $M_K$ 感度を32.41 sで完走。関連6 suite ALL PASS、acceptedへ移動。

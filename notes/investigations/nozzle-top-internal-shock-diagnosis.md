# TOP ベルノズル内部衝撃波と軸マッハ過大の診断 — forge は無罪 (SU2 2×2 で確証)

- **date**: 2026-08-13 (夜間自律調査)
- **area**: `verification / nozzle`
- **発端**: case/40 の Rao 80% ベル (rao80: θn33.8°/θe12.1°/L5.97, ε=9, Pt4MPa/Tt1500K/Pb20kPa)
  のマッハコンタに「不連続」があるとの指摘。別途の画像ベース分析 (別 Claude インスタンス) は
  「軸上出口 M=4.7 は 1D 値 3.81 の +22% で物理的に不可能 → 実装バグの徴候」と判定した。
  本調査はこれを実データ + 独立ソルバで検証した。
- **関連**: [tooling-nozzle-moo-loop](../../plans/accepted/tooling-nozzle-moo-loop.md) (Phase 2 評価器),
  [procedures/su2-cross-check.md](../../procedures/su2-cross-check.md) (比較手順),
  case/40 README run 一覧 (run_0075–0086)

## 結論

**forge のコーディングは無罪。観測された構造は 2 つとも実在の物理である。**

1. **内部斜め衝撃波** — スロート下流円弧→放物線接合部 ($P_a$, 曲率不連続) から発する圧縮波が
   下流で合流して成長する、TOP (Rao 放物線近似) ノズルの**教科書的な内部衝撃波**。
   トレース (各 x 断面の最大 Mach 降下位置) は $P_a$ から立ち上がり r≈13–15mm を漂い、
   出口に向かって軸側へ収束する。強度は下流ほど成長 (ΔM −0.1 → −0.4〜0.5)。
2. **軸上 Mach 過大 (出口 4.6–4.7 vs 1D 3.81)** — 短尺 (80% 長)・大初期膨張角 (θn=33.8°) の
   kernel 過膨張。軸は円弧膨張波を全部受けてから壁转向の圧縮 (=内部衝撃波) が届くまで
   膨張し続けるため、1D 値を大きく超える。**準 1D 論法を強い 2D 流に適用した「+22% は不可能」
   という判定が誤り** (質量保存は面平均で成立: ṁ/ṁ_1D = 0.986)。

## 証拠マトリクス (全 run は case/40, 同一メッシュ 221×65)

| # | 構成 | 軸上出口 M | 内部ジャンプ max ΔM | 判定 |
| --- | --- | --- | --- | --- |
| rao80_p1 | forge node RANS-SST (基準) | 4.70 | −0.43 | 発端 |
| run_0077 | forge, **Rd=0.382 (Huzel–Huang 正規値)** | 4.73 | −0.49 | Rd 0.4/0.382 の差は無関係 |
| run_0078 | **SU2 v8.5 RANS-SST, 同一メッシュ・同一 BC** | **4.63** | **−0.46** | **コード非依存を確証** |
| run_0079 | forge **cell** RANS | 4.72 | −0.44 | node 離散化固有でない |
| run_0085 | forge **Euler** (cell, slip, visc=0) | 4.72 | −0.42 | 乱流/壁関数/粘性は無関係 |
| run_0086 | **SU2 Euler** (rms[Rho] 10⁻¹¹ 深収束) | **4.62** | **−0.42** | 非粘性 2×2 完成 |
| run_0081/0082 | forge, Rd=0.8 / 1.5 | 4.56 / 4.35 | −0.40 / −0.40 | **Rd↑ で過膨張が緩和** = 円弧曲率起因と整合 |
| run_0083/0084 | forge, Ru=0.75 / 3.0 | 4.70 / 4.69 | −0.43 / −0.51 | **上流円弧は無関係** (理論通り) |
| run_0080 | forge, axisymMethod=1 | (4.13) | −0.44 | 陰解法深収束未達 (既知) のため証拠から除外 |

- forge↔SU2 の軸上 M 差は RANS +1.5% / Euler +2.3%、衝撃波位置差 0.3–0.7mm。
- y+1 low-Re 検証は不要と判断 — Euler (壁関数なし) が構造を完全再現しており、壁処理は
  無罪がより強い形で示されたため。
- 図: `case/40.nozzle_design_tool/run_0078_su2_rao80/cross_compare_2x2.png` (2×2 コンタ +
  衝撃波トレース重ね)・`cross_compare_profiles.png` (軸プロファイル 4 者一致 + 出口半径
  プロファイルの二層構造)。再生成スクリプト: 本ノート末尾参照。

## 副次知見

- **可視化の罠**: 散布セル中心の Delaunay 三角形分割 (`tricontourf`) は開いた出口境界近傍に
  縞状アーティファクトを作る。構造格子への `griddata` 線形補間 + 解析壁マスクで根治
  (`case/40.nozzle_design_tool/compare_shapes_mach.py`)。
- **画像ベースの外部診断は準 1D 論法に注意**: 「軸 M > 1D 値 → 非物理」は 2D 性の強い
  短尺ベルでは成立しない。軸過膨張は TIC/TOP で既知 (センターライン M オーバーシュート)。
- Ru は η には効く (0.958@Ru0.75 → 0.969@Ru3.0) が内部衝撃波構造には効かない —
  Ru は遷音速域/収縮損失側の変数。
- SU2 Euler は本ケースで完全収束する (8.4k iter で rms[Rho] 10⁻¹¹) — 今後の非粘性
  クロスチェックの基準として安価。

## 再現手順

- forge 側バッチ: `night_batch.py` (セッション scratchpad; 要旨 = runner で prepare 後、
  Euler は `viscMethod: 0, visc: 0.0, thermCond: 0.0` + `kind: slip` に post-edit、全 run
  soft 段 3000 step → 本段 12000 step)。診断は各 run の `diag.json` / `diag.txt`。
- SU2 側: forge の `nozzle.msh` を `gmsh -2 -format su2` で変換 (マーカー名維持)、
  cfg は `run_0078_su2_rao80/bell_sst.cfg` / `run_0086_su2_euler/bell_euler.cfg`。
- 比較図: `case/40.nozzle_design_tool/cross_compare_rao80.py`
  (SU2 バイナリ restart の直接パース込み — ヘッダ 20B + 33B×nf 名前 + f64 行列)。

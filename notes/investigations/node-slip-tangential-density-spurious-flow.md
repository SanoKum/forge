# node slip 境界の接線密度勾配スプリアス流 (未修正バグ)

- **発見**: 2026-07-20、等温壁検証 (純伝導ボックス) の副産物
- **証拠 run**: `case/24.laminar_channel_bl/run_isoT_condL_node` (発生) / `run_isoT_condP_node` (slip→periodic で消失)
- **状態**: 未修正。切り分けのみ完了

## 症状

静止流体 + 線形温度場 (下壁 350K / 上壁 300K、ρ=P/RT が境界接線方向に変化) の厳密解を
IC に張った純伝導ボックスで、**slip 側面境界の列に市松状 (符号交互) の接線速度 max 0.47 m/s が
発生し、誤った定常解として定在する** (残差プラトー)。Pe~250 の移流で温度場が最大 2.2 K 汚染される。

- 速度の分布: slip 境界上のノード列に集中 (max\|U\| 0.475)、内部 0.10、壁 (Dirichlet) は厳密 0。
  左右の slip 境界でミラー対称、境界に沿って 1 ノードごとに符号が反転する市松。
- **slip → periodic に変えるだけで 0.47 → 3.6e-4 m/s** (≈1300 分の 1)・T ドリフト 2.2 → 0.11 K。
  等温壁・内部スキームは無罪で slip 閉包が犯人。
- cell モードは同条件で max\|U\| ~1e-4 (無害)。**node 固有**。

## 含意

- node で slip 境界と「境界接線方向の密度 (温度) 勾配」が同居するケースは要注意:
  例) ノズル遠方場/対称面を slip にした温度成層のある流れ、`case/29` 系のプルーム境界。
- 一様密度なら顕在化しない (これまでの node slip 検証 (Sod slip 等) が素通りした理由)。

## 未調査 (次の一手)

- `slip_d` の node 閉包 (bvar ro/roe の構成) と境界半割面の運動量流束のどこで接線運動量が
  湧くか (ρ の接線差分と法線反転の交差項が疑わしい)。
- ghostless 化 ([discretization-node-boundary-ghostless plan](../../plans/active/discretization-node-boundary-ghostless.md))
  で slip を作り直す際に本ケース (`run_isoT_condL_node` 構成) を回帰テストに使うこと
  (厳密解シード + 「静止が保てるか」は 10k step で判定できる)。

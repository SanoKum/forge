# TP 擬似種 + H₂O 独立種 (split_h2o) で燃焼ガスノズルの H₂O 凝縮を解く

## メタ

- **area**: `tooling / condensation / thermophysics`
- **status**: `done`
- **related_docs**:
  - `methods/design/overview.md` § ガスモデル (split_h2o の定義)
  - `methods/condensation.md` §7b (carrier を擬似種に畳む運用)
  - `design/forge_design/gas/semiperfect.py` (`mixture_pseudo_species_split`)
  - `design/forge_design/evaluate/ic.py` (`paste_isentropic_ic(species_Y=)`, `cpg_field_to_tp`)
  - `design/forge_design/evaluate/runner_axismach.py` (`evaluate.tp_species`, `evaluate.condensation`)
- **related_plans**:
  - 前提: [tooling-nozzle-semiperfect-gas.md](../accepted/tooling-nozzle-semiperfect-gas.md) (単一擬似種 TP、thermoHrefTemp)
  - 参照: [tooling-nozzle-axismach-knot-law.md](../accepted/tooling-nozzle-axismach-knot-law.md) (M4.2 ベスト run_0035)
- **created**: `2026-08-17`
- **owner**: `Claude (autonomous)`

## 1. 目的

イソブタン燃焼ガス風洞 (M4.2 ベスト `run_0035_ib_R2_LU6_Lc8_tp`) を凝縮計算へ引き継げるようにする。
現行 TP CFD は全組成を単一擬似種 `MIX` に畳むため凝縮種 (`condGasSpecies`) を指せない。
**H₂O 以外を擬似種 `MIXDRY`、H₂O を独立種**とする 2 種 TP で dry と同じ軸 M が出ることを確認し、
既存の H₂O 凝縮モデル (`condModel 1`, Kantrowitz+Hertz–Knudsen) を有効にして Wyslouzil ノズルで
既存結果を再現、イソブタン M4.2 (H₂O 質量分率 5 %) で凝縮計算を行う。液相は既存モーメント方程式
($g$=液相質量分率 β, $Q_0..Q_2$) で解く。

## 2. スコープ

- **やる**: `mixture_pseudo_species_split`、IC の `roY{s}` 出力と CPG→TP 場変換ツール、runner の
  `tp_species: split_h2o` / `condensation` 通し、Wyslouzil fig3 2D SST 再現 (TP split vs 既存 CPG
  `run_0050`)、イソブタン M4.2 H₂O 5 % の dry (pseudo vs split 一致) と凝縮 ON。
- **やらない**: forge ソルバ本体の凝縮/多成分実装変更 (既存機能の組合せで済む範囲)。二温度
  (`condTwoTemp`) や成長則の再検証。3D。

## 3. 前提

- forge: `thermalMethod 2` 多成分 TP (`species`, `speciesDBFile` で擬似種を追加/上書き、`Tlo` クランプ)、
  凝縮 `condensation/condModel 1/condGasSpecies` は実装済 (case/16 は CPG+species で検証済、TP 分岐は
  unit 検証のみ)。IC は `VALUE/roY{s}` を読む (無ければ Y0=1)。
- 既知: 多成分 TP × 陰解法は cfl ≤2 (thermoHrefTemp で緩和)、TP `inlet_Pressure` は rf=0 修正済。

## 4. 設計方針

- **DB**: `species_db.yaml` に `MIXDRY` (H₂O 以外の NASA-9 質量分率線形混合、MW_mix、LJ は質量平均) と
  `H2O` (内蔵と同じ係数を明示的に書き出す = 自己完結) の 2 エントリ。`species: [MIXDRY, H2O]`、
  入口 `Y0=1−Y_H2O`, `Y1=Y_H2O`、`condGasSpecies: 1`。
- **設計**: MOC は全組成 frozen semi-perfect (変更なし)。IC は同じガスで roe (datum 298.15) を作り
  `roY0/roY1` を追加 (組成一様)。
- **検証 1 (Wyslouzil)**: 既存 `run_0050_fig3_2d_sst_kwhk` (CPG, `[N2,H2O]`, Kw+HK) と同メッシュ・同 BC で、
  `MIXDRY`=N₂ のみの擬似種 + H₂O、TP、`run_0048` (dry SST 発達場) を CPG→TP 変換した IC から Kw+HK を
  回し、中心線 p/p0 と g を比較する。TP と CPG の差は N₂ cp が 100–290 K でほぼ一定なので小さい見込み。
- **検証 2 (イソブタン M4.2)**: 組成を H₂O 5 % (残りは元の dry 比率で 95 %) にした問題定義で
  (a) pseudo 単一種、(b) split dry、(c) split + 凝縮 (Kw+HK) を run し、(a)(b) の軸 M 一致と (c) の
  g・p・T 変化を報告する。

## 5. 実装ステップ

1. `gas/semiperfect.py`: `mixture_pseudo_species_split(Y, keep=("H2O",))`。
2. `evaluate/ic.py`: `paste_isentropic_ic(..., species_Y=None)` で `roY{s}`、`cpg_field_to_tp(...)`。
3. `evaluate/runner_axismach.py`: `_apply_gas_to_config` の split 分岐、bcond に `Y0/Y1`、
   `evaluate.condensation` ブロックの通し。
4. `case/16.nozzle_wys/run_tp_split_wys.py` (Wyslouzil 検証投入)、`case/42.isobutane_wt/problem_ib_R2_LU6_Lc8_*h2o5*.yaml`。
5. `design/tests/run_gas_tests.py` に split の等価性テスト追加。

## 6. 検証

- **単体**: `run_gas_tests.py` に split の等価性 (混合 cp/h が全種混合と 1e-9/1e-6 で一致、H₂O 係数が内蔵と同一) を追加、
  ALL PASS。他 6 スイート回帰なし。
- **Wyslouzil (case/16, H₂O 1.095 %)**: `run_0170` (cell, 既存 SST メッシュ, TP split, Kw+HK) は既存 CPG `run_0050` を再現
  (g_max 0.0106 vs 0.0107、onset 2.45 vs 2.35 cm、中心線 p/p0 差 ≤3.1 %)。**node** は平面 2D 一様メッシュ・非粘性で
  cell と一致 (`run_0184` cell CPG / `run_0185` node CPG / `run_0186` node TP split: g 0.0109/0.0109/0.0108、onset 1.75/1.75/1.85 cm、
  p/p0 差 ≤2 %)。node rms_ro 3e-9。
- **node の限界 (申し送り)**: 平面・壁クラスタ (y₁≈0.5–5 µm, AR 728) の no-slip メッシュでは node が
  壁ノード T > Tt (出口コーナー起点) → 衝撃列遡上 → unstart (SST/laminar、1 次/2 次、Dirichlet 0/1、出口 Pt=Ps、
  1D 等エントロピー IC いずれも; `run_0171/0179/0181/0183`)。slip なら健全 (`run_0182`)。cell の同メッシュ SST は問題なし。
  → node の粘性壁 (平面, 極薄壁セル) の問題としてソルバ側へ申し送り。
- **イソブタン M4.2 (case/42, H₂O 5 %, node Euler TP)**: `run_0063` pseudo と `run_0064` split dry は完全一致
  (‖ΔM‖∞ 0.062 %, ε_M 0.0117 %, ε_θ 0.0194°; 軸 M 差 max 1.7e-5)。`run_0065` split + Kw+HK: onset x≈6.7 r_t (T≈250 K)、
  軸 g_max 0.0126 (H₂O の 25 %)、出口軸 T 245→278 K、P +4.5 %、出口軸 M 4.20→3.94 (ε_M 7 %)。品質 PASS・NaN 0、
  残差 warm 床 (NOT CONVERGED 判定、0065 は 2.7 桁降下)。

## 7. 影響範囲

`design/forge_design/{gas/semiperfect.py, evaluate/ic.py, evaluate/runner_axismach.py}`, テスト、
case/16・case/42 の投入スクリプト/yaml/README。既定 (`tp_species: pseudo`) は無変更。

## 8. 完了条件

- [x] split DB の等価性テスト
- [x] Wyslouzil TP split (Kw+HK) が既存 CPG run_0050 を再現 (cell)、node は非粘性で cell と一致
- [x] イソブタン M4.2 H₂O 5 %: pseudo vs split dry 一致、split+凝縮 完走・報告
- [x] status done → accepted、`plans/README.md` 同期

## 変更ログ

- 2026-08-17: 起票・実装 (`mixture_pseudo_species_split`, `paste_isentropic_ic(species_Y=)`, `cpg_field_to_tp`,
  `zero_wall_velocity_ic`, runner `tp_species`/`condensation`, `run_staged(mesh_h5=)`, `interp_field` の凝縮モーメント移植)・
  検証 (Wyslouzil cell/node、イソブタン H₂O 5 %)。**結論: split_h2o は dry で pseudo と同一、凝縮は既存モデルで動く。
  node は非粘性で cell と一致、平面 no-slip 壁クラスタメッシュの node 粘性は申し送り。** status done → accepted。

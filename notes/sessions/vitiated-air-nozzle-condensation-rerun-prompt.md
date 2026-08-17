# 引き継ぎプロンプト: 加熱空気風洞ノズル (case/44) を新しい入口条件で「同じ計算」をやり直す

- **起票**: 2026-08-18。対象 case: `case/44.vitiated_air_wt/` (README が正本、run 一覧あり)。
- **前提知識は不要**: このプロンプトと case/44 README の「新条件」節だけで再現できる。理論は [`methods/design/overview.md`](../../methods/design/overview.md)
  (axis-Mach チェーン)、[`methods/condensation.md`](../../methods/condensation.md) (非平衡/平衡凝縮)。

## 依頼内容 (そのまま貼って使う)

> `case/44.vitiated_air_wt` の「新条件 (va2 → 正本は va2c `run_0079`–`0090` + dry `run_0049`…)」と**同じ手順で**、入口条件だけ次に変えて計算してください:
> 全圧 Pt = ____ Pa、全温 Tt = ____ K、組成 (モル分率) H2O ____ / N2 ____ / O2 ____ / Ar ____ / CO2 ____ (合計 1)。
> 出口径 φ1600 mm、入口配管 φ1600 mm、R=2 / L_U=6 は変えない。M_d = {4.19, 4.3, 4.4, 4.5, 4.75, 5.0} の 6 形状 ×
> {dry, 非平衡凝縮 (Kw+HK), 平衡凝縮} = 18 run。各形状の壁点列 (Euler 版: CSV + **mm 単位 "x r" の .dat** + gmsh .geo を `geometry/` に)、
> 全 run の軸上物理量 CSV (h5 の VALUE 全部 + Tsat 後処理列)、dry 1 点の NASA-CEA 凍結流照合 (`cea/`) を出し、
> README の run 一覧と結果表・図を更新して commit/push まで。

## 手順 (エージェント向け)

1. **組成の質量分率化とガス参照値**: モル分率 → 質量分率 (MW: H2O 18.015, N2 28.013, O2 31.999, Ar 39.948, CO2 44.010)。
   `gas.species` に質量分率、`gas.gamma` はスロート γ\* (`GasSemiPerfect(Y, Tt).gamma_throat()`)、`gas.cp` は cp(Tt)。
   これらは Hall 遷音速の参照値で、MOC/CFD 本体は NASA-9 (semi-perfect) が使われる。CPG は使わない (出口径を 9–12 % 誤る)。
2. **M_d ごとの r_t と背圧**: `r_t = 0.8 / sqrt(A/A*(M_d))` (semi-perfect の面積比、`design/forge_design/gas/GasSemiPerfect.area_ratio`)、
   `spec.p_ambient` は dry 等エントロピー出口静圧の **0.3 倍** (超音速流出保証; 実際は読まれない)。参考実装: 前回の YAML 生成コード
   (`case/44` の履歴、`problem_va2lp_M*_{dry,noneq,eq}.yaml` の冒頭コメント)。L_c は M_d 4.19/4.3/4.4→8、4.5→9、4.75→10、5.0→11
   (単一 quintic の設計可能域; 壁 QA REJECT なら +1)。
3. **YAML**: `problem_va2lp_M<Md>_<kind>.yaml` をコピーし `spec.Pt/Tt/p_ambient/r_throat/M_design`、`gas.species/gamma/cp`、`dv.L_c.value` を書き換える。
   `evaluate` は共通: `nStepOuter 12000, cfl_main 2.0, mid_stage true, tp_species: split_h2o`、noneq は
   `condensation: {condensation: 1, nCondSpecies: 1, condModel: 1, condKantrowitz: 1, condGrowthModel: 0}`、eq は
   `condensation: {condensation: 1, nCondSpecies: 1, condModel: 1, condEquilibrium: 1}`。dry は condensation ブロック無し。
   **H2O が組成に無い場合は split_h2o/凝縮は不可** (dry のみ)。
4. **投入**: `case/44.vitiated_air_wt/run_va2c_batch.py` (凝縮 12) と `run_va2lp_batch.py` (dry を含む 18) を写して run 番号を README の一覧の続きから振り、
   `cd design && ./.venv-opt/bin/python <driver>`。1 run ≈ 1.5 分 (RTX 3060)。既存 run ディレクトリは使い回さない。
5. **後処理**: `python3 axis_csv_va.py run_00NN_*/` (軸 CSV: h5 全 VALUE + `pH2O_post/Tsat_post/subcool_post/S_post`)、
   `design/.venv-opt/bin/python exit_table_va.py "run_00NN_*" --out study_<tag>_exit.json` (物理出口 2 r_t 上流の軸値・質量平均)、
   `summarize_va2.py <tag>` (図: 軸 T/Tsat, S, g) と表、`export_wall_va.py --problem problem_<tag>_M<Md>_dry.yaml --kind euler --tag <tag>_M<Md>` (点列 CSV + .geo)、
   点列 → mm .dat は README「成果物の所在」の `geometry/` 生成スニペット (x[mm] r[mm] 2 列、ヘッダ `#`)。
   CEA 照合: `cea/va2_cea.inp` を写して p,bar / t,k / moles / supar (=A/A\*(M_d)) を書き換え `echo <name> | ../../../.venv-cea/nasa_cea/FCEA2`、
   `cea_check_va2.py` の PT/TT/S/CEA_OUT/CEA_CASE/SUPAR を合わせて実行 (dry の M4.19 相当 1 点で十分; T/ρ/u/a/M ≤0.05 % 一致が期待値)。
6. **必須チェック** (AGENTS.md): 全 run `MESH_QUALITY.txt` PASS、`res_12000.h5` の NaN 0、`check_convergence` VERDICT を貼る (warm 床 plateau で
   NOT CONVERGED が通常; 「収束」とは書かない)、軸 M が 8k→12k で凍結。
7. **README**: 「新条件」節に倣って条件表・結果表・図・run 一覧行を追記。commit は英語命令形、`res_*.h5` は含めない (config/CSV/json/png は `git add -f`)。

## 既知の罠・注意

- **forge binary は 2026-08-18 の修正 (SLAU 二相面温度 [エネルギー保存]、平衡凝縮 `condEquilibrium`、`condTsat`、CEA 由来の潜熱 L=h_v−h_l) を含むものを使う** (commit 5b55a2a7 以降) (`solver_density_cuda/build/forge`、
  ソースが変わっていたら `ninja -C solver_density_cuda/build`; solverConfig.hpp を触ったら clean rebuild)。
  それ以前の binary は凝縮帯で全エンタルピーが +0.6 % 非保存 (T +4–5 K)。**凝縮 run は軸の全エンタルピー h0=(roe+P)/ro が一定 (±0.1 %) であることを報告に添える**。
- 出口直前 1–2 r_t で平衡 run の g が下がる/S<1 になるのは平衡流れ固有の集束圧縮波 (BC 無関係、背圧 0.3 倍で不変を確認済み)。物理出口の議論は 2 r_t 上流断面。
- H2O 以外の凝縮種は未対応 (CO2 等は不可)。組成に H2O が無ければ dry のみ。
- 凝縮の判定は軸上の `condS_0` (過飽和度) と `g_0/Y_H2O`。非平衡の Wilson 点は概ね T≈205–210 K (S≈100–300)、平衡は S=1 到達点 (T≈260 K)。
- 「出口軸 M」は `metrics.json` の M(x_E) (設計終点) で、凝縮帯の途中になり得る。**出口の議論は `exit_table_va.py` の物理出口断面 (質量平均) を使う**。
- Tsat は凝縮 ON の run にだけ h5 (`condTsat_0`) に出る。dry は `axis_csv_va.py` の `Tsat_post` (同じ Murphy–Koop、0.01 K 一致)。
- 出口 T が 200 K を切る M_d では NASA-9 の下限 (200 K) 外挿になる (M5.0 dry で 190 K)。結果は連続だが注記する。
- 別セッションが同じ作業ツリーで forge を再ビルドすることがある。run 失敗 (rc 127) はそれ。`run_0019` のような回帰確認 (dry 1 本) を挟むと安全。

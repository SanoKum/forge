# codex レビュー (2): 2026-09-05 午後の対応レビュー

codex-cli 0.153.4 (read-only)。対象: 前回レビュー対応 (機構検証・check_quasisteady 反応流量・Cabra 出口逆流・文書同期)。

## (1) 総評

午後の作業は、問題を「残差」ではなく着火位置・出口流量の時系列で捉え直した点で有益です。特に、Cabra の出口変動と環状逆流を発見したことは重要です。

ただし、現時点で次の2点は P0 級です。

- 0-D 機構比較の coflow 組成に誤りがあり、掲載した着火遅れ数値は再計算が必要です。
- 「付着の支配因子は TCI 欠如と確定」は、非収束計算と Li 1機構だけの A/B からは言い切れません。「TCI/着火統計クロージャ不足という仮説を強く支持」までが妥当です。

`check_quasisteady.py` の Cabra/node/軸対称での積分値自体は概ね正しい一方、NaN の黙殺や cell モードの出口抽出には誤判定余地があります。

コード・文書は変更していません。

## (2) 観点別レビュー

### (a) `check_quasisteady.py` の抽出子

良い点：

- Cabra のような node・直交出口では、`x >= x_max - tol` により出口境界177節点が正しく選ばれています。
- 台形重みは、端点が半セル、内部点が左右幅の半和になっており正しいです。
- 軸対称の `2πr` 重みも、`r≥0` の Cabra では正しいです。
- `h0=(roe+P)/ro` は対流全エンタルピーとして正しいです。
- 実データで主場の節点積分と `res_outlet` の境界値を照合すると、例えば `run_0072` step 30000 は 0.04487 対 0.04481 kg/s で一致しました。したがって今回検出した大きな質量流束変動は、積分コードだけが作った偽信号ではありません。
- `exit_y_存在しない種` は `chem_bad` で error になり、unknown-quantity ガードとの整合も取れています。[check_quasisteady.py](/home/sano/work/forge-chem/solver_density_cuda/tools/check_quasisteady.py:396)

問題点：

1. **NaN/例外を黙って除外するため、偽の `STEADY` があり得ます。**

   抽出中の全例外を `nan` に変換し、`classify()` が非有限値を除外しています。[check_quasisteady.py](/home/sano/work/forge-chem/solver_density_cuda/tools/check_quasisteady.py:417)

   特に `ignx` は、火炎が途中で消えて末尾がすべて NaN でも、過去の有限区間だけから `STEADY` を返し得ます。`tmax`、出口量も、破損・欠落・ゼロ除算・NaN を含む最新 snapshot を無視して過去だけで判定できます。少なくとも末尾 snapshot の非有限値、抽出例外、有限率不足は error または `TRANSIENT-UNSETTLED` にすべきです。

2. **`exit_y_<種>` は「主要種流束」ではありません。**

   実装は符号付き正味質量流束で重み付けした平均質量分率です。[check_quasisteady.py](/home/sano/work/forge-chem/solver_density_cuda/tools/check_quasisteady.py:367)

   逆流が強いと分母がゼロに近づき、範囲外の値や巨大値になり得ます。また、P0-2 が要求した「主要種流束」は本来 $\int\rho uY_s\,dA$ です。以下を分ける方が安全です。

   - 符号付き正味種質量流量
   - 正方向流出だけの組成
   - 逆流量
   - 正味質量流量

3. **cell/unstructured モードの出口選択は一般には頑健でありません。**

   `x_max` 近傍のセル重心は必ずしも出口面全体を構成しません。`tol` を広げると複数の x 列や同一 y の重複点を一つの台形列として積分してしまいます。また node/cell 判定を座標数の一致だけに頼るため、偶然同数のメッシュでは誤選択します。[check_quasisteady.py](/home/sano/work/forge-chem/solver_density_cuda/tools/check_quasisteady.py:321)

   `mesh.discretization` と出口 `physID`、可能なら `res_outlet_*` を使うべきです。

4. `exit_hflux` は対流エンタルピー流束だけで、熱伝導・種拡散エンタルピー・粘性仕事を含む保存則監査用の「全エネルギー流束」ではありません。コメントには限定が書かれていますが、plan 側でも区別が必要です。

5. `read_solver_config()` が YAML 読込失敗を `{}` に落とすため、軸対称指定を読み損ねても平面積分を続けます。ここも error にすべきです。

### (b) 機構検証と TCI 結論

最大の問題は coflow 組成です。

[ign_delay_mech_compare.py](/home/sano/work/forge-chem/case/48.cabra_h2n2/ign_delay_mech_compare.py:29) は、

```python
O2:0.1474, H2O:0.0989, N2:0.8537
```

で、合計が 1.1 です。Cantera が正規化するため、実際は概ね `O2=0.1340, H2O=0.0899, N2=0.7761` となり、CFD の `O2=0.15, H2O=0.099, N2=0.751` と違います。したがって README の 1.32/1.54/2.73/14 ms と比率は、正しい条件で再計算するまで確定値として扱えません。

そのほかの交絡：

- Cantera比較は各機構付属の NASA 熱力学を使用します。一方 `run_0072` は Jachimowski と同じ CEA `species_db.yaml` を使い、反応表だけ Li に交換しています。前進反応は Li でも、逆反応と反応熱は「Li付属熱力学の Cantera」と同一ではありません。
- Forge の Li+Troe 経路について、同じ CEA熱力学を使った 0-D Cantera/Forge 一致がまだありません。`run_0072` の完走は「読めて動いた」証拠ですが、Troe を含む反応率の定量一致証明ではありません。
- 0-D は `max dT/dt`、CFD は `Y_OH>2e-4`、実験は OH 発光です。機構差を火炎基部差へ結びつけるには閾値感度と定義の対応が要ります。
- Li と Jachimowski は両方とも、今回の平均場では付着解へ到達できるほど速い可能性があります。これだけでは、TCI、coflow温度、入口乱流、管内発達、リップ熱条件、格子拡散の寄与を順位付けできません。
- 両 run は定常擬似時間で NOT CONVERGED です。step 20k/26k の比較は物理的な伝播時間比較ではありません。

したがって妥当な結論は、

> Li 2004への交換だけでは付着解を解消できず、平均場反応率／着火統計クロージャ不足という仮説を支持する。

です。「支配因子と確定」は過大です。

### (c) Cabra 出口環状逆流

環状逆流の存在自体は信頼できます。軸対称では外側環状部が `2πr` により強く効くため、y=54–78 mm の逆流が総質量流量を大きく変えるのも整合的です。

ただし、「反応が引き金」は相関までは確認できていますが、原因特定ではありません。反応 ON はいずれも非収束で、定常擬似時間の振動は物理的な非定常現象ではありません。

レビュー時点の既存出力では、`case/48.cabra_h2n2/run_0073_react_li_cont/` は step 20000 まであり、`check_quasisteady` は次のままです。

- `exit_massflux`: **DRIFTING**, drift 41%
- `exit_hflux`: **DRIFTING**, drift 36%
- `tmax`, `ignx`: STEADY

したがって、まだ頭打ちではありません。

次に疑う順序は以下がよいです。

1. **外周 slip による閉じ込め**

   実験の加熱 coflow は開放噴流なので、r=100 mm の slip は半径方向の entrainment を禁止します。外周を開放圧力境界にする A/B は妥当です。ただし現行 `outlet_statPress` は逆流時の組成・エントロピーを内部から外挿するため、外部の既知 coflow/ambient 組成を与える物理的 farfield BC にはなっていません。[boundaryCond_d.cu](/home/sano/work/forge-chem/solver_density_cuda/cuda_forge/boundaryCond_d.cu:599)

2. **出口距離・領域半径感度**

   出口で大規模逆流が横切るなら、境界が内部解へ干渉しています。x方向延長と半径拡張を確認すべきです。

3. **低Mach前処理／出口特性**

   この経路では同日に3件のバグが見つかっています。`lowMachPrecond: 0` の低CFL比較、precond 2、dual-time の3者比較が必要です。また TP出口は局所 γ を使うよう修正済みですが、`P/ρ^γ` は可変比熱気体の厳密な等エントロピー関係ではありません。

4. **dual-time**

   妥当ですが、正しい外周境界の代替にはなりません。物理時間で複数 flow-through time を回し、各 physical step の inner convergence と時間平均±振幅を確認して初めて物理的振動と呼べます。

5. **格子・軸対称**

   品質 PASS は確認済みですが、せん断層・外周・出口の格子感度は未確認です。中半径の現象なので軸特異点が第一容疑ではありませんが、同一メッシュ/BCの別解法や SU2 比較は有効です。

### (d) 前回 P0/P1 の仕分け

| 指摘 | 判定 | コメント |
|---|---|---|
| P0-1 設計QoI・合格基準 | 開 | ユーザ判断待ち |
| P0-2 SST/準定常性 | 一部対応 | 抽出子と問題発見は前進。ただしツールのNaN耐性、Cabra/BKの非収束、dual-time確認は未完 |
| P0-3 機構着火遅れ | 対応不十分 | coflow組成バグ、実加熱器圧力なし、Mueller/Keromnesなし、Forge Li 0-Dクロスチェックなし。「done」は不可 |
| P0-4 TCI別plan | 開 | plan未作成、ユーザ判断待ち |
| P1-5 BK不確かさ | 開 | 出口再確認のみ。機構A/B、壁温、入口乱流、3D blockage整理は未実施 |
| P1-6 Cabra入口・管・リップ・格子 | 開 | 未実施 |
| P1-7 end-to-end収支 | 開／残作業表から脱落 | 元レビューの質量・元素・全エンタルピー収支等が §5.1 に明示されていない |
| docs/index同期 | 概ね閉 | 古いPhase 0記述は更新されたが、新記述に過大主張あり |
| run表 1 run=1行 | 開 | 「ユーザ判断待ち」ではなく AGENTS.md 上は必須 |
| OH閾値統一 | 一部対応 | 新規ツールは2e-4だが、BK READMEの1e-4値と旧解析値が注記なしで混在 |
| `chemistry.enabled:0` 回帰 | 開 | エビデンス未追加 |
| 0-D自動テスト | 開 | 手動のまま |
| repo内一次成果物 | 一部対応 | スクリプトは追加されたが、引用するCSV/PNG/logはgit管理されておらず、外部Claude artifactも残存 |

### (e) 文書上の誤認・過大主張

- plan §5.1 の「P0-3 done」「TCI欠如と確定」は撤回が必要です。[plan](/home/sano/work/forge-chem/plans/active/chemistry-finite-rate-h2.md:66)
- §5.1 は Cabra変動を「37–48%」としますが、§9 と実ツールでは `run_0072` が 62.9%です。「37–63%」が整合します。
- P1-5 は「P0-3 A/B待ち」のままで、同じ表では A/B done とされており矛盾しています。
- plan のメタ情報は依然「Phase 3 BK一次結果 2026-09-04」で止まり、Cabraが §6 の検証ケース・合格基準へ追加されていません。
- README の `run_0072` は「付着で安定」と書きつつ、同じ行で NOT CONVERGED、出口流束とTmaxは DRIFTINGです。「`ignx` が閾値上0で維持」程度に限定すべきです。[README](/home/sano/work/forge-chem/case/48.cabra_h2n2/README.md:51)
- README の機構比較節末尾は「forge用変換は未」「次: Li A/B」のままで、直上の `run_0072` と矛盾します。[README](/home/sano/work/forge-chem/case/48.cabra_h2n2/README.md:65)
- methods は BK/Cabra双方で TCI欠如が支配としていますが、Li A/B は Cabraのみです。[chemistry.md](/home/sano/work/forge-chem/methods/chemistry.md:9)
- methods の「Li 2004, Burke 2012をそのまま読めることを `run_0072` で確認」は誤りです。`run_0072` が実行したのは Li のみです。[chemistry.md](/home/sano/work/forge-chem/methods/chemistry.md:147)
- 「Q1DはCanteraと一致」は範囲が広すぎます。OH/NOの部分的一致であり、`check_convergence` は NOT CONVERGED、T/Mの絶対比較にも制約があります。
- 「H/d再現には PDF/CMC級が要る」も、文献調査自身の「最も強い実績であり、のみとは断定できない」に合わせて弱めるべきです。

## (3) 優先度付き推奨アクション

### P0

1. `X_COF` を正しい balance に直して全機構比較を再実行し、README・planの数値を更新する。
2. Liについて、ForgeとCanteraを同じ CEA熱力学・同じ反応表で0-D照合し、Troe、逆反応、着火遅れを確認する。
3. P0-3を reopen し、少なくとも正しい Cabra条件、実加熱器圧力、閾値感度を終える。「TCI確定」は「仮説を支持」へ修正する。
4. `check_quasisteady.py` は抽出例外・末尾NaN・ゼロ流量を error にし、`exit_species_massflux` と正方向流出組成を分離する。
5. `run_0073` 完了後も DRIFTING なら、外周 open/pressure A/Bを優先し、次に領域延長、precond 0、dual-timeの順で切り分ける。

### P1

6. node/cell、平面/軸対称、逆流、NaN、複数出口を含む小さな合成HDF5テストを抽出子へ追加する。
7. TCI planでは「CMC/PDFを実装すること」より先に、`H(T_c)`、OH閾値、平均/RMS、質量・元素・エンタルピー収支を合格基準として固定する。
8. plan・methods・READMEの過大主張、37–48/63%不整合、P1-5状態、Li/Burke記述、旧1e-4閾値を同期する。
9. run表の1 run=1行化、非反応ビット不変回帰、0-D自動テスト登録を残作業として継続する。

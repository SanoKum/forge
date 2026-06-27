# 分析依頼プロンプト: 多成分 TP の陰解法 (block-DPLUR) で安定 CFL が ~1 に頭打ちする原因

> このファイルは **根本原因の分析を依頼するためのプロンプト**である(実装依頼ではない)。
> forge リポジトリにアクセスできる AI にも、できない AI にも渡せるよう、観測事実・現状アーキ
> テクチャ・切り分け済みの手がかりを自己完結的にまとめてある。**最初にこの全文を読み、
> 「何が安定 CFL を律速しているか」を分析・仮説化し、検証実験と根拠を提示してほしい。**
> コードを直す前段の「原因特定」が成果物。リポジトリにアクセスできる場合は該当ファイルを
> 読んで根拠を示すこと。

---

## 0. 1 行サマリ

forge(密度ベース圧縮性 CFD, `solver_density_cuda/`)の定常陰解法 `block-DPLUR`
(`timeIntegration: 11`, `blockDPLUR: 1`)で、**多成分 thermally-perfect gas (TP, `thermalMethod==2`)
+ RANS-SST** のケースだけ安定 `cfl_pseudo` が **~1.0** に頭打ちする。**単成分 TP や多成分 CPG では
もっと高い CFL が回る**ことは確認済み。なぜ「多成分 × TP」の組み合わせだけ CFL が上がらないのか、
律速機構を特定してほしい。

## 1. 観測事実(切り分け済み — ここが最重要の手がかり)

検証ケース: `case/28.cutler_coaxial_jet`(軸対称・超音速同軸ジェット。中心 He-O2 ジェット
Mach1.8 を空気 coflow が囲む。非反応。化学種 `[He, O2, N2]`)。数値は全て 2 次風上 +
Venkatakrishnan, SLAU 対流, block-DPLUR(5×5), `nStepInner=5`, `implicitRelax=0.7`,
`speciesDiffusionMethod=1`, `viscMethod=2`(kinetic theory 混合輸送)。実効定常 CFL は
`cfl_pseudo`(後述)。発達場からの「継続」で `cfl_pseudo` を振った結果:

| 設定 | 安定上限 `cfl_pseudo` | 備考 |
| --- | --- | --- |
| **多成分 TP + SST** | **~1.0** | 1.0 安定 / 1.2 で step247 発散 / 1.5→step116 / 2→step50 / 4→step1 / 8→step9 |
| **多成分 CPG + SST** | **≥2.0** で安定(発散せず) | 同一メッシュ・同一 BC・同一スキーム。**EOS だけ CPG(単一 γ/cp)** に変えたもの |
| **単成分 TP**(別ケースでユーザ確認済) | 多成分 TP より明確に高い | 化学種なし |

**この対比が決定的**: CPG↔TP の唯一の差は **EOS が組成 `Y` に依存するか否か**。
- CPG(`thermalMethod==0`): 化学種 `roY_s` は移流スカラとして輸送されるが、**EOS(`T,P` の決定)は組成に依存しない**(単一 γ/cp)。→ cfl=2 でも安定。
- TP(`thermalMethod==2`): `T=T(e,Y)`(NASA-9 で `e` から `T` を反転)・`P=ρ R_mix(Y) T`・`h_s(T)` と、**EOS が組成に強く依存**。→ cfl~1 で頭打ち。

→ 化学種「輸送」自体や SST 単体ではなく、**`roe`(エネルギー)と `roY`(組成)の EOS を介した結合**が
律速の第一容疑。発散は CFL が高いほど早く崩壊する典型的な数値不安定(物理発散でなく LHS 演算子の問題)。

**He 特有の事情**: He は単原子で `cp≈5193 J/kgK`(空気の ~5 倍)。生成エンタルピーは He/O2/N2 とも
≈0(本ケースに H2O のような生成熱の大きい種は無い)。つまり **「組成差による顕熱エンタルピー差」が
大きい**(`h_He(300K)≈1.56 MJ/kg` vs `h_N2≈0.31 MJ/kg`)。後述の既存知見は生成エンタルピーを犯人と
したが、本ケースは生成熱ゼロなので **顕熱 cp コントラストが増幅器**になっている可能性が高い(要検証)。

## 2. 現状アーキテクチャ(正確なコードポインタ — `implicitNonlinearUpdate`)

`solver_density_cuda/main.cpp` の `implicitNonlinearUpdate`(≈L900)が 1 擬似時間ステップで行う更新は
**完全に segregated(3 つの独立陰解法)**で、相互の off-diagonal 結合が無い:

1. **流れ 5 変数** `[ro, roUx, roUy, roUz, roe]`: `blockDPLURSolve` → `applyBlockImplicitCorrection`。
   5×5 block-DPLUR(`cuda_forge/timeIntegration_d.cu` の `namespace block_dplur`)。対流は**真の分離 A±**
   固有値ヤコビアン、粘性はスペクトル半径を対角に加算。**EOS ヤコビアン**(`∂P/∂(保存量)`, 可変 γ)は
   一般 EOS 閉形式 `eos_split_jacobian_general_closed` / `build_jacobian_split`(precond=2)で評価。
   **ただしこの block 内では組成 `Y` は固定(frozen)** — `roY` はこの解法では動かさない。
2. **SST** `[k, ω]`: `applySSTPointImplicit`。流れとは分離(segregated)、**生産項 Pk/Pω は陽的(lagged)**、
   消散項のみ対角に陰的(`ransSource_d.cu` の `src_jac_*`)。壁近傍 ω~10^? で stiff。
3. **化学種** `roY_s`: `speciesImplicitCoupled` が真(config `time.deltaT.speciesImplicitCoupling: 1`)なら
   `speciesImplicitDPLURSolve_d_wrapper`(`cuda_forge/speciesTransport_d.cu` L491)。
   **緩和整合 scalar-DPLUR**: `dq=0` から `nStepInner` 回 Jacobi sweep、流れ block と**同一 `dt_local`・
   同一 `implicitRelax`・同一 sweep 回数**。ただし**スカラー独立**で `roe↔roY` のクロス項は無い。
   その後 `speciesRenormalize_d`(`ρY_s≥0`, `ΣρY_s=ρ` を強制再正規化)→ `speciesPrimitive`。

残差(RHS)側は整合している: SLAU 対流流束は**被移流エンタルピーをセル組成で再構成**しており、
化学種エンタルピー輸送は対流流束と二重計上しない(`docs`/README 参照)。したがって
**不安定は RHS(残差)でなく LHS(陰解演算子)の分離**にあると見られる。

## 3. 主仮説(検証して棄却/採択してほしい)

**H1(本命): `roe`↔`roY` の擬似時間緩和ミスマッチ + EOS クロス項欠如。**
`roe` は 5×5 block で組成固定のまま緩和され、`roY` は別 scalar-DPLUR で緩和される。両者の更新が
同一線形系に入っていない(off-diagonal `∂(エネルギー残差)/∂roY`, `∂(種残差)/∂roe` が無い)。
非線形 commit 後の `T=T(e,Y)` 反転で、`e` と `Y` の更新量が不整合だと **温度ジャンプ**が出て増幅する。
増幅率は組成差による顕熱エンタルピー差(He の高 cp)に比例 → 多成分 TP 特有。CPG は EOS が `Y` 非依存
なのでこの経路が消え、安定上限が上がる(観測と整合)。
- **検証案**: (a) `speciesImplicitCoupling` を 0/1 で安定上限が変わるか。(b) `nStepInner` を増やしても
  上限が上がらない(=外側分離の問題で内側 sweep 不足でない)ことの確認。(c) He を O2/N2 に置換した
  「cp コントラストの小さい多成分 TP」で上限が回復するか(顕熱コントラストが効くなら回復するはず)。
  (d) 既存の `physProp.thermoHrefTemp`(生成エンタルピー除去オフセット)は**本ケースに効かないはず**
  (生成熱ゼロのため)— これが効かないことを示せれば「犯人は生成熱でなく顕熱コントラスト」を裏付ける。

**H2: 一般 EOS block ヤコビアンの不正確さ(可変 γ/混合)。**
`eos_split_jacobian_general_closed` が混合・組成勾配の大きい contact で `∂P/∂(ρ,ρe)` や音速・固有値を
取りこぼし、5×5 の `|A|` が真のスペクトルを下回って局所 dt 過大になる可能性。
- **検証案**: 単成分 TP と多成分 TP で、組成一様領域では block ヤコビアンが一致するか(構成上)。
  最初に発散するセル(下表 detectNaN)が **contact/せん断層の組成勾配が最大の場所**かを確認。

**H3: SST segregation の上乗せ。**
SST 自体が `cfl_pseudo` を律速している可能性。生産項 lagged + k-ω stiff。
- **検証案**: 同一 TP 多成分で **SST off(層流)** にして安定上限を測る。SST on/off で上限が変われば SST も寄与。
  (本ケースは自由ジェットで壁 no-slip 無し → 壁 ω stiff は軽いはず。だが μ_t 経由の結合はある。)

**H4: 実現可能性再正規化(`speciesRenormalize`)の非線形キック。**
`ΣρY_s=ρ`/`ρY_s≥0` の強制クリップが、陰解で得た `δroY` を毎ステップ非線形に書き換え、
`roe` 側と不整合を作る可能性。
- **検証案**: クリップ発動頻度・発動セルと最初の発散セルの相関を見る。

**H5: 軸対称ソース項の陰解扱い。** 近軸 `1/r` ソースが TP × 組成で増幅する可能性(優先度低だが念のため)。

## 4. 依頼する具体作業(分析の出力)

1. 上記 H1–H5 を **証拠付きで採択/棄却**し、律速機構を 1–2 個に絞る。コードを読める場合は
   `main.cpp` / `timeIntegration_d.cu`(`block_dplur`)/ `speciesTransport_d.cu`
   (`speciesImplicitDPLURSolve`)/ `dependentVariables_d.cu`(`T=T(e,Y)` 反転)の該当箇所を引用して根拠を示す。
2. **最初に発散するセル**を `time.deltaT.detectNaN: 1` で特定し(`res_nan_<step>.h5` をダンプ)、
   その場所の物理(組成勾配・温度・音速・μ_t)を述べる。CFL を 1.2/2/4 と変えて発散セル・発散ステップが
   どう動くかを表にする。
3. **線形安定性の観点**: segregated な 3 分割反復(flow block / SST / species)を 1 つの結合反復と見たとき、
   どの off-diagonal を落としていることが増幅につながるかを、簡単なモデル(`roe`-`roY` 2×2 緩和系など)で説明する。
4. **最小コストの決定実験**を提案(理想は既存バイナリ + config 変更だけで切り分くもの)。
   例: ①`speciesImplicitCoupling` 0/1、②SST on/off、③He↔追加 N2 で cp コントラスト除去、
   ④`thermoHrefTemp` on/off、⑤`nStepInner` sweep、⑥`implicitRelax` を流れ/種で変える。
   それぞれが H1–H5 のどれを支持/棄却するかを事前に明示する。
5. 律速が H1 なら、恒久修正の方向性(化学種を 5×5→`5+nSpecies` block に同梱する full coupling=既存
   `species-in-dplur-session-prompt.md` の案A、または `roe↔roY` の 2-way 整合補正)を**二重計上の罠**
   (下記)に触れつつ評価する。**実装はしない、方向性の妥当性だけ述べる。**

## 5. 既存知見(前セッション。本ケースとの相違に注意)

`.github/plans/species-in-dplur-session-prompt.md` に、多成分 TP 不安定の前分析がある。要旨:
- 主因 = `roe`(block)と `roY`(別 point-implicit)の緩和ミスマッチ → Newton 温度ジャンプ。**ただし
  そこでは増幅器を「生成込み絶対エンタルピー」とし、`h_H2O(298K)≈−13.4 MJ/kg` の H2O を犯人とした**
  (N2/O2 は `h≈0` で無害と結論)。対策 `physProp.thermoHrefTemp`(sensible-enthalpy datum offset, commit
  `280d11d`)で安定上限 1→2 に改善したが**2 倍止まり・残差プラトー**、本質未解消。
- **本ケースとの重要な相違**: Cutler は **He/O2/N2 で生成熱ほぼゼロ**。それでも上限 ~1。よって
  「生成熱が犯人」では説明できず、**顕熱 cp コントラスト(He 5×)が増幅器**という仮説(H1-c/d)を
  検証する価値が高い。前分析の `thermoHrefTemp` が本ケースで効かないことは、その傍証になる。
- **死蔵カーネルの罠**: `speciesTransport_d.cu` の `species_energy_correction_kernel`
  (`roe += Σ(roY−roYN)h_s`)は未配線(意図的)。配線すると対流流束が既に整合輸送している化学種
  エンタルピーと**二重計上**して O2 まで発散する(過去 run で確認済)。修正方向を論じる際はこの轍を踏まないこと。

## 6. 環境・規律(リポジトリにアクセスする場合)

- **実効 CFL**: 定常(`unsteady:0`)は `cfl` でなく **`cfl_pseudo`** が局所刻みを決める
  (`cuda_forge/setDT_d.cu` `setDT_d_wrapper`)。`cfl` は表示用。`procedures/solver-settings.md` 参照。
- **build**: native `solver_density_cuda/.build-native/relwithdebinfo/forge`(nvcc 12, sm_86)。
  `tools/build_native_wsl.sh`。**構造体(`solverConfig.hpp` 等)変更後は full rebuild**(stale build trap)。
- **収束判定の規律**: `rms_ro` だけで判断しない。`tools/check_convergence.py <run>` の VERDICT を必ず使い、
  全残差列(`rms_ro..roe`, `rms_roK/roOmega`)・NaN・`ΣY=1`・`Y∈[0,1]` を確認する。
- **run 運用**: 既存 `run_*` を使い回さず複製した新 `run_NNNN_<slug>` で実行。`case/28.cutler_coaxial_jet/README.md`
  の「## 計算 run 一覧」を同期。本分析の元データは `run_0037`(CPG-SST)/`run_0038`(TP-SST cfl0.5)/
  `run_0039`(TP-SST cfl1)。

## 7. 期待する出力フォーマット

1. **結論**(律速機構 1–2 個、確信度)。
2. **根拠**(コード引用 + 観測事実 + 発散セルの物理)。
3. **H1–H5 の採否表**(各仮説: 支持/棄却/未判定 + 理由)。
4. **決定実験リスト**(config 差分・期待結果・どの仮説を切り分けるか)。
5. **修正方向の評価**(実装はしない。full block 同梱 vs 2-way 整合補正 vs SST 結合、コスト/リスク)。

> 注意: 推測で断定しない。リポジトリにアクセスできるなら必ずコードを読んで裏を取り、
> できる実験は `check_convergence.py` の VERDICT を貼って示すこと。

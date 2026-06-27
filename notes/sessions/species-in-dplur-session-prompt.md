# 引継ぎプロンプト: 化学種 `roY_s` を block-DPLUR に結合する (roe↔roY coupling)

> 次セッションの起点。冒頭にこのファイルを参照させること。
> **方針確定 (2026-06-20)**: 安定 CFL 律速は `roe`(block)と `roY`(別解法)の結合欠如。
> 緩和整合 (案B) は実装済みだが不十分。**恒久対応として線形系を結合する (II) で進める**
> — まず **案C(EOS クロス項だけ足す partial coupling、診断兼用)**、不足なら **案A(full block)**。

## ゴール

多成分 thermally-perfect (`thermalMethod==2`) の陰解法で、エネルギー `roe` と組成 `roY_s` を
**同一線形系で結合**し、EOS (`T=T(e,Y)`, `P=ρR_mix(Y)T`, `h_s`) を介した `roe↔roY` の擬似時間
緩和ミスマッチを**構造的に解消**する。目標は安定 `cfl_pseudo` 上限を現状(冷流 He/air で ~1、
hot N2/H2O で ~2)から引き上げ、残差プラトーを破って真の定常収束に到達すること。

## 確定した律速機構 (2 セッションぶんの結論)

`roe` を 5×5 block-DPLUR で、`roY` を**別の線形系**で解いているため、両者の更新が同一連立に
入っておらず off-diagonal `∂R_energy/∂(ρY_s)`(組成→圧力・被移流エンタルピー)が欠落している。
非線形 commit 後の `T=T(e,Y)` 反転で `e` と `Y` の更新量が不整合だと**温度ジャンプ**が出て増幅する。

**増幅器は「組成を変えたとき `T,P,h` がどれだけ動くか」**。これは 2 形態ある:
- **生成エンタルピー**(反応系/H2O 等): `h_H2O(298K)≈−13.4 MJ/kg`。hot N2+H2O ケースの主犯。
- **顕熱 cp コントラスト**(非反応でも効く): He は `cp≈5193`(空気の ~5 倍)。**He/air は生成熱ほぼ
  ゼロなのに上限 ~1** → 顕熱コントラストが増幅器。`case/28.cutler_coaxial_jet` で確認(下記)。

**決定的な切り分け (case/28, 2026-06-20)**: 同一メッシュ・BC・スキームで EOS だけ変えると、
- **多成分 CPG**(EOS が `Y` 非依存)→ `cfl_pseudo=2` で安定。
- **多成分 TP**(EOS が `Y` 依存)→ `cfl_pseudo~1` で頭打ち(1.2→step247, 2→step50, 4→step1 発散)。
- 単成分 TP はさらに高 CFL 可。

→ 律速は化学種「輸送」でも SST 単体でもなく、**`roe↔roY` の EOS 結合**で確定的に説明できる。
分析の全文・H1–H5 の切り分けは [`cutler-tp-multispecies-cfl-analysis-prompt.md`](cutler-tp-multispecies-cfl-analysis-prompt.md)。

## 既に打った手 (ここから先へ進む)

- **案B = 緩和整合 scalar-DPLUR(実装済・不十分)**: `time.deltaT.speciesImplicitCoupling: 1` で
  `speciesImplicitDPLURSolve_d_wrapper`(`cuda_forge/speciesTransport_d.cu` L491)。`roY` を `dq=0` から
  流れ block と**同一 `dt_local`・`implicitRelax`・`nStepInner` sweep** で緩和。**緩和速度は揃うが、
  線形系は別のまま(off-diagonal 無し)** → case/28 TP-SST で上限は依然 ~1。緩和整合だけでは不足、を実証。
- **sensible-enthalpy datum offset(実装済・部分的)**: `physProp.thermoHrefTemp`(commit `280d11d`)。
  生成エンタルピーを除き hot H2O 系の上限を 1→2 に改善。**だが case/28(生成熱ゼロ)には効かないはず** —
  顕熱コントラストが犯人なので。これが効かないことは「犯人は生成熱でなく結合欠如」の傍証になる(検証項目)。
- **死蔵カーネルの罠**: `speciesTransport_d.cu` の `species_energy_correction_kernel`
  (`roe += Σ(roY−roYN)h_s`)は未配線(意図的)。配線すると対流流束が既に整合輸送している化学種
  エンタルピーと**二重計上**して O2 まで発散(過去 run で確認)。**結合は LHS(ヤコビアン)で行い、
  RHS にエンタルピー補正を足し込む方向には行かない。**

## 実装・検証結果 (2026-06-20): 案C は CFL 上限を上げない (decisive negative)

案C (block-triangular roe↔roY coupling) を実装・検証した。結論: **案C は実装は正しいが、継続 CFL
上限を上げない**。これは「deferred (RHS 移項) 結合は陰解ブロックの線形安定限界を動かせない」ことを示す。

- **実装** (`speciesImplicitCoupling: 2`): `cuda_forge/species_eos_coupling_d.cuh`(解析 JVP)、
  `cuda_forge/speciesTransport_d.cu`(射影 `species_eos_project_tangent_d` / δp_Y `species_eos_dp_cell_d`
  / クロス流束 `species_eos_cross_flux_d` / 最終 commit `species_eos_final_commit_d` と各 wrapper)、
  `main.cpp` `implicitNonlinearUpdate` の mode-2 並び替え(species 予測→射影→δp_Y→res_roe 移項→flow block
  →最終 commit `ρY=ρY^N+z+Y^N δρ`)。単成分/CPG/mode 0,1 は predicate でゲートし不変。
- **解析 JVP 検証 (host FD)**: `tests/unit/test_species_eos_cross.cpp`。δp_Y/δT_Y を中央差分と比較し
  He/O2/N2(300/800K)・hot N2/H2O(1500/500K)・N2/O2 で **rel ~1e-11 で一致**。式・datum 整合は正しい。
- **injection が活きていることの確認**: mode-2 (run_0048) vs mode-1 (run_0039) を同一 restart・cfl=1 で
  比較 → 残差軌跡は ~5-6 桁目で相違(no-op でない)。ただし効果は **rms_roe で ~0.1% と微小**。
- **CFL 上限 (case/28 TP-SST, run_0038 発達場から継続)**: mode-2 の発散ステップは **mode-1 とほぼ同一**
  (cfl1.2→step254 vs 247 / cfl2→50 vs 50 / cfl4→1 vs 1)。**上限 ~1.0 のまま、改善せず**。

**解釈**: 案C は block-triangular = `A_QY δY` を予測 δY で RHS へ移項する deferred (Gauss-Seidel) 結合。
これは定常解の**整合性 (残差フロア)** は改善しうるが、高 CFL の発散を決める**flow block 自身 (`A_QQ`) の
線形安定性は変えられない**(LHS が不変だから)。よって:
- **H1 (segregated roe↔roY relaxation mismatch) は本ケースの CFL 律速ではない**、または、結合を
  **LHS に入れない限り (案A) CFL は上がらない**、のいずれか。案C だけでは両者を完全には分離できないが、
  「deferred では不十分」は確定。
- CPG(cfl2 安定)↔TP(cfl~1)の差は segregated 結合では説明できず → **flow block 内の一般 EOS ヤコビアン
  (H2) か SST/flow 剛性 (H3)** が律速の可能性が高い(案C は `A_QQ` を触らないため効かない)。

### freeze 切り分け + 発散セル特定 (2026-06-20, 案C 後の追加診断)

案C が効かなかった原因を freeze 試験で切り分けた (case/28 TP-SST, run_0038 発達場, cfl=2; mode-1 は step50 発散)。

| 試験 | 初 NaN | 示唆 |
| --- | --- | --- |
| SST off (laminar, μ_t=0) | step 40 | SST でない (μ_t 無いと逆に早い) |
| μ_t freeze (`FORGE_FREEZE_TURB=1`, k-ω 解かず) | step 50 | **H3 棄却** |
| species freeze (`FORGE_FREEZE_SPECIES=1`, ρY 凍結) | step 51 | **H1 棄却** (組成完全凍結でも発散) |
| both freeze | step 52 | — |

全て発散 → **律速は flow 5×5 block / 一般 EOS / 境界 / DPLUR 本体**。species 凍結でも不安定なので
**案A は空振り (ユーザ判定表どおり)。着手しない。**

**発散セルの特定** (`detectNaN:1`, frozen species, cfl=2): 'ro' が step52 で非有限。NaN セル 50 個は全て
**軸近傍 (r≈0.04–0.6mm = 最内セル列) かつ下流 x/D≈18.7–19.4**。max cfl が step51 で **2e9 に暴走**。
壁でも入口/出口境界でも He/air contact (r≈5mm) でもない → **軸対称 near-axis (r→0) 領域**。

その場所の物理 (発達場): **超音速 He コア on axis** = Mach~1.8, T~150K (冷), **sonic~608 m/s (He-rich で高い)**,
最内セル半径 r=0.04mm。停滞域ではない。CPG (cfl2 安定) 比較: 同位置で **sonic~437 m/s** (T~475K に誤熱化,
γ=1.4)。**TP は物理的に正しい He-rich 冷コアゆえ音速が ~40% 高く**、微小半径の軸セルで局所安定限界を
先に超える → cfl~1 で発散。CPG は音速が低いので cfl2 まで保つ。

**最有力仮説 (要確認, 確定ではない)**: 多成分 TP の CFL 律速は **軸対称 near-axis (r→0) の flow-block
安定性** (既知の forge 弱点, [[architecture-axisym-axis-singularity]] / [[mixed-precision-axisym-refuted]] の
block-DPLUR 近軸ジオメトリ)。**species coupling でも SST でも EOS Jacobian の正しさでもない** (freeze で確定)。

ただし「TP の高音速 (608 vs CPG 437, 1.39倍) が CFL を 2→1 に半減させた」までは**断定しない**。Δτ が正しく
`Δτ_i = CFL·V_i / [Σ_f(|U_n|+c)S_f + …]` で計算されていれば音速上昇は Δτ を自動で縮めるので、1.39倍では
2→1 の説明として弱い。**本質はおそらく「setDT が見積もる近軸スペクトル半径」と「実際の軸対称残差・DPLUR
演算子の強さ」の幾何学的不整合**。`max cfl=2e9` (設定 CFL=2 に対し) はその傍証 (単なる高音速では出ない)。
TP の高音速は不整合を低 CFL で露呈させる増幅器、という位置づけ。

**次に点検すべき場所 (near-axis, ②; 別 commit で)**:
1. **軸対称残差と時間刻みの幾何一貫性**: RHS 流束 (`r_f S_f`)・時間項 (`r_i V_i`)・setDT (planar?)・DPLUR 対角の
   重みが同一の軸対称幾何定義か。`D_i ∼ V_i^axi/Δτ_i·I + Σ A_f^+ S_f^axi + J_axis·V_i^axi` と
   `Δτ_i = CFL·V_i^axi/[Σ(|U_n|+c)S_f^axi + λ_axis V_i^axi]` が同じ幾何で揃っているか。
2. **軸面 r=0 の扱い**: 最内セル中心 r_c≈0.04mm → 軸面はほぼ r_f=0。面積をゼロ扱いしているか、planar 面積の
   まま流束/ghost 法線流束を二重計上していないか、圧力・粘性の軸対称ソースと相殺しているか。
3. **1/r ソースが LHS に入っているか**: 軸対称運動量/エネルギーの幾何ソース `S_axis(Q)/r` が RHS にあるのに
   DPLUR 対角に `∂S_axis/∂Q` が抜けていると r→0 で過大な陽的項 → 「低 CFL なら動く・高 CFL で近軸から崩れる」
   = 今回の症状と一致。径方向運動量・エネルギー・軸対称粘性ソース (swirl あれば周方向) を見る。完全な
   axis-source Jacobian でなくても `D_i ← D_i + α_i I` の対角近似で安定性が変わるか診断できる (案C と同手法)。
4. **max cfl=2e9 の意味**: NaN 直前 step の軸セルで `V_i, r_i, ΣS_f, Σ(|U_n|+c)S_f, Δτ_i, V_i/Δτ_i` を出し、
   `V_i/Δτ_i` が DPLUR 対角の流束スペクトル半径と同じ桁か比較。Δτ 異常増大/スペクトル半径ゼロ・負/体積・面積
   不正のいずれかを特定する。

### 機構確定 (発散セル時系列トレース) + α 診断 (2026-06-20)

近軸セル (r=0.04mm, x/D19) を step45→51 で追跡 (frozen species, cfl=2):
- **半径方向運動量 Uy が振動・増大** (−16→+103→+215→+300 m/s, step45-48) → **密度崩壊** (ro 0.46→0.28→step50 で負)
  → catastrophe。駆動 = 軸対称圧力ソース `res_roUy += P·A_pl` の半径加速度 **P/(ρr)≈7e9 m/s²**
  (dt_local=1.2e-7 で ~800 m/s/step のキック)。`max cfl=2e9` は ro<0 後の garbage (発散後値, 幾何 dt エラーでない)。
- **TP が CPG より早く死ぬのは音速でなく密度**: 物理的に正しい冷 He コア ρ~0.42 < CPG 誤熱化 0.68 →
  P/(ρr) が ~1.6倍。音速は setDT で dt を縮めるので逆に安定化側 (先の「音速 40% 説」は撤回)。
- 根本: dt_local が **近軸半径音響タイムスケール r/c≈6.7e-8 を超過**。block 対角の軸対称ソース Jacobian は
  ∝(γ-1)u_y で小さく不足。revolved 軸面積 (r_f·S→0) が半径音響スペクトル半径を落としている。

**α 診断 (env `FORGE_AXIS_DIAG_ALPHA`, timeIntegration_d.cu, 既定0=不変)**: roUy 対角に `α·A_planar·c` を追加:

| CFL | baseline | α=1 | α=5/10 |
| --- | --- | --- | --- |
| 2 | step50 死 | step278 | **安定** (limit cycle) |
| 4 | step1 死 | step115 | **安定** (6000步, ro>0/He コア保持) |
| 8 | step9 死 | — | step160 (不足) |

→ **TP-SST 安定 CFL が ~1 → 4+** に上昇。機構を確定。ただし crude 対角は roUy を per-equation 過減衰し
残差フロアが高め (cfl4 で rms_ro~4e-6 vs cfl1 の 3e-7)。

### production 形 = setDT 軸スペクトル半径 (実装済 commit, config `axisTimestepBeta`)

**実装** (`time.deltaT.axisTimestepBeta`, 既定0=ビット不変): `setCFL_cell_d` で軸対称セルの cfl に
`dt·β·(|u_r|+c)·A_planar/V` を加算 (face 項と同無次元化)。V=r·A_planar より近軸で Δτ∝CFL·r/(|u_r|+c)。
DPLUR 時間項 v/Δτ が全保存式で自動的に強まる (per-equation 過減衰なし)。**結果** (run_axdt_*, α off):

| CFL | β=1 | β=2 |
| --- | --- | --- |
| 2 | **安定** | — |
| 4 | step88 死 | **安定** (6000步級, physical) |

β は CFL に応じスケール (Θ=cfl_pseudo/β 近軸; cfl4 は β≈2)。検証: 近軸 **Θ=Δτ(|u_r|+c)/r ≈ 0.5–1.3 (mean 1.02)**
= 設計目標 (半径音響 CFL≈1)。**遠方セル dt 不変** (near-axis dt~6.9e-8≈r/c vs far-field 5.6e-7, 軸項は 1/r で
近軸のみ効く)。残差フロアは α 対角版とほぼ同 (~4e-6 @ cfl4 = 高 CFL の limit-cycle floor, fix 由来でない)。
CPG β=0 既定はビット不変 (gating)・CPG cfl2 安定を確認。field 物理 (ro>0, He コア保持, ΣY=1)。

> production 既定: `axisTimestepBeta=0` (opt-in)。軸対称ケースで CFL を上げたいとき β≈cfl_pseudo/2 を設定。
> α 対角診断 (`FORGE_AXIS_DIAG_ALPHA`) は既定0/デバッグ専用に残置 (production には使わない)。

### 最終ステータス (2026-06-20): production 確定 + contact 混合層モードを別 issue に分離

**production 形を additive fixed-β に確定** (`time.deltaT.axisTimestepBeta`, 既定0)。`axisAcousticCFLMax`
(min-clamp / 目標Θ) は試したが棄却: min-clamp は最内セルしか触れず CFL≳4 で不足、目標Θ形 (係数=cfl_pseudo/C)
は face 項のぶん実効 Θ が C/2 になり C=1 で過減衰発散。**設定値と実効 Θ が一対一でない**ため名称を約束しない。
docs に「`axisTimestepBeta` は数値スペクトル半径の重み・Θ≈CFL/β は厳密保証でない・β≈CFL/2 は case/28 経験則」を明記
([docs/axisymmetric/theory.md](../../docs/axisymmetric/theory.md) / [implementation.md](../../docs/axisymmetric/implementation.md))。

**CFL4 長時間 + down-test (検証完了)**: β=2 で CFL4 安定 (8000步, 残差 5.4e-6→4.5e-6 緩やか低下)。
質量流量 in/out 比 1.362 (baseline 1.363 と一致)、ΣY=1, He コア保持。**CFL4→1 down-test で残差再低下**
(rms_ro 4.5e-6→5e-7) し近軸 |u_r| 174→4.4・軸上 P スパイク消失で **baseline 低 CFL 解に一致** → 高 CFL 残差床は
擬似時間 limit cycle、**定常解は不変**。

**catastrophic divergence と residual floor は別機構 (確定)**:
- catastrophic = near-axis 半径運動量不安定 → `axisTimestepBeta` で CFL4 まで解消。
- residual floor (~4e-6 @ CFL4) = **He/air contact 混合層モード**。スナップショット差分で contact の
  |Δu_r|・|ΔY_He|・|Δp| が軸近傍より大、最大は軸と組成勾配が交わる下流。
- **Test 4 (pseudo-CFL 依存性) で数値 limit cycle と確定**: contact 振幅が CFL に強依存
  (A_YHe 1.7e-4→9.6e-3→2.5e-2 @ cfl1/2/4, ~100×)。CFL1 でほぼ消滅 → **物理非定常でなく定常反復の数値 limit cycle**。
- 化学種移流は既に 1 次風上 (`scalar_advection_first_order_d`) なので **species 高次再構成ではない** (Test1 無効)。
  flow 2次再構成と低次 LHS の不整合 (defect-correction limit cycle) / limiter chatter が候補。flow 1次化は
  別の発散を起こし切り分け不能だった (transient 再平衡)。

**case/28 への影響 (engineering 量)**: massflux・軸上 p,T・組成場は baseline 低 CFL 解と一致 (down-test)。
field 物理。**CFL8 は現 production 要件としない** (まだ near-axis inlet で発散)。

**案A (5+nSpecies full block) には着手しない**: species/SST freeze でも catastrophic 継続・案C JVP は FD 一致だが
CFL 無影響・catastrophic は軸 flow モード・floor は contact モード、が確定済み。contact が `A_YQ`/同期 LHS 不足と
明示されたときのみ再判断。

**別 issue として分離**: 「多成分 TP 高 CFL の contact 混合層 数値 limit cycle」(flow 高次再構成 / limiter / 高次RHS-低次LHS)。

**merge 可否**: near-axis 安定化 (`axisTimestepBeta`) は production 候補として確定・既定0でビット不変・docs 整備済。
merge 可。contact 混合層 limit cycle は独立 issue として継続。

#### (旧メモ) production 形 = setDT 軸スペクトル半径 (設計)

人工対角でなく、**通常のスペクトル半径に軸項を足す** (単一式の人工項→全保存式整合の擬似時間スケールへ):

  λ_i^total = λ_i^faces + λ_i^axis,  λ_i^axis = β_axis·(|u_r|+c)·A_planar
  Δτ_i = CFL·V_i^axi / (λ_i^faces + β_axis·(|u_r|+c)·A_planar)

近軸では V_i^axi = r_i·A_planar (per-radian) なので Δτ_i ∝ CFL·r_i/(|u_r|+c) が自然に出る。利点: 通常セルは
軸項が相対的に小さく不変、近軸だけ制限、CFL の意味維持、face spectral radius と同構造。
**DPLUR 時間項が V_i/Δτ_i·I なら setDT 短縮で LHS 対角が全式について自動で強まる** (α 対角の per-equation
過減衰を避けられる)。**二重計上注意**: setDT の軸補正と DPLUR 対角の軸補正を同時に入れない (DPLUR が
setDT とは別の独立スペクトル半径を使っていないか確認)。
- 診断値 `Θ_i = Δτ_i(|u_r|+c)/r_i` を出し、近軸で 1 を大きく超えないことを確認。
- 検証: 対角補強版 vs setDT 版を CFL sweep (1,2,4,8) で比較 — 初 NaN step / 最小 ρ / max|u_r| / Θ /
  rms_ro / rms_roUy / 残差フロア / massflux / 軸上 P,T / He コア半径対称性。
- 期待: setDT 版で cfl2,4 安定 & ρu_r 振動消失 & α 対角版より残差フロア低 & 遠方セル dt/収束率ほぼ不変
  → production 候補。config 化 (β_axis)、α 診断は既定0/デバッグ専用に。

**次の手 (要判断)**:
- **本筋**: 軸対称 near-axis の flow-block 安定化 (局所 dt 制限 / 近軸 block-DPLUR ジオメトリ補正)。
  既存の [[architecture-axisym-axis-singularity]] 系の課題に合流する。TP は音速が高いぶん厳しいだけ。
- 確認したければ: ① CPG を cfl↑ して同じ軸位置で発散するか (geometry 起因の裏取り)、② setDT の近軸
  局所 dt と block-DPLUR 対角の r-weight を軸セルで点検、③ 5×5 ヤコビアン FD は優先度低 (局所化済み)。
- **案A/案C 系 (species-in-DPLUR) は CFL 律速ではないので棚上げ。** 案C コードは診断モードとして残置。

#### (以下は案C 着手前の方針メモ。上記診断で species coupling は CFL 律速でないと判明したため棚上げ。)
**当初の次の手 (要判断)**:
- **案A (full LHS coupling)**: 化学種を 5×5→`5+nSpecies` block に同梱し `A_QY` を **LHS** に入れる。
  deferred でないので安定性を変えうる。ただし H2 が真因なら案A でも上がらない可能性。実装重い。
- **先に H2/H3 を安く切り分け**: ① TP 多成分を SST off(層流)で cfl 上限を測る(H3 切り分け)、
  ② block の一般 EOS ヤコビアン `eos_split_jacobian_general_closed` の音速・固有値を単成分 TP と比較、
  ③ 最初の発散セル(`detectNaN:1`)が contact/組成勾配最大か(H1/H2 切り分け)。
- 推奨: **案A に着手する前に ②③(+①)で H2 を確認**。H2 が真因なら修正対象は block EOS ヤコビアンになる。

> 案C のコード (`speciesImplicitCoupling: 2`) は正しく動く診断モード/案A の足場 (JVP・射影・最終 commit
> は流用可) として残置。本番の既定は引き続き mode 0/1。

---

## 実装方針 (II): 結合ブロック — 段階的に (案C 実装済・上記結果を反映して案A/H2 へ)

`reordering / sweep の同時化だけでは無意味`(ブロック対角系の Jacobi/GS は順序に依らず同じ不動点)。
**LHS に十字を入れる**ことが要件。以下の 2 段階で進める。

### Step 1 = 案C: EOS クロス項だけの partial coupling(まず実装、診断兼用)

- 線形系を `5 + nSpecies` に拡張するが、**埋めるのは支配項だけ**:
  - **エネルギー行 × 化学種列** `∂R_energy/∂(ρY_s)`: 組成変化が `P=ρR_mix(Y)T` と被移流エンタルピー
    `h_mix(T,Y)` を通じてエネルギー対流流束へ与える感度。これが**落としている支配的な向き(Y→e)**。
  - 化学種行は既存の**対流対角(contact 波速 `V` のスペクトル半径)+ 緩和整合**を流用(固有ベクトル
    分解はやり直さない)。`∂R_Y/∂Qflow`(速度・密度依存)は二次的、必要なら upwind 対角で近似。
- **狙い**: 5×5 の固有系を温存したまま、`δ(ρe)` と `δ(ρY)` を同一 sweep で整合させ温度ジャンプを消す。
  実装は block サイズ拡張 + エネルギー行へのクロス項追加に限定でき、案A より軽い。
- **これ自体が H1 の決定実験**: 案C で case/28 TP-SST の上限が上がれば「結合欠如が律速」を確証。

### Step 2 = 案A: full coupling((5+nSpecies) 真の結合、恒久)

- 案C で不足(なお上限が低い/プラトー)な場合に進む。`block_dplur`(`timeIntegration_d.cu`)を
  `N=5+nSpecies` に一般化。化学種対流は同一対流ヤコビアンの追加行(contact 波速)として同緩和で進め、
  `roe↔roY` 両方向のクロス項を入れる。**`lambda[5]`/`[5][5]` 等の 5×5 ハードコードを可変長 block に
  全面テンプレ化**するのが重い部分。

> 推奨順: 案C を最小実装 → case/28 と case/16 で上限を測る → 足りなければ案A。案C で十分なら案A は不要。

## 現状アーキテクチャ (正確なコードポインタ)

- **オーケストレーション**: `main.cpp` `implicitNonlinearUpdate`(≈L900)。1 アウター擬似ステップ =
  `assembleResidual` **1 回**(組成・流れ凍結)→ `blockDPLURSolve`(flow を `nStepInner` 回 Jacobi sweep,
  `cuda_forge/timeIntegration_d.cu` の `namespace block_dplur`)→ `applyBlockImplicitCorrection`(flow 5 変数
  commit)→ `applySSTPointImplicit`(SST 別 segregated, 生産項 lagged)→ `speciesImplicitDPLURSolve`
  (化学種 scalar-DPLUR)→ `speciesRenormalize`(`ρY_s≥0`,`ΣρY_s=ρ`)→ `speciesPrimitive`。**3 解法が
  完全に別連立**で off-diagonal 無し。
- **block-DPLUR コア**: `timeIntegration_d.cu` `block_dplur`。**5 変数ハードコード**(`ro,roUx,roUy,roUz,roe`)。
  対流は真の分離 A±(`build_jacobian_split`, precond=2)、EOS ヤコビアン `∂P/∂(保存量)`(可変 γ)は
  `eos_split_jacobian_general_closed`(`eos_jacobian_d.cuh`)。**block 内では組成 `Y` 固定**。
- **EOS 反転**: `cuda_forge/dependentVariables_d.cu`(`thermalMethod==2` 分岐, `thermo_T_from_e`)。
  ここで `T=T(e,Y)` を解くので、`e`/`Y` 不整合が温度ジャンプとして顕在化する。クロス項の素は
  `∂T/∂Y`, `∂P/∂Y`, `∂h/∂Y`(NASA-9, `thermo_d.cuh` の混合則)。
- **化学種残差/解法**: `cuda_forge/speciesTransport_d.cu`(`res_roY{s}`/`transport_diag_Y{s}`/`src_jac_Y{s}` を
  `assembleResidual` で確定、`speciesImplicitDPLURSolve_d_wrapper`)。device ポインタ配列
  `g_roY_dev`/`g_resroY_dev` 等。`registerSpecies`(`variables.cpp`)が `roY{s},Y{s},...` を動的登録。
- **SST 緩和分離**: `cfg.implicitRelaxSST`(既定 -1→`implicitRelax`)。SST 寄与は H3 として別途切り分け。

## 検証計画

- **主ドライバ = `case/28.cutler_coaxial_jet`(He/O2/N2, 非反応, sub-200K 無し, 生成熱ゼロ)**: 顕熱
  コントラスト機構を**最もきれいに**切り分けられる。元データ: `run_0037`(CPG-SST)/`run_0038`
  (TP-SST cfl0.5, 20000步)/`run_0039`(TP-SST cfl1)。継続は `run_0038/res_20000.h5` から。
  - **比較軸**: `cfl_pseudo` sweep(0.5,1,1.5,2,4)。**目標 = TP-SST で cfl>1 安定 & 残差プラトー割れ**。
  - **回帰として CPG-SST**(`run_0037`)が**不変**(上限・解とも)であること(EOS 非依存経路を壊さない)。
- **副ドライバ = `case/16.nozzle_wys`(hot N2+H2O 1%, >200K)**: 生成エンタルピー機構側。`thermoHrefTemp`
  on/off と案C/案A の効果を交差で見る(結合が効けば offset 無しでも上限が上がるはず)。run_0081 系/A-B sweep run_0101〜0109。
- **判定**: `solver_density_cuda/tools/check_convergence.py <run>` の VERDICT を必ず貼る。
  NaN/Inf・`ΣY=1`・`Y∈[0,1]`・T>200K(wys)を確認。最初の発散セルは `detectNaN:1` で特定。
- **回帰 (必須)**: ① 単成分 TP・CPG がビット不変(`nSpecies<2` で従来経路)、② case/28 CPG が不変、
  ③ `case/08.bump`/flat_plate など既存 implicit が不変。

## 進め方・落とし穴

- **開発フロー (AGENTS.md)**: 数値スキーム変更 → 4 ステップ厳守。①`docs/time_integration/`(block へのクロス
  項)theory/implementation 更新、②`docs/thermophysics/`(EOS の `∂T,∂P,∂h/∂Y`)更新、③本ファイル/
  `.github/plans/thermophysics-multicomponent-tpgas.md` §10 の plan を `status` 込みで更新、④実装。
- **stale build trap** ([[stale-build-struct-layout-trap]]): block サイズや `solverConfig.hpp` 変更後は
  `ninja -C .build-native/relwithdebinfo -t clean && ninja` で **full rebuild**(差分は step0 NaN 凍結)。
- **native build**: `.build-native/relwithdebinfo/forge`(nvcc 12, sm_86)。`tools/build_native_wsl.sh`。
- **CFL の定義**: 定常は `cfl_pseudo` が実効(`setDT_d.cu`)。`cfl` は表示用。`guide/solver-settings.md`。
- **収束報告の規律** ([[convergence-check-discipline]]): `rms_ro` 単独で判断しない・全残差列・ツール経由。
- **run 索引**: `case/28.cutler_coaxial_jet/README.md` の「## 計算 run 一覧」を同期。新 run は
  `run_NNNN_<slug>` 連番(現在 run_0039 まで使用、0040–0044 は CFL 上限探索で破棄済み)。
- 着手前に `git status`/`git diff` で無関係 WIP を巻き込まないこと。

## 参考

- 分析プロンプト(H1–H5 切り分け): `.github/plans/cutler-tp-multispecies-cfl-analysis-prompt.md`
- plan: `.github/plans/thermophysics-multicomponent-tpgas.md`(§10 多成分 implicit 結合)
- docs: `docs/thermophysics/theory.md`(sensible-enthalpy datum)/`implementation.md`, `docs/time_integration/`
- case README: `case/28.cutler_coaxial_jet/README.md`(CPG vs TP / CFL 上限節), `case/16.nozzle_wys/README.md`
- memory: [[cutler-cpg-vs-tp-dplur-sst]] / [[wys-tp-divergence-is-cold-not-multispecies]] /
  [[implicit-blockdplur-config]] / [[forge-implicit-inner-cfl-tuning]]

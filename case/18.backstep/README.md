# case/18.backstep — 後退段 (Backward-Facing Step)

ステップ後方の剥離・再循環・再付着を扱う基本剥離流ケース。RANS/LES の速度評価
(`SPEEDUP_CAMPAIGN.md`) と **SST-DES (DDES) の機能検証 (T1-B)** に用いる。

## メッシュ

| ファイル | 次元 | 要素数 | 用途 |
| --- | --- | --- | --- |
| `mesh/backstep_2d.msh` | 2D (厚さ0.1の1層hex押し出し=実体3D) | 32,712 | RANS・速度評価 (スパンなし) |
| `mesh/backstep.msh` | 3D (スパン 4H, 80 層) | 923,200 (→ 857,600 cells) | **DDES**: スパン方向解像が必要 |
| `mesh/backstep_planar.msh` | 2D 平面 (押し出しなし, `gmsh -2`) | 10,997 nodes | **node-centered 検証用** (median-dual は平面2Dのみ; 押し出しスラブは is3D=true で未対応) |

- 物理面: `inlet`(1) / `outlet`(2) / `top`(3, wall) / `bot`(4, wall) / `side2`(5) / `side1`(6)。
- DDES では `side1`/`side2` を **周期 BC** にしてスパン方向の解像乱流を許す (2D・slip では RANS 同然)。
- 3D メッシュ品質: `check_mesh_quality.py` VERDICT **PASS** (AR max 2.3, skew 0.333)。

## 流れ条件

低マッハ (M≈0.1 級) の非圧縮的後退段。NASA/TMR 2DBFS (Driver & Seegmiller 1985) を
定量参照とする (Re_H≈36,000, 再付着 x_R/H=6.26±0.10)。定量比較は Phase 1.5 以降。

## 計算 run 一覧

> 速度最適化キャンペーン (`run_0001`〜`run_0034`) の詳細は [`SPEEDUP_CAMPAIGN.md`](SPEEDUP_CAMPAIGN.md)。
> ここでは代表 run と DDES 検証 run を示す。

| run_* | 目的・主要設定 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_rans_implicit` | 2D RANS SST + block-DPLUR 陰解 (`backstep_2d.h5`) | 収束 RANS 場 (再循環泡)。DDES restart 元 | active |
| `run_0002_slau_rans_rk3` | 2D RANS RK3 陽解 | ⚠️ **元から発散** (自身の res_2000/4000/6000 が全 NaN, 06-07 作成時から)。RK3 explicit は backstep 低マッハ RANS で不安定。「速度評価基準」表記は誤り | broken (記録) |
| `run_0007_slau_2ndup_baseline` 〜 `run_0034_slau_K5_final` | 速度最適化キャンペーン (2D) | カーネル最適化の段階計測 ([SPEEDUP_CAMPAIGN.md](SPEEDUP_CAMPAIGN.md)) | ref |
| `run_0035_des3d_ddes` | **SST-DES T1-B 機能検証 (3D)**: `backstep.h5`(857k, 周期スパン)、`DESmode:1`、2D RANS 場から `interp_field` cross-mesh restart、500 step implicit | **NaN なし**。f_d: 付着近壁 BL≈0.02 (frac<0.1=0.96)、剥離せん断層 (x∈[3,9))≈0.67=LES 活性、再循環泡内部 (Ux<0,高νt)≈0=RANS、再付着~6H。Δmax≈0.13、rd∈[5e-4,10]。`residual_history.png` | active |
| `run_0036_node_rans_2d` | **node-centered (median-dual) 検証**: 新規平面メッシュ `mesh/backstep_planar.geo` (10997 nodes, `gmsh -2`, push 無し)、`discretization:node`+`gradLSQ:1`+`nodeWallViscGradFlux:1`、cell 収束場から最近傍 cross-mesh restart、RANS SST RK3 cfl0.3 | **DIVERGED @step18**。roOmega が step18 で 9.4e25→NaN。`residual_history.png` | 破棄予定 (記録) |
| `run_0037_node_laminar_2d` | run_0036 を **層流** (LESorRANS:0) に。turbulence 固有か切り分け | **同じく DIVERGED @step18** → turbulence 非依存 | 破棄予定 (記録) |
| `run_0038_node_lam_nolsq` | 層流 + `gradLSQ:0`+`nodeWallViscGradFlux:0` (検証済み M3 経路) | **同じく @step18** → 未検証フラグでも convMethod でもない | 破棄予定 (記録) |
| `run_0039_node_lam_diag` | 発散直前 (16 step, out interval 2) ダンプで空間切り分け | **種=段差直後 再循環域 (x≈2.26,y≈0.56-0.62) の隣接 CV odd-even (チェッカーボード) ro 振動**が 2step ごと倍々成長→step18 爆発。種域 Mach≈0.005-0.03 (淀み) | 破棄予定 (記録) |
| `run_0040_node_lam_dirichlet` | 働く壁レシピ `nodeWallDirichlet=1`+`bndFirstOrder=1`+cfl0.1 を適用 (case/29 由来) | **同じく DIVERGED @step18** → 壁レシピでは直らない (種は壁でなく内部淀み域) | 破棄予定 (記録) |
| `run_0041_node_lam_midfx` | + `nodeMidpointFx=1` (median-dual checkerboard 低減フラグ) | **同じく DIVERGED @step18** | 破棄予定 (記録) |
| `run_0042_node_lam_precond` | + `lowMachPrecond=2` | **step0→step1 即 NaN** (前処理が node で破綻) | 破棄予定 (記録) |
| `run_0043_cell_lam_planar` / `run_0045_cell_rans_planar` | **対照: 同じ planar 2D を cell で** (層流/RANS) | **どちらも DIVERGED** (step26/27、node と同じ再循環域 x≈2.46) → 発散は node 固有でない | 破棄予定 (記録) |
| `run_0047_rerun0002_rk3` / `run_0048_rerun0002_oldbuild` | run_0002 を現 build / 旧 build(06-14) で再計算 | **両方 step28 DIVERGED** (収束場 restart でも) → build リグレッションでなく RK3 explicit が元から不安定 | 破棄予定 (記録) |
| `run_0049_node_rans_implicit` | **node + implicit (block-DPLUR) 初検証**: run_0036 と全設定同一で積分法のみ RK3→`timeIntegration:11`(blockDPLUR=1, implicitRelax=0.5, cfl_pseudo=1)。5000 step | **step18 の壁を初突破・NaN なし・全 7 列 falling** (ro 2.5桁/roe 2.6桁低下)。explicit RK3 発散は node 固有でなく積分法の交絡だったと確定 | ref (PoC) |
| `run_0050_node_rans_imp_long` | run_0049 を **元 IC から 30000 step** に延長 (chain restart 回避)、cfl_pseudo=1 維持 | `check_convergence` VERDICT=**`NOT CONVERGED (stalled/plateau — needs scheme change)`**: rms_ro 3.9桁/roUx 3.3/roUy 3.3/roe 4.2/roK 2.2 桁低下するが **rms_roOmega が 1桁プラトー** (壁近傍 ω ソースの SST 既知挙動)。machmax/pmax は `check_quasisteady` で `STEADY`、再循環泡あり (Ux<0 115node, Ux_min=-2.17)、場は物理的 (T~288K, ro>0, NaN無)。**発散せず平均流は落ちるが omega 律速で未収束**。`CONVERGENCE_VERDICT.txt`/`residual_history.png` | active |
| `run_0051_node_imp_muscl` | run_0050 から **`gradLSQ`/`nodeWallViscGradFlux` を OFF**・`convMethod` 1次→**2(MUSCL)** (場の過拡散が 1次風上由来か検証) | VERDICT=**`NOT CONVERGED (plateau)`** (run_0050 同様 omega 律速)。**場は 1次とほぼ不変**: vis_turb/vis_lam 場平均 **734** (max 4304)、x_R≈4.75H、Ux[-2.2,34]。**過拡散の主因は 1次風上でなく SST の過大乱流粘性**と判明 (MUSCL でほぼ変わらず) | 記録 (診断) |
| `run_0052_node_imp_laminar` | run_0051 を **層流 (`LESorRANS:0`)** に (乱流粘性が犯人か対照) | VERDICT=**`NOT CONVERGED (rising/divergent)`**: vis_turb=0 にすると node 低マッハ不安定が復帰し残差 rising、Ux が入口の倍 68m/s へオーバーシュート (非物理)、x_R=22.8H。**乱流粘性 (~700×) が場を過拡散させると同時に node 不安定を覆い隠していた**ことを裏付け | 記録 (対照) |
| `run_0053_node_imp_autowall` | run_0051 + **`wallTreatmentSST:1` (automatic wall function)** (y+~27 で low-Re 前提が崩れるため) | crash 前は **明確に改善** (step4000 で vis_t/l 平均 237・x_R≈8.9、low-Re の 577 より低乱流) が **step5733 で DIVERGED**: omega 残差が 1 step で 21→6.4e5 と stiff 爆発 (cfl_pseudo=1, nStepInner=10 が抑えきれず) | 記録 (診断) |
| `run_0054_node_autowall_cfl04` | run_0053 の crash 対策で **cfl_pseudo 1.0→0.4・nStepInner→20** | **逆効果: step891 で DIVERGED** (5733 より早い)。node 低マッハは CFL を下げると stiffness が速く成長 (memory 既知) | 破棄 (記録) |
| `run_0055_node_autowall_in5` | run_0053 から **nStepInner 10→5** に (cfl_pseudo=1.0 維持) | **crash 回避・30000 step 完走 NaN無**。inner sweep 過多が stiff omega を overshoot させていた。VERDICT=`NOT CONVERGED (still converging)`、**omega 残差が plateau→falling に改善**、vis_t/l 平均 252 | ref→0056 |
| `run_0056_node_autowall_long` | run_0055 を **100000 step** に延長 | **定常解到達** (x_R step45k 以降一定、`check_quasisteady` machmax/pmax `STEADY`、残差 3.5–3.8桁 plateau=limit cycle)。⚠️**x_R の壁バンド法 9.59 は誤り**: `nodeWallDirichlet=0` で壁速度が負ドリフトしたのを拾ったアーティファクト。第一内層法では x_R=6.71 (run_0059 と同値)。`residual_history.png`/`CONVERGENCE_VERDICT.txt` | active |
| `run_0057_cell_autowall_imp` | **クロスチェック A**: 同一 planar mesh を **cell** で (node 0056 と数値設定一致: implicit+MUSCL+autowall+nStepInner5) | **完走・omega 残差 4.6桁** (node の 1.5桁より良収束)。x_R=**7.95** (第一層法)、vis_t/l 198。VERDICT=`NOT CONVERGED (plateau)` だが派生量定常。`CONVERGENCE_VERDICT.txt` | active |
| `run_0058_su2_sst` | **クロスチェック B**: SU2 v8.5 SST、同一 `.geo` から `.su2` 生成、同条件 (Pt/Tt/Ps/定数粘性/壁関数) | **forge cell と一致**: x_R=**7.89** (Cf 符号反転=第一層法とも一致)、eddy/lam 239、y+~62。残差は低マッハ剥離で limit-cycle プラトー (rms_density -2.95)。**forge cell≈SU2 (7.95 vs 7.89, 差0.8%) で forge-SST 較正は正しいと確定**。`case.cfg`/`su2.log`/`vol_solution.vtu` | active |
| `run_0059_node_autowall_dirichlet` | run_0056 + **`nodeWallDirichlet:1`** (壁速度ドリフト修正; user 指摘) | **壁ノード max\|Ux\|=0 (厳密ピン)**。x_R=**6.71** (定常, step75k 以降一定)、Ux_min=-8.05 (SU2 -8.8 と同等)、vis_t/l 424。VERDICT=`NOT CONVERGED (plateau)` だが派生量定常。`residual_history.png`/`CONVERGENCE_VERDICT.txt` | active |
| `run_0060_node_kbudget` / `run_0061_cell_kbudget` | **k 収支診断** (収束場から 200 step、`wf_pk`/`Pk_diag`/`src_jac_k` 出力): node 過剰乱流の真因切り分け | **確定**: 段差壁の角で **cell は第一セル列が全て wf_pk≥0 (壁関数処理→生産≈0→k=3-5)**、**node は壁ノード1点のみ処理・第一内層は wf_pk=-1 未処理→標準生産 Pk=166-228→k=52 暴走**。真因=**node 壁関数の近壁カバレッジ不足** (omega/decouple/strain でない)。[notes/investigations/turbulence-node-wall-function-coverage.md](../../notes/investigations/turbulence-node-wall-function-coverage.md) | 記録 (診断) |
| `run_0062_node_wfcov` / `run_0064_node_wfcov2` | **修正試行 (失敗)**: wf_pk と ω ピンの両方を第一内層 irep へ移設 | **悪化** (k 52→246-350): ① 凹コーナーで irep 共有→ω ピン値 race、② 壁ノード ω アンカ喪失→ω 崩壊。**ω ピンは irep へ移してはいけない**と判明 | 破棄 (記録) |
| `run_0065_node_wfcov3` / `run_0067_node_wfcov_long` | **修正 (成功)**: **生産 wf_pk のみ壁+第一内層に被覆、ω ピン/decouple は壁ノード据え置き** ([plan](../../plans/accepted/turbulence-node-wall-function-coverage.md)) | **段差コーナー vis_t/l 6789→77・k 52→3** (cell 同等)、場平均 424→**207** (cell 198)、**x_R 6.71→7.63** (cell 7.95/SU2 7.89 に接近)、machmax/pmax `STEADY`。残: 再付着で μt ピーク (vt/l~5800 vs cell 990、局所)。`CONVERGENCE_VERDICT.txt` | active |
| `run_0066_cell_invar` | **cell ビット不変確認**: 修正ビルドで cell を 2000 step | vis_t/l mean **197.71**・max **994** が baseline 一致 → **cell 経路は無影響** | 記録 (回帰) |
| `run_0068_node_kdirichlet` | **k Dirichlet 試行** (SU2 SetTurbVars_WF 流): 第一内層ノードで k=ω_w·μt,wall/ρ をハード Dirichlet | **再付着 μt ピーク除去** (vis_t/l max 5769→894=cell 994 並、再付着 vt/l 5834→5、場平均 207→183) **だが x_R 7.63→8.67 と伸びる** (非平衡再付着で u_τ→0→k_wf→0 と乱流抑えすぎ)。トレードオフ | 記録 (診断) |
| `run_0069_node_kwf_off` / `run_0069_node_kwf_on` | **`nodeKwfDirichlet` フラグ検証** (既定 0) | OFF=**prod-fix と完全一致** (mean 206.8/max 5769/x_R 7.63, k-pin 0 node)。ON=k-Dirichlet (max 893/x_R 8.54, k-pin 472 node)。フラグで opt-in 化、既定は prod-fix | 記録 (回帰) |
| `run_0071_node_lmp1_eps15` | **P/ρ 市松の打ち手検証**: node + `lowMachPrecond=1` (cfl1)、`run_0067` 収束場から restart | **DIVERGED (NaN)**: step0 rms_ro 0.0126 に跳ね倍々成長。RHS のみ前処理が node で不安定 | 破棄 (記録) |
| `run_0072_node_ctrl_lmp0` / `run_0074_cell_ctrl_lmp0` | control (lmp0) 10k step、収束場 restart | 残差プラトー据え置き (rms_roUx ~1e-2)。市松 P odd-even ~70 Pa 維持。A/B 基準 | ref |
| `run_0073_cell_lmp1_eps15` | cell + `lowMachPrecond=1` cfl1 | **DIVERGED (NaN)** (node と同様、RHS/LHS 不整合) | 破棄 (記録) |
| `run_0075_cell_lmp1_cfl03` | cell + `lowMachPrecond=1` **cfl0.3** | NaN回避だが高残差停滞 (rms_roUx 0.41)。P odd-even 78→21 Pa。前処理は市松を減らすが不整合で収束せず | 記録 (診断) |
| `run_0076_cell_lmp2_eps15` / `run_0077_cell_lmp2_long` / `run_0079_cell_lmp2_converge` | **打ち手確立: cell + `lowMachPrecond=2` cfl1** (完全前処理 block-DPLUR, [plan Phase4](../../plans/accepted/time_integration-lowmach-preconditioning.md))。`run_0079` は完全収束 (200k 相当) | **NaN無・全残差 falling**。**P 市松 78→18.4 Pa (−76%)**、同一 restart A/B で rms_roUx **0.29→5.5e-4**・roOmega 1.13→8.9e-4 でプラトー打破、quasisteady STEADY。vt/l 196 (baseline 198)。図 `run_0077.../corner_oddeven_lmp2_vs_baseline.png` | active |
| `run_0078_node_lmp2_eps15` | node + `lowMachPrecond=2` cfl1 | **発散せず** (lmp1 と対照)。**P 市松 96→22 Pa (−77%)**、平均流残差 −2.5〜2.8桁。rms_roOmega は ~9.5 プラトー据え置き (node 既知 ω 問題、市松とは別) | active |
| `run_0080_cell_lmp0_double` / `run_0081_cell_lmp2_double` | **精度切り分け (2×2)**: `build-double` (全 float64, `FORGE_CUDA_BLOCKSIZE=128`) で lmp0/lmp2 を再計算 | **倍精度は市松に全く効かない**: lmp0 double=float (P odd-even 69.9=69.9)、lmp2 double=float (18.5=18.5)。**改善は 100% `lowMachPrecond=2` の前処理物理、精度は無関係**と確定 (float32 桁落ち仮説を決定的に棄却) | 記録 (診断) |
| `run_0082_cell_lmp0_dsolve` | **c' 物理 vs double Jacobian 切り分け**: lmp0 + `implicitSolvePrecision=1` (block-DPLUR Jacobian+solve を double, c' なし) | **double Jacobian 単独は市松に無効**: P odd-even **69.9 (=lmp0 float, bit 一致)**。対して lmp1 (c' + float Jac) は 21.1。→ **市松を直すのは c' 前処理の物理、lmp2 の double Jacobian は前処理系の条件付け (収束安定化) 用で主因でない** | 記録 (診断) |
| `run_0100_cell_lmp0_ref` / `run_0101_cell_lmp2_ref` / `run_0102_cell_lmp3_lhs` | **`lowMachPrecond=3` (LHS-only) 検証 A/B**: 同一 restart・同一 build で 0/2/3 を 10k step (cell, cfl1, [plan Phase4](../../plans/accepted/time_integration-lowmach-preconditioning.md))。`=3` は LHS 完全前処理のみ・RHS 散逸 `c'` を行わない新モード | **意図通り分離**: ① **step0 RHS** lmp3=lmp0 (rms_ro 8.45e-4, roe 307.36 bit一致; atomicAdd ノイズ~1e-6 のみ) vs lmp2 は 6×/5× 差=RHS 改変。② **P 市松** lmp3 rms43.5 Pa≈lmp0 38.2 (市松残存=解不変) vs lmp2 5.6 Pa (解変更で市松除去)。③ **収束** lmp0=`NOT CONVERGED(plateau, roUx 0.29 stall)` に対し lmp3=lmp2=`NOT CONVERGED(still converging)` で**プラトー打破** (roUx 0.088)。→ **LHS 前処理=収束加速 / RHS c'=市松除去** と機能分離を実証。NaN無。図 `run_0102.../corner_oddeven_lmp0_lmp2_lmp3.png` | active |
| `run_0083_cell_omegabudget` / `run_0084_node_omegabudget` | **omega 残差プラトー局在** (収束場から 200 step、`res_roOmega`/`src_jac_omega`/`transport_diag_omega` 出力) | **node: 入口×壁コーナーノード (0,1.06)/(0,2.94) で \|res_roOmega\|=568** (入口内部 111 の5倍) = rms_roOmega ~9.5 プラトーの主因。**cell: 入口静か (0.14)、ホットスポットは段差リップ (max 79)**。真因=node 入口 Dirichlet が primitive omega のみピンし保存量 roOmega 放置 (roOmega/ro vs omega 不整合 888727)・res 除外なし。図 `run_0084.../res_roOmega_map.png`。[notes/investigations/backstep-lowmach-checkerboard.md](../../notes/investigations/backstep-lowmach-checkerboard.md) | 記録 (診断) |
| `run_0085_cell_kbudget2` / `run_0086_node_kbudget2` | **k 残差局在** (同上 + `res_roK` 出力) | **保存量不整合は k にも**: node 入口 roK/ro=24–7730 vs k=4 (最大 7729 不整合)。**ただし res_roK は物理側支配**: node max 9.39 は再付着 x≈6.7 (k 生産ピーク)、入口コーナー 2.31 は副次。cell は段差リップ (18)、入口≈0。→ omega プラトー=入口 BC 支配、k=再付着物理支配 | 記録 (診断) |
| `run_0060_node3d_periodic_rans` | **3D spanwise periodic 化検証** (median-dual M4 §4.5): 2D backstep を spanwise (z) 4H×8層で hex 押し出し、side1↔side2 を Cartesian periodic、run_0059 収束場を 2D→3D 最近傍移植 (一次量 k/ω 再構成) して restart、node+implicit SST | **periodic 機構は動作**: `buildPeriodicNodeGroups` が 10997 slave merge (spanwise 面ノード数)、`conv main planes periodic 0` で除外、NS 残差は step0 で有限。**だが `wallTreatmentSST=1` で rms_roK=inf at step0**。**slip 置換(periodic 無し)でも同 inf** → periodic 非依存の **3D node SST 壁関数バグ** (wf_pk が 3D 壁で inf)。`wallTreatmentSST=0` では inf 消失・NaN なし完走 (rms_roK 4.5e19→低下) だが coarse mesh で低-Re 不適。**k/ω periodic gather 実装済**。当初 `wallTreatmentSST=1` で omega が段差凹コーナー (x=2,y=0) の wall_dist=0 ノードで inf (2D は wall_dist=0.015 で偶々救われていた) → `compute_wall_y_eff_d` に 2-ring 探索を追加し修正 (commit `da22836`)。**修正後 3D RANS が安定収束** (rms_ro 1.57e-5)。**periodic 完全検証**: spanwise で x_R=5.163–5.170 (std 0.003) = z 非依存、かつ **3D x_R 5.16 == 2D(同コード)5.18 に 0.4% 一致** (run_0062) | **検証完結** |
| `run_0062_2d_withfix` | **2D 照合**: run_0059 の 2D を wall_y_eff fix 込みコードで収束場から restart | **x_R=5.18** (旧 fix 前 7.6 から変化=コーナー omega 補正の効果、y 0.015→0.0625=第一セル高さ)。3D periodic (5.16) と一致 → periodic 検証の参照 | ref (3D==2D 照合) |
| `run_node3d_iles_unsteady` | **3D node ILES** (MUSCL+乱流モデル無し LESorRANS=0, 分子粘性 visc0.001, explicit RK3 cfl0.5)。2D 収束場を spanwise 複製+剪断層摂動 seed。80000 step | **NaN無で乱流発達** (flow time 3.32s)。剥離せん断層/再付着の解像乱れ。SLAU は安定。KEEP+WALE 比較の参照 | ref (ILES) |
| `run_inlet_profile_test` | **入口分布プロファイル機能の動作例** (`ints:{inletProfile:1}` + `inlet_profile_1.csv` 壁法則, SLAU node)。30 step | ログ `[applyInletProfiles] set 3 quantities ... 297 faces`、入口 Ux(y) が壁法則分布 (壁0/中央ピーク)。機能検証 ([plan](../../plans/accepted/boundary-inlet-profile.md)) | ref (機能例) |
| `run_node3d_keep_wale` | **node + KEEP + WALE LES** (`solver:KEEP` 復活 + `LESorRANS:1 LESmodel:1` node WALE 有効化, MUSCL explicit RK3 cfl0.3)。res_80000 発達場 seed。40000 step | **安定**: NaN/cfl0 ゼロ、残差健全低下 (rms_roUx 36→0.3)、**WALE SGS 稼働** (vis_turb max0.31/mean0.05 vs vis_lam0.001=~50×)、node 弱形式パス確認。ただし **periodic 継ぎ目 (z=0/z=4) で P が ~90 Pa 振動** (内部 ~2 Pa)。[plan](../../plans/active/convection-keep-revive-node.md) | active |
| `run_node3d_keep_slip` | **切り分け①**: run_node3d_keep_wale の spanwise を **periodic→slip** に。12000 step | **継ぎ目振動が消えて clean** (P-mean vs z が滑らかな勾配 ±10 Pa・odd-even 無し)。→ 振動は **periodic 機構固有** (一般の KEEP 低マッハ checkerboard でない) | ref (診断) |
| `run_node3d_keep_pure` | **切り分け②**: `keepDissipation:0` で **純粋 KEEP** (Roe 散逸無し), periodic。12000 step | **全域で激しい spanwise odd-even checkerboard** (~130 Pa peak-to-peak, 教科書的市松)・残差上昇 (rms_roUx 1.14)。→ **KEEP 中心流束は非散逸で低マッハ圧力 checkerboard を抑えられず、Roe 散逸が必須**と確定 | ref (診断) |
| `run_node3d_keep_1storder` / `run_node3d_keep_euler` | **切り分け③④**: KEEP+Roe periodic を **1次 (convMethod0)** / **完全Euler (visc0,thermCond0,SGS off)** で。seam欠陥場から restart | **どちらも seam 欠陥が残存** (Euler は 1次LESと数値完全一致 -105.7/60.1)。→ **2次再構成・勾配ブレンド・粘性/WALE すべて棄却**。原因は periodic マージ機構そのものと確定 | ref (診断) |
| `run_node3d_keep_mirror` | **真因特定→解消**: master(z=0)/slave(z=4) の保存量 desync を実測 (**Uz 最大 14.9 m/s 差**, ro/roe/Ux は同期; 体積比=1.0で無罪)。**NS 保存量 root→member ミラー** (`periodicMirrorNSState`, 初期化+各RK stage) を実装し seam欠陥場から restart | **seam 欠陥 150→0.2 Pa・master/slave 差→機械ゼロ**で解消。NaN無・WALE稼働。[plan §4.5.9](../../plans/active/discretization-median-dual-3d.md) | **解消確認** |
| `run_0096_node3d_keep_wale_les_fine` | **細密メッシュ LES**: 新 `mesh/backstep_les.geo` (230×50×**80**span, Δz=0.05H, 933,606 nodes) で node+KEEP+Roe(`keepDissipation:1`)+WALE。旧粗メッシュ (span8) は WALE 過散逸で定常化した反例の解像版。MUSCL/venkata, explicit RK3, `initial:backstep`, 150k step, snap@1000。メッシュ品質 AR=2.67/skew≈0 PASS | **`cfl 1.0` は発散** (rms 1000×spike@step80 で非物理プラトー, NaN手前)→**`cfl 0.5` で安定** (rms_ro/roUx/roe 単調低下)。step34406 で停止・保存 (まだ roUz≈±0.007 で 2D 的)。run_0097 の restart 元 | active (停止・保存) |
| `run_0097_node3d_keep_pure_les_fine` | **pure KEEP** (`keepDissipation:0`, Roe 散逸無し低散逸 LES) + WALE。**run_0096 の res_34000 場から同一メッシュ restart** (保存量を直接 index コピー)。cfl 0.5, RK3, 150k step | **安定・clean** (rms_ro~8e-5/roUx~0.026/roe~21 で低位定常, NaN無, state-mirror 後 spanwise checkerboard 無し)。rms_roUz が 5e-5→1.7e-4 と漸増=3D モード覚醒中。**step11684 で停止・保存** (res_11000 が 0103-0105 の seed) | active (停止・保存) |
| `run_0103_node3d_keepes_iles` | KEEP-LES スタック投入 (のつもりが **keepDiss\* を `space:` 配下に誤記=黙って無視され実質 pure KEEP σ=0** + advGauge node + ILES)。run_0097 res_11000 seed、10186 step で中断 | **貴重な σ=0 対照になった**: spanwise P 市松は 0.05-0.4 Pa と無傷だが、**再循環域 (x2.2-5.5) に x 方向柱状 odd-even が seed 12.5→47.3 Pa mean (max 100) で飽和** (有界・発散せず)。3D 遷移進行 (roUz 77×)。図 `pfield_slice_compare.png`。VERDICT=NOT CONVERGED (非定常LESの想定) | 記録 (σ=0 対照) |
| `run_0106_es_probe_s005_cp1` / `run_0107_es_probe_s005_cp0` / `run_0108_es_probe_s002_cp0` | σ/c' プローブ (のつもりが同じ config 誤りで **3 本とも pure KEEP**) 4000 step | **run_0103 res_4000 と全数値一致** (x-oddeven 47.27/shear Uy 7.102) = node 決定論性の副次実証 + config 誤りの発覚経緯。keepDiss\* は **top-level キー** ([memory: keepdiss-keys-toplevel-trap]) | 記録 (誤記・決定論確認) |
| `run_0109_es_probe_s002_cp1` | **ES 散逸プローブ (正)** σ0.02+c' jump2, 4000 step, seed=run_0097 res_11000 | 再循環 x-市松 mean **47.3→28.6 Pa (−40%)**、せん断層 Uy rms 7.12 (不変=渦無傷)、roUz 4.7e-4 (遷移 1.8× 減速) | ref (σ/c' カーブ) |
| `run_0110_es_probe_s005_cp1` | 同 σ0.05+c' | x-市松 **22.1 Pa (−53%)** max 82、Uy rms 6.89 (−3%)、roUz 3.5e-4 | ref (σ/c' カーブ) |
| `run_0111_es_probe_s005_cp0` | 同 σ0.05+**full c** (keepDissCprime:0) | x-市松 **11.0 Pa (−77%)** max 37 と最小。ただし **Uy rms 6.20 (−13%)・3D 遷移 4× 減速** = 音響減衰はせん断層 2-4Δ 構造も食う。市松レバーとしては σ↑ (全モード一律) より c' の方が効率的 | ref (σ/c' カーブ) |
| `run_0112_les40k_iles` | 長尺 SGS 比較 (σ0.05+c') — **step1052 で中断** (前処理散逸の開発を優先、σ設定も再検討対象に) | VERDICT 記録済。`run_0113/0114` (WALE/σ) は未実行のままディレクトリ破棄 | 破棄予定 (中断記録) |
| `run_0115_es_probe_s010_cp1` / `run_0116_es_probe_s020_cp1` | σ0.10 / σ0.20 (+c') プローブ | x-市松 15.6 / 11.8 Pa、Uy rms 6.46 / 6.25 = σ 増は市松・物理を**比例して**削る | ref (σ/c' カーブ) |
| `run_0117_es_probe_s010_cp0` / `run_0118_es_probe_s002_cp0jst` | σ0.10+full c / σ0.02+full c (JST k₄ 相当の文献アナログ点) | 7.4 Pa (Uy 5.98) / 17.3 Pa (Uy 6.57)。**σ・c' はどちらを動かしても同一トレードオフ曲線上** (σ0.02fullc≈σ0.10c', σ0.05fullc≈σ0.20c') = 一様固有値スケーリングでは市松と 2-4Δ 物理を分離不能 → 前処理散逸 plan 起案の根拠 | ref (曲線縮退の実証) |
| `run_0119_es_regr_s005_cp1_newbin` | 回帰: `keepDissPrecond` 実装後バイナリで run_0110 再実行 (precond=0) | x-oddeven 22.14/54.14・Uy 6.885 = **run_0110 と一致 (既定経路無影響)** | ref (回帰) |
| `run_0120_es_probe_precond_s005` / `run_0121_es_probe_precond_s002` | **Turkel 前処理散逸** (`keepDissPrecond:1`, [plan](../../plans/active/convection-keep-diss-lowmach-precond.md)) σ0.05 / σ0.02 | 安定・NaN無。d2 メトリクス 22.3 / 27.5 Pa = c' 版と同水準**だが物理は改善** (Uy 6.94/7.16 ≥ c' 版 6.89/7.12)。case/35 L1 で真の市松は 142× 減衰 → **backstep の ~22 Pa 床は数値市松でなくリップ 2-4Δx の解像限界物理と確定** (根治はメッシュ細分化) | ref (precond 検証・床の帰属確定) |
| `run_0122_les40k_iles_precond` | **SGS 比較①: ILES** (precond σ0.05 jump2 + advGauge, seed=run_0097 res_11000, cfl1.0・出力1000毎)。step29847 でユーザー判断により区切り | **cfl1.0 安定** (run_0096 の cfl1 発散は旧散逸+一様IC 起因と切り分け)・NaN無。**rms_roUz 4e-4→0.163 (roUx/roUy 並み) = 3D 遷移完了、本物の 3D 乱流場に到達**。res_29000 が 0123/0124 の共通 seed。VERDICT=NOT CONVERGED rising (非定常LESの想定) | active (SGS 比較 seed) |
| `run_0123_les40k_wale_precond` | 同**②: WALE** (`LESmodel:1`)。**seed=run_0122 res_29000 (発達 3D 乱流場)**、30k 完走 | **安定・NaN無・vis_turb 稼働** (mean 0.0036≈3.6×分子, nonzero 96%)。せん断層 Uy/Uz rms 6.4/6.8 = 乱流維持。VERDICT=NOT CONVERGED plateau (非定常想定) | active (SGS 比較 WALE) |
| `run_0124_les40k_sigma_precond` | 同**③: σ-model** (`LESmodel:2`)。同 seed。step10636 でユーザー判断により終了 | 安定・NaN無・vis_turb mean **0.0060 = WALE の 1.7 倍散逸的** (TGV での σ>WALE 散逸傾向と整合)。res_10000 まで | active (SGS 比較 σ, 短縮) |

| `run_0125_ddes_keep_es` | **DDES × 単一スキーム KEEP+ES** ([iddes plan §4.8 設計更新](../../plans/active/turbulence-iddes-sst.md)): run_0035 と同一 (3D 857k, DESmode:1, implicit) で flux のみ KEEP+matrix ES (σ0.05 jump2 precond)、500 step 機能確認 | **NaN 無し・f_d zoning 正常** (せん断層 0.974 = SLAU 0.976、泡内 0.000)。付着 BL 帯 f_d 0.43 vs SLAU 0.28 (遮蔽保持、要観察)。**Phase 1.5 前提検証②合格** | ref (KEEP+ES DDES 機能確認) |

**SGS パイロットの総括 (2026-07-19)**: 目的 (確定版 KEEP-LES スタックの実地統合試験) は達成 —
node×KEEP×precond ES 散逸×advGauge×WALE/σ/ILES が cfl1.0 で全て安定・3D 乱流を維持・SGS 活性正常。
**定量的な再付着位置比較 (Driver–Seegmiller x_R=6.26) は本構成では狙わない (ユーザー判断)**:
①入口が一様 Pt/Tt で流入境界層の厚み・乱流度が合っていない (要 `inletProfile` + 合成乱流/リサイクリング)、
②near-wall 未解像 (y+~数十) で壁モデル無し (要 wall-modeled LES または IDDES = [turbulence-iddes-sst plan](../../plans/active/turbulence-iddes-sst.md))。
残存する x 方向縞 (~22 Pa) はリップ Δx=0.13 の解像限界物理 (run_0109-0121 で帰属確定)。

> **P・ρ 振動の正体と打ち手 (run_0071–0079, 2026-06-28)**: node/cell 共通の P・ρ「振動」は **低マッハ再循環域に
> 局在する odd-even チェッカーボード** (collocated 圧力-速度デカップリング) と確定。段差角を1セル過ぎた x≈2.13 の
> y 分布で P が節点ごと **±50 Pa 鋸歯** (2次差分 rms=全振幅の 48–54%)、Ux は滑らか。時間リミットサイクル (<10 ppm)・
> リミッタチャタリング・**float32 桁落ち (振幅 5e-4≫ε、かつ float32 のままスキーム変更で消える) は棄却**。機構は
> SLAU 質量流束の圧力散逸が `c_hat` スケールで M→0 縮退 ([plan §9](../../plans/accepted/time_integration-lowmach-preconditioning.md))。
> **打ち手 = `lowMachPrecond=2`** (完全前処理 block-DPLUR): P 市松 **−76〜77%** (cell 78→18, node 96→22 Pa)、cell は同一 restart A/B で
> **残差プラトー打破** (rms_roUx 0.29→5.5e-4・roOmega 1.13→8.9e-4)。**`lowMachPrecond=1` は LHS 非前処理で不整合→発散/停滞、使わない**。
> 詳細: [`notes/investigations/backstep-lowmach-checkerboard.md`](../../notes/investigations/backstep-lowmach-checkerboard.md)。

> **メッシュ/面ベクトル検証 (2026-06-21, user 指摘)**: 双対メッシュは正常 — 体積保存 88=88・closure 6.7e-6・負体積0・
> 体積比4.4・発散種域は均一セル。**面ベクトルの向きも正常** (内部双対面 反転0%/cos min0.936、境界半割面 内向き0/556、
> 段差コーナーも外向き整合)。→ 発散は**メッシュでも面ベクトルでもなくスキーム**。種域 Mach≈0.005-0.03 の低マッハ淀みで
> SLAU 風上散逸が消え **collocated median-dual の odd-even チェッカーボード**が減衰しない。CFL を下げると逆に速く成長
> (低マッハ stiffness)、restart 初期残差がノズル(超音速)7.95e-5 に対し backstep 0.055=700倍。`nodeWallDirichlet`/
> `nodeMidpointFx`/`gradLSQ`/CFL いずれも無効。`lowMachPrecond=2` (run_0042) も **step0 rms_ro=77.7→step1 NaN で逆効果**
> (node で前処理破綻、case/29 と同傾向)。**config レベルのレバーは全滅** → 低マッハ node 用の散逸/圧力-速度結合という
> スキーム改修が必要 (plan discretization-lsq-gradient / median-dual の低マッハ未対応領域)。ノズル(超音速)が回るのは
> 風上散逸が効くためで、低マッハ剥離は node-centered density-based の最難ケース。

> **SST-DES (DDES) T1-B 機能検証**: f_d が剥離域で立ち付着 BL で 0、NaN なし、を満たせば Phase 1
> 合格。**収束 VERDICT は NOT CONVERGED が想定どおり** (DDES は本質的に非定常で残差は下げ止まらない)。
> 解像乱流の定量化 (RMS・スペクトル −5/3・x_R/H) は Phase 1.5 (DES 用低散逸 flux) 後の仕事。
> 計画: [turbulence-iddes-sst.md](../../plans/active/turbulence-iddes-sst.md)。

> **node-centered 検証メモ (run_0036–0039, 2026-06-21)**: backstep を node モードで回す試み。
> ① 既存 `backstep_2d.msh` は厚さ0.1の押し出しスラブ=実体3Dで median-dual が未対応 (3D dual=M4 未実装) のため、
> 平面2D `backstep_planar.geo` を新規作成。② node モードは (当時) **RANS/層流/フラグ ON-OFF いずれも step18 で発散**。
> 種は **段差直後 再循環域の隣接 CV odd-even チェッカーボード** (median-dual 中点双対の既知弱点)。
> 当時の結論「node モードは backstep では未成立」は **下記 2026-06-27 の結果で覆った** (発散の真因は積分法の交絡)。
> 関連 plan: [discretization-median-dual.md](../../plans/active/discretization-median-dual.md) (M4),
> [discretization-lsq-gradient.md](../../plans/active/discretization-lsq-gradient.md) (checkerboard),
> [discretization-node-boundary-ghostless.md](../../plans/active/discretization-node-boundary-ghostless.md)。

> **node モードで backstep が発散しなくなった (run_0049/0050, 2026-06-27)**: run_0036–0042 の step18 発散は **node 固有でも
> median-dual checkerboard 固有でもなく、explicit RK3 という積分法の交絡**だったと確定。当時の node 実験は全て RK3
> explicit で、同 planar mesh の cell-RK3 も同様に発散していた (run_0043/0045)。**積分法だけを block-DPLUR implicit
> (`timeIntegration:11`) に替えると、他設定 (node, gradLSQ=1, nodeWallViscGradFlux=1, 1次, SST, cfl_pseudo=1) は
> run_0036 と同一のまま step18 を越えて完走** (run_0049=PoC 5000step / run_0050=30000step)。
> block-DPLUR の対角優位が低マッハ淀み域の odd-even を減衰させたため (cell-implicit run_0001 が cell-explicit より
> 安定だったのと同理)。**ただし収束はしていない**: run_0050 の `check_convergence` VERDICT は
> **`NOT CONVERGED (stalled/plateau — needs scheme change)`**。rms_ro 3.9桁・運動量 3.3桁・roe 4.2桁・roK 2.2桁は
> 落ちるが **rms_roOmega が 1桁プラトー** (壁近傍 ω ソースの SST 既知挙動、steps 増では解消しない)。場自体は物理的
> (NaN無・T~288K・再循環泡 Ux<0 再現)、machmax/pmax は `check_quasisteady` で `STEADY`。
> **残課題**: ① rms_roOmega プラトー脱却 (SST 共通)、② 2次/MUSCL での node 収束 (run_0049/0050 は
> 1次)、③ スパン periodic + 3D median-dual (M4) は依然未実装で 3D node DDES は別ゲート。

> **node/cell/SU2 三者クロスチェック (run_0057/0058/0059, 2026-06-28)**: 同一 planar mesh + 同一 BC + 壁関数 (autowall)
> で再付着長 x_R/H を三者比較 (x_R は**第一内層 Ux 符号反転**で統一; SU2 は Cf 符号反転とも一致)。
> **SU2=7.89 / forge cell=7.95 / forge node(`nodeWallDirichlet=1`)=6.71 / 実験(Driver-Seegmiller)=6.26**。
> 結論: ① **forge cell ≈ SU2 (差 0.8%) → forge-SST 較正は正しい**。② **SU2 でも 7.9 で実験 6.26 と不一致 → 差は
> forge バグでなく粗い y+~62 mesh+壁関数の限界** (実験値は低Re 壁解像 mesh が要)。③ node は cell/SU2 より ~1H 短い
> (6.71) が実験には最も近い; node-cell 差 ~1H は実在 (要追究)。
> **壁速度の落とし穴 (user 指摘で判明)**: node は **`nodeWallDirichlet=1` を明示しないと壁速度がゼロにならない**
> (既定0 は半割面弱形式のみで壁ノード DOF がドリフト、再循環域で負速度)。run_0056(Dir=0) の「x_R=9.59 壁バンド法」は
> この負ドリフトを拾った**測定アーティファクト**で、第一内層法では 6.71。壁関数 (`ransWallFunction_d.cu`) は速度を読むだけで
> Ux を書かないので壁すべりは作らない。SU2 cfg は [`run_0058_su2_sst/case.cfg`](run_0058_su2_sst/case.cfg)。

## 3D node ILES (非定常・MUSCL のみ・乱流モデルなし)

`run_node3d_iles_unsteady`: 3D spanwise periodic backstep を **node + MUSCL(convMethod2)+ 乱流モデルなし(LESorRANS=0)+ 分子粘性(visc0.001, Re_H~40000)+ 非定常 explicit RK3** で計算 (ILES: MUSCL の数値散逸が SGS 代替)。2D 収束場を spanwise 複製 + 剪断層に Uz±3 摂動で seed。
- **node 低マッハ非定常 explicit が安定動作** (従来「node 低マッハ explicit は発散」だったが、ILES 設定=MUSCL 散逸+発達場 restart で安定)。場健全 (ro 1.18..1.20, P/T 物理的, NaN なし)。
- **解像乱流が発達**: 摂動は最初 2000step (t0.08s=flow-through 14%) では減衰したが、**長く回す (40000step, t1.74s≈2.9 flow-through) と剪断層 KH 不安定から 3D 遷移→乱流再付着が発達** (rms_roUz 0.05→0.53 持続, |Uz|平均1.66, spanwise 渦構造・乱流エディが明瞭)。図 `iles_turbulent.png` (Ux/Uz/変動強度)。
- **教訓**: LES は数 flow-through 回さないと発達しない (14% で止めると未発達=減衰に見える)。

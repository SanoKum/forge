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
| `run_0083_cell_omegabudget` / `run_0084_node_omegabudget` | **omega 残差プラトー局在** (収束場から 200 step、`res_roOmega`/`src_jac_omega`/`transport_diag_omega` 出力) | **node: 入口×壁コーナーノード (0,1.06)/(0,2.94) で \|res_roOmega\|=568** (入口内部 111 の5倍) = rms_roOmega ~9.5 プラトーの主因。**cell: 入口静か (0.14)、ホットスポットは段差リップ (max 79)**。真因=node 入口 Dirichlet が primitive omega のみピンし保存量 roOmega 放置 (roOmega/ro vs omega 不整合 888727)・res 除外なし。図 `run_0084.../res_roOmega_map.png`。[notes/investigations/backstep-lowmach-checkerboard.md](../../notes/investigations/backstep-lowmach-checkerboard.md) | 記録 (診断) |
| `run_0085_cell_kbudget2` / `run_0086_node_kbudget2` | **k 残差局在** (同上 + `res_roK` 出力) | **保存量不整合は k にも**: node 入口 roK/ro=24–7730 vs k=4 (最大 7729 不整合)。**ただし res_roK は物理側支配**: node max 9.39 は再付着 x≈6.7 (k 生産ピーク)、入口コーナー 2.31 は副次。cell は段差リップ (18)、入口≈0。→ omega プラトー=入口 BC 支配、k=再付着物理支配 | 記録 (診断) |

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

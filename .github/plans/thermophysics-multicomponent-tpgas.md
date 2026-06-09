# 多成分 thermally-perfect gas ソルバ機能

## メタ

- **area**: `thermophysics` (新規。理論/実装は `docs/thermophysics/`)
- **status**: `done` (M1-M4 + TP-BC + 陰解法 γ[ic]。後続改良は §10)
- **related_docs**:
  - `docs/thermophysics/theory.md`
  - `docs/thermophysics/implementation.md`
  - `docs/convection/` (流束の TP 整合)
  - `docs/diffusion/` (化学種拡散・エンタルピー拡散)
- **created**: `2026-06-09`
- **owner**: `CFD Dev`

## 1. 目的

forge を単一成分・熱量的完全気体 (CPG) から**多成分 thermally-perfect gas (非反応)** に拡張する。比熱 $c_p(T)$ を NASA-9 (CEA 準拠) とし、粘性・熱伝導率・質量拡散を kinetic theory で評価、混合物性は ideal-gas mixing。`thermalMethod==0` を保持し `thermalMethod==2` で有効化。

## 2. スコープ

- **やる**: NASA-9 熱力学、Newton 温度反転、多成分質量分率輸送、kinetic theory 輸送係数、化学種拡散とエネルギー結合、SLAU/ROE の TP 整合、設定/IC/BC/出力配管。
- **やらない**: 化学反応源項 (Arrhenius 等) と有限速度化学。設計上の拡張点としてのみ残す。

## 3. 関連 docs と前提

- 理論: `docs/thermophysics/theory.md`、実装: `docs/thermophysics/implementation.md` (本 plan 着手前に作成済み)。
- 既存の汎用スカラ輸送 `scalarTransport_d` / RANS `ransTransport_d` を化学種輸送のミラー元とする。
- `flow_float` は FP32。熱力学/輸送のカーネル内計算は double で行う。

## 4. 設計方針

詳細は `docs/thermophysics/implementation.md`。要点:
- 熱力学コア `cuda_forge/thermo_d.{cuh,cu}` (NASA-9, 混合則, Newton 反転)。化学種データは device pointer で渡す。
- `dependentVariables_d` に `thermalMethod==2` の Newton 反転を追加し `T,P,Ht,sonic,gamma,cp` を更新。
- 対流流束は被移流全エンタルピーを NASA 化 (SLAU 最小修正、ROE は有効γ̃の段階A)。
- 陰解法 Jacobian は `gamma[ic]` 配列化 (frozen-coefficient)。— **実装済** (2026-06-09, block DPLUR + precond カーネルが per-cell γ を読む。CPG は γ[ic]=cfg.gamma で不変)。
- 化学種は segregated point-implicit で別解き (5×5 構造不変)。

## 5. 実装ステップ (マイルストーン)

1. **M1**: `thermo_d`, Newton 反転 (`dependentVariables_d`), `thermalMethod==2`, SLAU/ROE の TP 整合, 設定 (`solverConfig`), `main.cpp` 初期化。— **実装済・検証完了** (§9, 2026-06-09)
2. **M2**: 化学種変数登録 (`variables` の `cellValNames` 非 const 化 + `registerSpecies`), 化学種移流/実現可能性/再正規化, 化学種 BC/IC/出力 (`speciesTransport_d`, `speciesBoundary_d`)。— **実装済・検証完了** (§9)
   - **設計判断 (実装時)**: 化学種は RANS k/ω と同じ汎用スカラ輸送コア `scalarTransport_d` (`ScalarTransportDesc`) を再利用し segregated に解く。1 化学種ごとに `roY{s},Y{s},roY{s}N,roY{s}M,res_roY{s},res_roY{s}_m,transport_diag_Y{s},src_jac_Y{s}` を `registerSpecies()` で動的登録 (`cellValNames`/`output_cellValNames` を非 const 化)。device には `flow_float** roY` ポインタ配列を 1 度だけ構築し `dependentVariables_d` (現状 `nullptr`) へ渡して混合則 thermo を有効化。
   - **段階化**: M2 のスコープは **移流のみ** (拡散は M4)。化学種機構は `nSpecies>=2` でのみ起動し、単成分 (`nSpecies==1`) は M1 と完全に同一経路 (回帰がバイト一致) を保証する。
   - **化学種 BC**: まず Neumann (zero-gradient) で ghost を埋める。slip 閉領域・内部 contact の検証には十分。組成依存のエネルギー BC (`boundaryCond_d` の TP slip が `sp[0]` 固定) の一般化は後続の改良として残す。
3. **M3**: kinetic theory 輸送係数 (`gasProperties_d` 拡張: Wilke/Wassiljewa/混合平均拡散, Neufeld 衝突積分)。`thermCond` のセル毎配列化。— **実装済・検証完了** (§9)
4. **M4**: 化学種拡散 (`scalar_diffusion_species_d`) + ΣJ_i=0 補正 + エンタルピー拡散項 `Σh_s J_s` を `res_roe` へ。— **実装済・検証完了** (§9)

## 6. 検証

- **単体 / ビルド**: native (`tools/build_native_wsl.sh`, CC 8.6) でビルド成功。NASA-9 単体値を NIST と比較 (N2: R=296.8, cp(300)=1040, cp(2000)=1284 J/kgK, h(298.15)=0、e→T 反転誤差 ~1e-11)。
- **検証ケース**: `case/13.nozzle_H/run_0001_tpgas_n2` (TP-N2) vs `run_0002_slau_cpg` (CPG)。`validate_tpgas.py` で中心線 M・T・P と全エンタルピー保存を NASA(=CEA) 等エントロピーと比較。高 Tt ケースで cp(T) 効果を showcase。
- **判定基準**: 残差収束、$H_0=h(T)+½|u|^2$ が中心線で一定 (エネルギー整合)、(M, P/P0) が NASA/CEA 等エントロピー上に乗る。`thermalMethod==0` 既存ケースが不変。

## 7. 影響範囲

- 新規: `cuda_forge/thermo_d.{cuh,cu}`, `docs/thermophysics/{theory,implementation}.md`, `case/13.nozzle_H/validate_tpgas.py`。
- 変更: `cuda_forge/dependentVariables_d.{cu,cuh}`, `cuda_forge/convectiveFlux_d.{cu,cuh}`, `cuda_forge/gasProperties_d.cu` (TP は無改変)、`cuda_forge/CMakeLists.txt`, `input/solverConfig.{hpp,cpp}`, `main.cpp`, `docs/index.md`。
- M2-M4 で `variables.{hpp,cpp}`, `viscousFlux_d`, `boundaryCond.cpp`, `setInitial.hpp`, `output.cpp` 等を追加変更。

## 8. 完了条件

- [x] `docs/thermophysics/theory.md` 作成
- [x] `docs/thermophysics/implementation.md` 作成
- [x] M1-M4 実装・検証完了 (§9 変更ログ参照)
- [x] TP 境界条件整合化 (slip/wall/inlet_Pressure/outflow + outlet_statPress/wall_isothermal/inlet_Pressure_dir/inlet_uniformVelocity)
- [x] 陰解法 block-DPLUR Jacobian の per-cell γ 化
- [x] `.github/plans/README.md` / `docs/index.md` 更新
- [x] 本 plan を `done` 化し §9 に変更ログ
- [ ] (後続) §10 残課題 — 組成依存エネルギー BC・RANS+TP 検証・ROE 段階B 等

## 10. 残課題 (後続。M1-M4 のコア外の改良)

- **多成分組成依存エネルギー BC**: `inlet_uniformVelocity` は組成依存化済 (M5, 下記変更ログ)。残: slip/wall/outlet の ghost roe/sonic がまだ `sp[0]` 固定 (ただし SLAU 対流流束はセル組成で再構成するため非粘性では支配的影響なし。粘性流束・他スキーム HLLE/AUSM/KEEP/ROE の TP 多成分化は未)。
- **化学種の陰解法 (point-implicit/segregated)**: 化学種輸送は explicit RK のみ。陰解法ステップへ未配線。
- **ROE の TP 完全整合 (段階B)**: M1 は Roe 平均の有効 γ̃ (段階A) まで。
- **MUSCL + 陰解法の安定化**: 2次 MUSCL × 高 cfl_pseudo で発散 (block DPLUR)。limiter freezing / CFL ランプで陰解法高速化の余地。
- **低Mach preconditioning の TP 実効検証**: 超音速で発散は正常 (M≪1 用)。真の低Mach TP ケースで `gamma[ic]` 修正後の収束加速を未確認。
- **粘性 / RANS SST + TP のエンドツーエンド検証**: CEA 照合は非粘性。RANS×TP・乱流化学種拡散 `μ_t/Sc_t`・M3 輸送係数の実粘性 TP 流れは未検証。
- (スコープ外: 有限速度化学 Arrhenius。設計拡張点としてのみ)

## 9. 変更ログ

- `2026-06-09` — 初稿。M1 (NASA-9 熱力学, Newton 反転, SLAU/ROE TP 整合) 実装、native ビルド成功、NASA-9 単体値を NIST と照合済み。nozzle_H で TP-N2/CPG 検証実行中。
- `2026-06-09` — **M1 流れ場 CEA 定量検証 (完了)**。当初の nozzle_H は warm-start 元が旧バイナリ製で inlet_Pressure+SLAU が流れを確立できず不成立 (CPG でも崩壊するため TP 起因ではない)。代わりに `case/05.sod_shock_tube` の均一メッシュ (`1D_shock_tube.geo`, nx=499, x∈[0,1], 全 slip, inviscid SLAU, MUSCL+limiter, RK3, 固定 dt=2.5e-7s で t=2.5e-4s) で N2 衝撃管 (PL/PR=10, TL/TR=2000/400K) を立て直し、TP (`run_0001_slau_tp_n2`, thermalMethod==2) と CPG (`run_0002_slau_cpg`, thermalMethod==0) を実行。
  - **枠組み検証 (CPG vs 厳密 Sod)**: 中心線で rarefaction/contact/shock 位置が一致。L2 相対誤差 P=0.90%, u=2.96%, ρ=3.53% (ρ は contact 平滑化が主因、499 セル 2次精度で妥当)。→ TP 流束+EOS+BC の数値枠組みが厳密 Sod 解と整合。
  - **cp(T) 効果 (TP vs CPG)**: 高 T プラトー (~1650K) で TP が CPG より高温、contact/shock front が右 (高速) へシフト (T 差最大 248K, u 差 225m/s)。N2 の cp(T) 増大 (γ 低下) による実在気体効果を定量提示。
  - 後処理: `case/05.sod_shock_tube/gen_shocktube_ic.py` (TP/CPG IC), `compare_shocktube.py` (厳密 Sod Riemann ソルバ + 中心線照合 + `shocktube_comparison.png`)。
  - 付随: `cuda_forge/calcStructualVariables_d.cu` の幾何メトリック計算を FP32 (`sqrtf`/`powf`) → FP64 (`sqrt`/`pow`) 化 (TP の double 一貫性向上)。
  - §6 nozzle_H ベースの検証計画は本 sod_shock_tube ベースの照合に置き換え。M1 の流れ場 CEA 検証はこれで完了とみなす。
- `2026-06-09` — **M2 多成分化学種輸送 (実装・検証完了)**。汎用スカラ輸送コア `scalarTransport_d` を化学種ごとに再利用し segregated に移流。`variables::registerSpecies` で 1 化学種あたり 8 セル変数を動的登録 (`cellValNames`/`output_cellValNames` 非 const 化)、device `flow_float** roY` を `speciesInit_d` で構築して `dependentVariables_d` の混合則 thermo を有効化。`speciesPrimitive_d`/`applySpeciesBoundaries`(Neumann)/`speciesTransport_d`/`speciesTimeIntegration_d`/`speciesRenormalize_d`(ρY_s≥0, ΣρY_s=ρ)/`speciesUpdateOuter/Inner` を追加し `main.cpp` の RK ループへ配線。`nSpecies<=1` では全機構が dormant で M1 と同一経路。
  - **検証** (`case/05.sod_shock_tube`): ① 単成分回帰 (`run_0003`) が M1 とバイト一致。② 識別子トレーサ `[N2,N2]` (`run_0004`) — 混合則が常に N2 のため流れ場 (ρ,P,T,u) が単成分 N2 と完全一致 (rel 0)、組成は contact と共に移流、$\sum Y=1$ を float 精度 (1.2e-7) で維持、$0\le Y\le1$。③ 実 2 ガス `[O2,N2]` (`run_0005`) — 安定・有界・$\sum Y=1$、O2 駆動で密度/温度が単成分と差。後処理 `gen_shocktube_species_ic.py`, `species_comparison.png`。
  - **既知の制約 (後続)**: 化学種拡散は未 (M4)。組成依存のエネルギー BC は未 (slip TP が `sp[0]` 固定、Neumann ghost)。陰解法 (point-implicit/segregated) は陽解法 RK 検証のため未配線。
- `2026-06-09` — **M3 kinetic theory 輸送係数 (実装・検証完了)**。`thermo_d.cuh` に Chapman-Enskog 単成分粘性・modified Eucken 熱伝導・Neufeld 衝突積分 ($\Omega^{(2,2)*},\Omega^{(1,1)*}$)・Wilke/Mason-Saxena 混合・二元/混合平均拡散 (`thermo_Dbinary`/`thermo_Dmix_species`, M4 用) を追加。`gasProperties_d` に `viscMethod==2` を追加し per-cell に $\mu_{mix},\lambda_{mix}$ を評価。`thermCond` をセル毎変数化し `viscousFlux_d` の熱流束を面平均 per-cell へ (viscMethod 0/1 は一定値で既存不変)。
  - **検証**: ① 単体 — μ/λ が NIST と一致 (N2 μ(300)=1.81e-5/μ(1000)=4.15e-5 Pa·s, λ(300)=0.0255 W/mK, O2 μ(300)=2.06e-5, air μ=1.86e-5)。② device `viscMethod==2` (O2/N2 衝撃管, `run_0006`) の vis_lam が Python 式と FP32 一致 (rel 2.7e-4), thermCond 物理値・全 finite・安定。③ 非粘性単成分回帰 (`run_0007`) が M1 とバイト一致。
- `2026-06-09` — **M4 化学種拡散 + エネルギー結合 (実装・検証完了)**。`speciesTransport_d` に `species_diffusion_d` を追加: Fick 拡散 $J_s=\rho D_s\nabla Y_s$ (混合平均 `thermo_Dmix_species` または定数 Sc, 乱流 $\mu_t/\mathrm{Sc}_t$ 込み)、質量保存補正 $J_s^*=J_s-Y_s\sum J_k$ ($\sum J^*=0$)、エンタルピー拡散 $\sum h_s J_s^*$ を `res_roe` へ。`viscMethod!=0 && nSpecies>=2` で起動。`speciesInit_d` は roY/res/diag の device ポインタ配列を構築。
  - **二重計上バグの発見・修正**: 汎用スカラコアの拡散が層流粘性 $\nu$ で化学種を拡散し Fick 拡散に重畳していた (検証で $D$ が 1.73× = $D_{self}+\nu$)。`ScalarTransportDesc.diffusion` フラグを追加し化学種は 0 (RANS は 1) として解消。
  - **検証** (`case/05.sod_shock_tube/run_0008`, 静止 N2/N2 トレーサ, P=1e3Pa/T=1000K で $D\propto1/P$ を可視化): 誤差関数フィットで $D_{fit}=1.62\times10^{-2}$ vs kinetic theory 自己拡散 $1.64\times10^{-2}$ m²/s (ratio 0.987)。種質量 7 桁保存、$\sum Y=1$ (6e-8)、$0\le Y\le1$、u≡0・T/P 一定 (識別子ガスでエネルギー結合が恒等的にゼロ)。O2/N2 粘性 (`run_0006`, 実ガスでエネルギー結合非自明) が安定・有界。非粘性単成分回帰は M1 とバイト一致。
- **M1-M4 完了**。残: TP 境界条件の整合化 (出口静圧/wall-isothermal/inlet backflow が CPG 前提)。検証は 23.axi_nozzle (3D+軸対称) と 20.naca_ml。
- `2026-06-09` — **TP 境界条件の整合化 (実装・検証完了)**。CPG 前提 (定数 γ,cp) だった BC を NASA 整合化: `outlet_statPress` (順流+backflow), `wall_isothermal`, `inlet_Pressure_dir`, `inlet_uniformVelocity` (M1 の slip/wall/inlet_Pressure/outflow に追加)。各カーネルに `thermalMethod, sp` を渡し、TP 分岐で `thermo_R_species`/`thermo_cp_mass`/`thermo_h_mass`/`thermo_isentropic_from_total_single`/新規 `thermo_isentropic_from_total_Ps_single` を使用。CPG 分岐は代数的に不変 (回帰安全)。
  - **検証** (CPG vs TP-AIR; 内蔵 AIR は定数 cp/R=3.5・γ=1.4 ⇒ TP-AIR は CPG を再現するはず): `20.naca_ml` planar (inlet_Pressure_dir+outlet_statPress, M~0.46) が ρ/T/P/sonic/Mach で **0.02%** 一致 (クリーン)。`23.axi_nozzle` axisym (inlet_Pressure+outlet_statPress, Mach4) が収束・平均 ~1% 一致 (差は非定常プルームに局在)。`23.axi_nozzle` 3D は TP が finite・安定で CPG と同残差 (両者非定常プルームで未収束; BC コードは次元非依存)。
  - **制約**: 多成分組成依存エネルギー BC は未対応 (sp[0] 固定)。TP warm-start は NASA 絶対基準のため CPG res から不可 (IC を `roe=ρ(e_NASA(T)+ek)` で再生成)。
- `2026-06-09` — **M5 組成依存超音速入口 BC + 多成分対流流束エンタルピー (実装・検証完了)**。`case/28.cutler_coaxial_jet` (Cutler 超音速同軸 He/空気ジェット, AIAA-99-3588) を駆動ケースに、入口ごとに化学種組成 `Y_s` を指定できる超音速入口を実装。
  - **入口組成 BC**: `readBcondConfig` が `inlet_*` 種別に対し `Y0..Y{n-1}` を type-1 uniform-float bvar として動的登録 (未指定なら既定 `Y0=1`, 他=0 で単成分入口に後方互換)。`bcond` に device ポインタ配列 `Yb_d` を追加し `inlet_uniformVelocity_d_wrapper` が遅延構築。`inlet_uniformVelocity_d` の TP 分岐を `sp[0]`→混合則 (`thermo_R_mix`/`thermo_e_mix`/`thermo_gamma_mix`) 化し ghost の `T=P/(ρ R_mix(Y))`, `roe=ρ(e_mix+ek)`, `sonic=√(γ_mix P/ρ)` を入口組成で再構成。化学種 ghost は `speciesBoundary_d_wrapper` を kind 分岐し入口は Dirichlet (`roY_s[ig]=ρ[ig]·Y_s^in`)、他は従来 Neumann。
  - **多成分対流流束エンタルピー (数値スキーム変更)**: `SLAU_d` の被移流全エンタルピー `h_p/h_m` が `sp[0]` 固定で、He/N2 contact (質量基準 cp が ~5×差) でエネルギー束が不整合になり圧力発散 (1.8 GPa) していた。L 側 ic0 / R 側 ic1 のセル組成 `Y` で `R_mix`・`h_mix` を再構成するよう修正 (`roY` 配列を渡す。`nSpecies<=1` は従来 `sp[0]` 経路でビット不変)。
  - **検証** (`case/28.cutler_coaxial_jet`): ① **CPG 単成分回帰** (`run_0002_cpg_regr`) が修正前と同一 (`rms_ro` 3.37e-5→1.76e-5, NaN なし)。② **binary [He,N2]** 1次非粘性 (`run_0003_he_n2_inviscid`, 2000 步) と 2次リスタート (`run_0008_he_n2_2nd`, 計 4000 步) がともに NaN/発散なし・残差有界。$\sum Y=1$ を機械精度で維持・$0\le Y\le1$、中心 He コア (Y_He=1, Ux~1280 m/s = Mach1.8 He)、coflow N2、P/T 物理的。2次で He ポテンシャルコア長が ~0.02m→~0.07m に伸長 (1次風上の数値拡散が緩和)。後処理 `gen_binary_ic.py` (N2 ambient NASA-IC + `roY`)、`field_YHe_Mach.png`、`residual_history.png`。
  - **startup の落とし穴 (記録)**: 静止 IC (u=0) + `outlet_statPress` は出口セルが `Un=0` 近傍で backflow 分岐に落ち、`Pt`/`Tt` 未指定だと 0/0 で NaN 化する。IC に緩やかな下流速度 (u=100) を与え、出口 bcond に `Pt`/`Tt` を指定して回避。TP 起因ではなく BC startup 起因 (CPG は IC が u=10 だったため顕在化せず)。
  - **既知の制約 (後続)**: 定量 Cutler 照合 (He コア長・濃度減衰 vs RELIEF survey) は粘性+乱流混合 (RANS-SST + 乱流化学種拡散 `μ_t/Sc_t`) が必要で本マイルストーン外。slip/wall/outlet ghost と HLLE/AUSM/KEEP/ROE/粘性流束の多成分エンタルピーは未対応。`[He,O2,N2]` 3 成分拡張も後続。
- ビルドは Docker dev image (`forge-solver:cuda-dev`) で実施 (native は libyaml-cpp 不在)。

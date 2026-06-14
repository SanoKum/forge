# 多成分 thermally-perfect gas — 実装解説

理論は [`theory.md`](theory.md) を参照。本章は `solver_density_cuda` 上の実装対応・データ構造・カーネル責務を述べる。`thermalMethod==0` (CPG) 経路は完全保持し、`thermalMethod==2` で本モデルを有効化する。

## 1. 熱力学モジュール `cuda_forge/thermo_d.{cuh,cu}`

### データ構造
`SpeciesThermo` (POD): 分子量 `MW` [kg/mol]、Lennard-Jones パラメータ `sigma_LJ` [Å]・`eps_kB` [K]、温度域 `Tlo/Tmid/Thi`、NASA-9 係数 `low[9]/high[9]`。

### device/host 共通関数 (`__host__ __device__`, 内部 double)
- `thermo_cp_molar / thermo_h_molar` — NASA-9 評価。範囲クランプ + エンタルピー線形外挿。
- `thermo_cp_mass / thermo_h_mass / thermo_R_species` — 質量基準・比気体定数。
- 混合則 `thermo_inv_Mmix / thermo_R_mix / thermo_cp_mix / thermo_h_mix / thermo_e_mix / thermo_gamma_mix` — 質量分率 `Y[]` を引数に取る (`__host__ __device__` 共通なので初期条件 host 計算とカーネルで同一)。
- `thermo_T_from_e` — 比内部エネルギーからの Newton 温度反転 (厳密微分 $c_v$、毎反復クランプ、最大 20 反復)。
- `thermo_T_from_h` — 比エンタルピーからの反転 (Roe 平均状態の有効 γ 用)。

### DB 構築・アップロード (`thermo_d.cu`)
- 内蔵 DB に代表化学種 (N2/O2/Ar/CO2/He/AIR) の NASA-9 + LJ を保持。
- `thermo_init_db(cfg)`: `cfg.speciesNames` 順に host 配列を構築し、`cfg.speciesDBFile` (yaml) があれば上書き、device global memory (`cudaMalloc`) へアップロード。
- アクセサ `thermo_species_device_ptr()` / `thermo_num_species()` / `thermo_species_host()`。化学種データは translation unit をまたぐため `__constant__` ではなく device pointer をカーネル引数で渡す (`-rdc` 非依存)。
- `main.cpp` の `initializeSimulation` で `cfg.read()` 直後に `thermo_init_db(cfg)` を呼ぶ。

## 2. 従属変数と温度反転 `cuda_forge/dependentVariables_d.cu`

`dependentVariables_d` に `thermalMethod`・化学種ポインタ・`gam_array/cp_array` を追加。`thermalMethod==2` では:
- 組成 `Y` を構築 (`nSpecies==1` は `Y={1}`、多成分は `roY_s/ρ` をフロア+正規化)。
- `thermo_T_from_e` で温度反転 (ウォームスタートに前ステップ `T[ic]`)。
- `P=ρ R_mix T`、`Ht=h_mix+ek`、`sonic=√(γ_mix R_mix T)`、`gamma[ic]=γ_mix`、`cp[ic]=cp_mix` を更新し、`roe` を (floor 済 ρ, 反転 T) と整合再構成。

`thermalMethod==0` は従来式 (`T=e/(cp/γ)`, `P=(γ-1)(ρe-ρek)`) をビット不変で保持。CPU フォールバック (`dependentVariables.cpp`, `gpu==0`) は CPG のみ対応 (GPU 経路が主)。

物性更新 `gasProperties_d.cu` は `thermalMethod==0` のときだけ `gamma/cp` を定数で埋める。`thermalMethod==2` では `dependentVariables` が毎ステップ `gamma/cp` を設定するため触らない (粘性は viscMethod に従う。kinetic theory は M3)。

## 3. 対流流束の TP 整合 `cuda_forge/convectiveFlux_d.cu`

`SLAU_d` / `ROE_d` に `int thermalMethod, const SpeciesThermo* sp, int nSpecies` を追加。ラッパーから `cfg.thermalMethod, thermo_species_device_ptr(), cfg.nSpecies` を渡す。

### SLAU (最小修正)
被移流全エンタルピー `h_p/h_m` のみ TP 化:面温度 `T_face=P_face/(ρ_face R)` → `h(T_face)`(NASA) + `½|u|²`。音速 `c_hat=½(sonic[i]+sonic[j])` は既に `sonic[]` 配列 (混合整合) のため無改造。圧力束は `P` 直接で EOS 非依存。

**多成分 (M5)**: `h_p/h_m` を L 側 ic0 / R 側 ic1 の**セル組成 `Y`** で再構成 (`thermo_R_mix`/`thermo_h_mix`)。`sp[0]` 固定だと He/N2 のように質量基準 cp が大きく異なる contact でエネルギー束が不整合になり圧力発散する。`SLAU_d` に `flow_float** roY` (= `species_roY_device_ptr()`) を追加。`nSpecies<=1` は従来 `sp[0]` 経路でビット不変。組成は 1 次 (セル中心値) で再構成 (高次は将来)。

### ROE (段階A)
L/R 状態の `roe_L/Ht_L/ca_L` (および R 側) を NASA で再構成。Roe 平均状態は `thermo_T_from_h` で `h̃=Ha-ek` から `T̃` を反転し、有効比熱比 `ga_eff=cp(T̃)/(cp(T̃)-R)` と音速 `ca=√(ga_eff R T̃)` を用いる。圧力微分ブロック `P_ro,P_roe,...` は `κ=ga_eff-1` を採用 ($\chi=\partial P/\partial\rho|_{\rho e}\ne0$ 項は stage-A で無視。defect-correction なので収束解は不変)。厳密化は Vinokur-Montagné (stage-B, 将来)。

> SLAU は多成分対応済 (M5, セル組成再構成)。ROE/HLLE/AUSM/KEEP は単成分 (`sp[0]`) のまま (多成分化は後続)。

## 4. 陰解法 Jacobian (予定)

`timeIntegration_d.cu` の `build_jacobian_split/build_abs_jacobian` は `chi=(γ-1)/sonic` で EOS をパラメータ化する。TP では `gamma` をスカラ→`gamma[ic]` (=γ_mix) 配列読みに変えれば `κ=γ_mix-1` が厳密一致し、`sonic[ic]` は既に混合整合。frozen-coefficient (毎反復凍結再構成) で収束解は不変。M1 の検証は陽解法 RK で行うため Jacobian 変更は後続。

## 5. 設定 `input/solverConfig.{hpp,cpp}`

`physProp` に追加 (いずれも任意、`thermalMethod==2` で使用):
- `species: [N2, O2, ...]` — 混合構成 (順序が index $s$)。
- `speciesDBFile` — 外部 NASA-9/LJ DB (yaml)。空なら内蔵 DB。
- `speciesDiffusionMethod` (0:定数Sc / 1:kinetic 混合平均)、`Sc`、`Sc_t`。
  乱流シュミット数 `Sc_t` は `physProp.Sc_t` または `turbulence.turbulentSchmidt` で設定する
  (両方あれば後者を優先。`turbulence.turbulentPrandtl` と同じブロックで揃えられる)。既定 0.7。
- `nSpecies` は `species` の要素数。未指定で `thermalMethod==2` なら既定 N2 単成分。

## 5b. 多成分化学種輸送 (M2) `cuda_forge/speciesTransport_d.{cuh,cu}`

化学種は NS の 5 元ブロックには結合させず、RANS k/ω と同じ汎用スカラ輸送コア
`scalarTransport_d` (`ScalarTransportDesc`) を化学種ごとに再利用して **segregated** に解く。
保存質量分率 $\rho Y_s$ の移流方程式 $\partial_t(\rho Y_s)+\nabla\!\cdot(\rho Y_s\mathbf{u})=0$ を、
対流流束カーネルが面で確定した `massflux` (SLAU の $\dot m$) による 1 次風上で組み立てる
(M2 は移流のみ。拡散は M4)。

- **変数登録**: `variables::registerSpecies(nSpecies)` が `cellValNames`/`output_cellValNames`
  (非 const 化) へ 1 化学種あたり `roY{s},Y{s},roY{s}N,roY{s}M,res_roY{s},res_roY{s}_m,`
  `transport_diag_Y{s},src_jac_Y{s}` を追加し、`c`/`c_d` マップにも空エントリを作る。
  `allocVariables` 前に 1 度だけ呼ぶ。**`nSpecies<=1` では一切登録せず M1 と同一経路**
  (回帰がバイト一致)。
- **device roY 配列**: `speciesInit_d` が `flow_float*[nSpecies]` を 1 度構築し、
  `dependentVariables_d` (従来 `nullptr`) へ渡して混合則 thermo ($R_{mix},c_{p,mix},\gamma_{mix}$) を有効化。
- **原始量・BC**: `speciesPrimitive_d` が $Y_s=\rho Y_s/\rho$ を更新。`applySpeciesBoundaries` は
  全境界を Neumann (zero-gradient ghost) で埋める (組成依存エネルギー BC は後続)。
- **時間積分・実現可能性**: `scalarTimeIntegration_d` (floor=0) で更新後、
  `speciesRenormalize_d` が $\rho Y_s\ge0$ にクランプし $\sum_s\rho Y_s=\rho$ へ再スケール
  ($\sum_s Y_s=1$)。`roY{s}N/M` は `speciesUpdateOuter/Inner` が D2D copy で NS の N/M に同期。

## 5c. 輸送係数 — kinetic theory (M3) `cuda_forge/thermo_d.cuh` + `gasProperties_d.cu`

`thermalMethod==2` で `viscMethod==2` を選ぶと、混合粘性 $\mu$ と熱伝導率 $\lambda$ を
Chapman-Enskog + Wilke/Mason-Saxena で per-cell に評価する (LJ パラメータ `sigma_LJ`,`eps_kB` を使用)。

- **単成分** (`thermo_d.cuh`): $\mu_s=2.6693\times10^{-6}\sqrt{M_s[\mathrm{g/mol}]\,T}/(\sigma_s^2\,\Omega^{(2,2)*})$ [Pa·s]、
  $\lambda_s=\mu_s(c_{p,s}+1.25R_s)$ (modified Eucken)。衝突積分 $\Omega^{(2,2)*},\Omega^{(1,1)*}$ は
  Neufeld et al. (1972) 閉形式。
- **混合**: Wilke $\mu_{mix}=\sum_i X_i\mu_i/\sum_j X_j\phi_{ij}$、Mason-Saxena $\lambda_{mix}$ (同 $\phi_{ij}$)。
  モル分率 $X$ は `thermo_X_from_Y` で質量分率から変換。
- **per-cell `thermCond`**: `variables` に `thermCond` を登録し、`gasProperties_d` が全 viscMethod で
  値を書く (0/1 は一定 `cfg.thermCond`、2 は $\lambda_{mix}$)。`viscousFlux_d` の熱流束はスカラ→
  面平均 per-cell 配列に変更。viscMethod 0/1 では一定値のため面平均は同値で**既存挙動不変**。
- **二元/混合平均拡散** `thermo_Dbinary`/`thermo_Dmix_species` も追加 (M4 の化学種拡散で使用)。

## 5d. 化学種拡散 + エネルギー結合 (M4) `cuda_forge/speciesTransport_d.cu`

粘性ケース (`viscMethod!=0`, `nSpecies>=2`) で `species_diffusion_d` が面ループで化学種拡散を加える。

- **Fick 拡散**: $J_s=\rho D_s\,\nabla Y_s$ (over-relaxed 法線, `viscousFlux_d` と同じ幾何)。
  $D_s$ は `speciesDiffusionMethod==1` で混合平均 (`thermo_Dmix_species`)、`==0` で定数 Schmidt
  $D=\mu/(\rho\,\mathrm{Sc})$。乱流寄与 $\mu_t/(\rho\,\mathrm{Sc}_t)$ を全種に加算。
- **質量保存補正**: $J_s^\* = J_s - Y_s\sum_k J_k$ で $\sum_s J_s^\*=0$ を厳密化 (混合平均拡散は
  そのままでは $\sum J\neq0$)。`res_roY_s += J_s^\*`。
- **エネルギー結合**: 化学種拡散が運ぶエンタルピー $\sum_s h_s(T)\,J_s^\*$ を `res_roe` に加える
  ($h_s$ は NASA 絶対基準で `roe` と整合)。
- **汎用拡散との二重計上回避**: 汎用スカラコア `scalarTransport_d` の拡散 (μ ベース) は
  `ScalarTransportDesc.diffusion` フラグで制御し、化学種は `0` (移流のみ) とする。RANS k/ω は `1`。
  これを怠ると層流粘性由来の $\nu$ 拡散が Fick 拡散に重畳する (検証で発見・修正)。

## 5e. 境界条件の TP 整合 `cuda_forge/boundaryCond_d.cu`

`thermalMethod==2` で各 BC の ghost / 境界状態を NASA 熱力学で再構成する (CPG 前提の定数 γ,cp を排除)。
各カーネルは `int thermalMethod, const SpeciesThermo* sp` を受け取り、`thermalMethod==2` 分岐で
`thermo_R_species`/`thermo_cp_mass`/`thermo_h_mass`/`thermo_isentropic_from_total_single`/
`thermo_isentropic_from_total_Ps_single` を用い、`γ_mix(T)`・`R_mix`・`e_NASA(T)` で
`roe`/`sonic`/`T` を整合させる。CPG 分岐は代数的に不変。

- **M1 整合済**: `slip` / `wall`(断熱) / `inlet_Pressure` / `outflow`。
- **本節で追加**: `outlet_statPress`(順流 + backflow 等エントロピー), `wall_isothermal`(等温壁 ghost),
  `inlet_Pressure_dir`(全状態と外挿静圧から NASA 等エントロピー反転), `inlet_uniformVelocity`。
- **新規 thermo**: `thermo_isentropic_from_total_Ps_single` (全温・全圧 + 静圧 Ps から
  `s0(Ts)=s0(Tt)+R ln(Ps/Pt)` の Newton 反転で Ts, ρ, |u|)。
- **組成依存超音速入口 (M5)**: `inlet_uniformVelocity` を多成分化。`readBcondConfig` が `inlet_*`
  種別に `Y0..Y{n-1}` を type-1 uniform-float bvar として動的登録 (未指定なら `Y0=1`・他 0 で単成分後方互換)。
  `bcond::Yb_d` (device `flow_float*[nSpecies]`) を `inlet_uniformVelocity_d_wrapper` が遅延構築し、
  カーネル TP 分岐が `thermo_R_mix`/`thermo_e_mix`/`thermo_gamma_mix` で `T=P/(ρR_mix(Y))`,
  `roe=ρ(e_mix+ek)`, `sonic=√(γ_mix P/ρ)` を入口組成で再構成。化学種 ghost は
  `speciesBoundary_d_wrapper` を kind 分岐し**入口は Dirichlet** (`roY_s[ig]=ρ[ig]·Y_s^in`)、他は Neumann。
- **制約**: slip/wall/outlet の ghost roe/sonic はまだ `sp[0]` 固定 (非粘性では SLAU 流束がセル組成で
  再構成するため支配的影響なし)。`inlet_uniformVelocity` 以外の入口種別の組成指定も後続。
- **注意**: TP の `roe` は NASA 絶対基準。CPG res からの warm-start 不可。IC は `roe=ρ(e_NASA(T)+ek)` で生成する。
- **startup 注意**: 静止 IC (u=0) + `outlet_statPress` は出口が backflow 分岐に落ち、`Pt`/`Tt` 未指定で
  0/0 NaN 化する。IC に下流速度を与え出口に `Pt`/`Tt` を指定する (`case/28` 参照)。

## 6. 検証 (M1)

`case/13.nozzle_H/run_0001_tpgas_n2` (TP-N2) と `run_0002_slau_cpg` (CPG 参照) を SLAU・非粘性・warm-start で実行。`validate_tpgas.py` が中心線の M・T・P と全エンタルピー $H_0=h(T)+½|u|^2$ を抽出し、NASA(=CEA) 等エントロピー曲線と比較する。単体の NASA-9 値は NIST と一致 (N2: $R=296.8$, $c_p(300)=1040$, $c_p(2000)=1284$ J/kgK, $h(298.15)=0$、$e\to T$ 反転誤差 $\sim10^{-11}$)。

## 6b. 性能最適化 (consumer GPU の FP64 律速対策, M6)

TP-gas 経路は CPG に比べ大幅に遅い。プロファイル前の分析で**主因はアルゴリズムではなく演算精度**と判明した。

### 背景: GeForce (consumer Ampere/Ada) は FP64 が FP32 の 1/64

開発機 RTX 3060 (GA106, CC 8.6) は SM あたり FP64 コアが 2 基のみで、**倍精度スループットが単精度の 1/64**。一方 `thermo_d.cuh` の TP コアは桁落ち回避のため**内部計算が全て `double`**。CPG 経路は全て `flow_float` (FP32) の代数演算なので、同一演算でも TP は 1op あたり 32〜64 倍遅い土俵で回る。さらに `log/pow/exp` を毎セル・毎面・毎 Newton 反復で評価するため、ホットパス (`dependentVariables_d` 温度反転 / `SLAU_d` 面エンタルピー / `gasProperties_d` kinetic 輸送) で TP が支配的に重くなる。

### 施策 (M6)

精度方針を「**桁落ちのある量だけ double、それ以外は FP32**」に改める。NASA 絶対エンタルピーは生成エンタルピーを含み大きさ ~$10^7$ J/kg で増分が小さく FP32 では桁落ちするため double を保持する。一方、輸送係数・cp・混合則・面エンタルピーの増分は FP32 で安全。

1. **輸送係数の FP32 化 (①)**: `thermo_omega22/11`・`thermo_mu_species`・`thermo_lambda_species`・`thermo_wilke_phi`・`thermo_mu_mix`・`thermo_lambda_mix`・`thermo_Dbinary`・`thermo_Dmix_species` を FP32 (`expf/powf/sqrtf`) 化。これらは桁落ちが無く、`gasProperties_d` (viscMethod==2) の $O(n_{sp}^2)$ pow/exp/sqrt と `species_diffusion_d` の拡散係数評価という最重量パスを FP32 へ移す。非粘性回帰 (viscMethod==0) は不変。粘性ケースは μ/λ を NIST 値と再照合 (FP32 で十分な精度)。
2. **Newton ループの cp+h 融合評価 (②)**: `thermo_T_from_e`/`thermo_T_from_h` は毎反復 `thermo_h_mix` と `thermo_cp_mix` を**別々の種ループ**で呼び、係数選択・$T$ 冪・`log(T)` を二重計算していた。両者を 1 スイープで返す `thermo_cph_mass`/`thermo_cph_mix` を追加し、係数・$T$ 冪・`lnT` を共有。反復あたりの多項式仕事がほぼ半減。演算の並べ替えのみで結果はビット同等。`dependentVariables_d` の反転後ブロック (cp/h/gamma 再計算) にも適用。
3. **per-cell R_mix のキャッシュ (③)**: `SLAU_d` の面エンタルピーで `thermo_R_mix(YL/YR)` を**毎面で種ループ再計算**していたが、これは再構成前のセル組成だけの定数。`dependentVariables_d` で `Rmix[ic]` を 1 回計算して `gamma/cp` と並べて格納し (CPG は $R=(\gamma-1)c_p/\gamma$)、`SLAU_d` は読むだけにする。面数 ≫ セル数なので有効。

> 残レバー (M6 後続): NASA エンタルピー反転自体の FP32 化は、生成エンタルピー ≈ 0 の不活性種 (N2/He/air) では安全だが、燃焼種では sensible enthalpy ($h-h(T_{ref})$) への再定式化が必要。最大の追加効果は温度ルックアップテーブル (cp/h/s を温度グリッドで事前計算し線形補間、`log/pow/Newton` を全廃) だが結果が微変するため別マイルストーン。

## 7. マイルストーン状況

- M1 (本実装): 単成分 TP, NASA-9, Newton 反転, SLAU/ROE TP 整合 — 実装済・検証完了 (衝撃管 CEA 照合)。
- M2: 多成分輸送 (`speciesTransport_d`), 実現可能性+再正規化, 化学種 BC/IC/出力 — 実装済・検証完了。
  単成分回帰がバイト一致、[N2,N2] 識別子トレーサで流れが単成分と完全一致、[O2,N2] 実 2 ガスが安定・$\sum Y=1$。
- M3: kinetic theory 輸送係数 (`gasProperties_d` 拡張) — 実装済・検証完了。μ/λ が NIST と一致 (N2: μ(300)=1.81e-5, μ(1000)=4.15e-5, λ(300)=0.0255; air μ=1.86e-5)、device viscMethod=2 が式と FP32 一致、非粘性回帰不変。
- M4: 化学種拡散 + エンタルピー拡散エネルギー結合 — 実装済・検証完了。静止 N2/N2 トレーサの拡散が kinetic theory 自己拡散係数と 1.3% 一致、ΣJ=0 で種質量保存 (7 桁)、O2/N2 粘性で安定。
- M5: 組成依存超音速入口 BC (`inlet_uniformVelocity` 多成分化) + SLAU 対流流束の多成分エンタルピー — 実装済・検証完了。`case/28.cutler_coaxial_jet` の binary [He,N2] 同軸ジェット (Mach1.8) が 1次/2次でともに NaN/発散なし、$\sum Y=1$、中心 He コア (Mach1.8 He, Ux~1280 m/s)・coflow N2 を再現。CPG 単成分回帰不変。定量 Cutler 照合 (粘性 + RANS-SST + 乱流化学種拡散) は後続。
- M6: TP-gas 性能最適化 (§6b, 輸送係数 FP32 化 / cp+h 融合反転 / per-cell `Rmix` キャッシュ) — 実装済・検証完了。① CPG 単成分回帰は M6 前後バイナリの差が CUDA atomicAdd 非決定性 (旧×旧の run 間差) と同等で系統的変化なし。② FP32 輸送が double と相対 ~1e-8 一致 (host 単体テスト `solver_density_cuda/tools/test_m6_transport.cpp`)、NIST と N2 μ 1.5%/λ 1.9%・AIR μ 0.3%/λ 4.3%・He μ 1.7%。③ TP binary [He,N2] スモーク (`case/28` `run_0026_he_n2_m6smoke`, 500 步) が NaN なし・$\sum Y=1$ (1.2e-7)・$0\le Y\le1$・ro/P/T 物理的。速度 (native RTX 3060, 76k セル発達場 per-step): 3 成分粘性 **1.56×** (46.07→29.56 ms)、binary 粘性 1.31×、binary 非粘性 1.02× — 支配項は FP32 輸送 (viscMethod=2 × $O(n_{sp}^2)$)。詳細は plan 変更ログ 2026-06-10。
- M7: zero-gradient/反射系 BC の組成依存エネルギー化 (`slip`/`wall`/`wall_isothermal`/`outlet_statPress` の ghost `R/e/cp/γ` を `sp[0]` 固定→内部セル混合則) — 実装済・検証完了。`bc_cell_Y` ヘルパ + `thermo_isentropic_from_total_mix`/`thermo_s0_mix` 追加。これまで Cutler (species[0]=He) で空気/coflow 領域の壁・出口が He 物性で熱再構成されるバグを解消。単成分 TP 回帰はカオス的非決定性内で不変、Cutler binary [He,N2] スモークが NaN なし・$\sum Y=1$。`inlet_Pressure`/`inlet_Pressure_dir` (指定組成 `Yb` 要・多成分未使用) は後続。詳細は plan 変更ログ 2026-06-10。

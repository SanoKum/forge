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

### ROE (段階A)
L/R 状態の `roe_L/Ht_L/ca_L` (および R 側) を NASA で再構成。Roe 平均状態は `thermo_T_from_h` で `h̃=Ha-ek` から `T̃` を反転し、有効比熱比 `ga_eff=cp(T̃)/(cp(T̃)-R)` と音速 `ca=√(ga_eff R T̃)` を用いる。圧力微分ブロック `P_ro,P_roe,...` は `κ=ga_eff-1` を採用 ($\chi=\partial P/\partial\rho|_{\rho e}\ne0$ 項は stage-A で無視。defect-correction なので収束解は不変)。厳密化は Vinokur-Montagné (stage-B, 将来)。

> 現状 M1 は単成分 (`sp[0]`)。多成分は面で `Y_s` を再構成し `R_mix`/`h_mix` を評価する形へ拡張する。

## 4. 陰解法 Jacobian (予定)

`timeIntegration_d.cu` の `build_jacobian_split/build_abs_jacobian` は `chi=(γ-1)/sonic` で EOS をパラメータ化する。TP では `gamma` をスカラ→`gamma[ic]` (=γ_mix) 配列読みに変えれば `κ=γ_mix-1` が厳密一致し、`sonic[ic]` は既に混合整合。frozen-coefficient (毎反復凍結再構成) で収束解は不変。M1 の検証は陽解法 RK で行うため Jacobian 変更は後続。

## 5. 設定 `input/solverConfig.{hpp,cpp}`

`physProp` に追加 (いずれも任意、`thermalMethod==2` で使用):
- `species: [N2, O2, ...]` — 混合構成 (順序が index $s$)。
- `speciesDBFile` — 外部 NASA-9/LJ DB (yaml)。空なら内蔵 DB。
- `speciesDiffusionMethod` (0:定数Sc / 1:kinetic 混合平均)、`Sc`、`Sc_t`。
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

## 6. 検証 (M1)

`case/13.nozzle_H/run_0001_tpgas_n2` (TP-N2) と `run_0002_slau_cpg` (CPG 参照) を SLAU・非粘性・warm-start で実行。`validate_tpgas.py` が中心線の M・T・P と全エンタルピー $H_0=h(T)+½|u|^2$ を抽出し、NASA(=CEA) 等エントロピー曲線と比較する。単体の NASA-9 値は NIST と一致 (N2: $R=296.8$, $c_p(300)=1040$, $c_p(2000)=1284$ J/kgK, $h(298.15)=0$、$e\to T$ 反転誤差 $\sim10^{-11}$)。

## 7. マイルストーン状況

- M1 (本実装): 単成分 TP, NASA-9, Newton 反転, SLAU/ROE TP 整合 — 実装済・検証完了 (衝撃管 CEA 照合)。
- M2: 多成分輸送 (`speciesTransport_d`), 実現可能性+再正規化, 化学種 BC/IC/出力 — 実装済・検証完了。
  単成分回帰がバイト一致、[N2,N2] 識別子トレーサで流れが単成分と完全一致、[O2,N2] 実 2 ガスが安定・$\sum Y=1$。
- M3: kinetic theory 輸送係数 (`gasProperties_d` 拡張)。
- M4: 化学種拡散 + エンタルピー拡散エネルギー結合。

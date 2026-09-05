# 有限速度化学 (finite-rate chemistry)

多成分 thermally-perfect gas ([thermophysics.md](thermophysics.md)) に**化学反応ソース項**を加え、H₂ 燃焼 (燃焼加熱器) と
ノズル膨張中の化学非平衡 (ラジカル再結合と凍結) を解く。本ドキュメントは理論 (反応速度・平衡定数・反応熱・Jacobian) と
実装 (カーネル・陰解法結合・入力配管) をまとめる。設計判断と進捗は
[`plans/active/chemistry-finite-rate-h2.md`](../plans/active/chemistry-finite-rate-h2.md)、
文献調査は [`notes/investigations/chemistry-finite-rate-h2-survey.md`](../notes/investigations/chemistry-finite-rate-h2-survey.md)。

**状態 (2026-09-05)**: Phase 0–2 (熱力学 DB ツール・機構ファイル・反応ソース項・解析 Jacobian・種ブロック point-implicit・
反応熱の陰的注入・falloff (Troe)・Strang 分離・Jacobian 凍結・PaSR) を実装・検証済。0-D 着火 (case/35) と Q1D ノズル再結合
(case/46) は $Y_{OH}$/$Y_{NO}$ が Cantera PFR と一致 (T/M の絶対比較と `check_convergence` PASS は未達)。Phase 3 (燃焼器 RANS) は Burrows–Kurkov (case/47) で出口組成が実験と一致する一方、着火位置は
上流に過大 (5–9 cm vs 18–25 cm)、Cabra (case/48) は laminar-chemistry / PaSR ともリップ付着 (実験 H/d≈10)。Cabra は
機構交換 (Li 2004) でも付着が変わらず、**TCI 欠如が主因という仮説を強く支持** (plan §9 2026-09-05; BK の機構 A/B は未)。TCI の方針は
[`notes/investigations/cabra-liftoff-model-fidelity-survey.md`](../notes/investigations/cabra-liftoff-model-fidelity-survey.md)。

---

## 理論

### 1. 支配方程式

化学種 $s$ の保存形輸送方程式 (非反応形は thermophysics.md §3) に生成率 $\dot\omega_s$ [kg/m³/s] を加える:

$$
\frac{\partial \rho Y_s}{\partial t}+\nabla\!\cdot(\rho Y_s\mathbf u+\mathbf J_s)=\dot\omega_s,\qquad \sum_s\dot\omega_s=0
$$

全質量・運動量は不変。エネルギー式は絶対エンタルピー基準 (NASA-9 の $a_7$ が生成エンタルピーを含む) では
ソース項を持たないが、forge は陰解法安定化のため **sensible-enthalpy datum** (`thermoHrefTemp`, 各種 $h_s(T_{\rm ref})=0$ に平行移動) を
使う。平行移動 $c_s=h^{\rm abs}_s(T_{\rm ref})$ は反応が無ければ厳密に相殺する (thermophysics.md §2.2) が、反応があると
$\sum_s c_s\,\partial_t(\rho Y_s)$ の相殺が崩れ、残差として反応熱が現れる:

$$
\frac{\partial \rho e_t}{\partial t}+\nabla\!\cdot(\dots)=\dot Q,\qquad
\dot Q=-\sum_s c_s\,\dot\omega_s=-\sum_s h^{\rm abs}_s(T_{\rm ref})\,\dot\omega_s
$$

$c_s$ は生成エンタルピー $\Delta h^\circ_{f,s}$ と顕熱 $h_s(T_{\rm ref})-h_s(298.15)$ の和で、$T_{\rm ref}=298.15$ K なら $c_s=\Delta h^\circ_{f,s}$。
`thermoHrefTemp: 0` (絶対基準) では $c_s=0$ で $\dot Q=0$ となり、同じコードで両基準が整合する。
**反応流では `thermoHrefTemp: 298.15` を標準**とし、IC/BC 生成器も同じ datum で `roe` を作る (非反応と同じ規約)。

### 2. 反応速度

反応 $r$: $\sum_s\nu'_{sr}X_s\rightleftharpoons\sum_s\nu''_{sr}X_s$ の正味進行率 [mol/m³/s]

$$
q_r=k_{f,r}\prod_s[X_s]^{\nu'_{sr}}-k_{b,r}\prod_s[X_s]^{\nu''_{sr}},\qquad [X_s]=\frac{\rho Y_s}{W_s},\qquad
\dot\omega_s=W_s\sum_r(\nu''_{sr}-\nu'_{sr})\,q_r
$$

- **修正 Arrhenius**: $k_f=A\,T^{b}\exp(-E_a/R_uT)$。機構ファイルの単位 (cm³, mol, s, cal/mol) は読込時に SI へ換算する。
- **三体反応** (`type: three-body`): $q_r$ に $[M]=\sum_s\alpha_{sr}[X_s]$ を乗じる ($\alpha_{sr}$ は種別効率、既定 1)。
- **Falloff** (`type: falloff`, Lindemann / Troe): $k=k_\infty\,\dfrac{P_r}{1+P_r}F$, $P_r=k_0[M]/k_\infty$。Phase 2 で追加
  (Jachimowski 機構は falloff を含まない。Burke 2012 の H+O₂(+M) で必要)。
- **逆反応**: 平衡定数から $k_b=k_f/K_c$、
  $$K_c=K_p\Big(\frac{p^\circ}{R_uT}\Big)^{\sum_s(\nu''_{sr}-\nu'_{sr})},\qquad
  K_p=\exp\!\Big(-\frac{\Delta G^\circ_r}{R_uT}\Big)=\exp\!\Big(\sum_s(\nu''_{sr}-\nu'_{sr})\big(\tfrac{S^\circ_s}{R_u}-\tfrac{H_s}{R_uT}\big)\Big)$$
  $H_s/R_uT$, $S^\circ_s/R_u$ は NASA-9 (thermophysics.md §2.2) から。**$\Delta G^\circ$ は絶対基準の $a_7$ で評価する**
  (sensible datum を焼き込んだ $a_7$ を使うと $\Delta H^\circ_r$ から生成エンタルピー差が消えて $K_c$ が狂う)。
  そのため `SpeciesThermo::h_datum` に焼き込みで除いた $h^{abs}_s(T_{ref})$ [J/mol] を保持し、$H_s^{abs}=h_s+h_{datum,s}$ で $\Delta G^\circ$ を組む。

### 3. 化学 Jacobian (解析)

保存量 $U=(\rho Y_1,\dots,\rho Y_{n_s},\rho e_t)$ に対し、セル内の $(\rho,\,\mathbf u)$ を固定して

$$
\frac{\partial\dot\omega_s}{\partial(\rho Y_k)}=
\underbrace{\frac{\partial\dot\omega_s}{\partial[X_k]}\frac{1}{W_k}}_{\text{濃度項}}
+\underbrace{\frac{\partial\dot\omega_s}{\partial T}\frac{\partial T}{\partial(\rho Y_k)}}_{\text{温度項}},\qquad
\frac{\partial T}{\partial(\rho Y_k)}=-\frac{e_k(T)}{\rho c_{v,\rm mix}},\quad
\frac{\partial T}{\partial(\rho e_t)}=\frac{1}{\rho c_{v,\rm mix}}
$$

濃度項は質量作用則の微分 ($\partial q_r/\partial[X_k]=k_f\nu'_{kr}[X_k]^{-1}\prod[X]^{\nu'}-k_b\nu''_{kr}[X_k]^{-1}\prod[X]^{\nu''}$ + 三体 $[M]$ の $\alpha_{kr}$ 項)、
温度項は $\partial k_f/\partial T=k_f(b/T+E_a/R_uT^2)$ と van 't Hoff $\partial\ln K_c/\partial T=\Delta H^\circ_r/R_uT^2-\Delta\nu_r/T$。
$\partial T/\partial(\rho Y_k)$ は既存の EOS クロス応答 (`species_eos_cross_response`, thermophysics.md 実装 §4b) と同式。
反応熱の温度応答 $\partial\dot Q/\partial(\rho e_t)=-\sum_s c_s\,\partial\dot\omega_s/\partial T\,/(\rho c_v)$ は流れブロックの $(5,5)$ 対角に入れる。

Jacobian は $n_s\times n_s$ 密 (H₂ 系で $n_s\le13$)。有限差分は検証ハーネス専用とする。

### 4. 剛性と時間積分

化学時間 $\tau_c=1/\max|\lambda(\partial\dot\omega/\partial U)|$ は燃焼域で $10^{-9}$–$10^{-8}$ s、擬似時間 $\Delta t$ (`cfl_pseudo` 2–8) は $10^{-7}$–$10^{-6}$ s。
剛性比 $10^1$–$10^3$ のため、ソース項は**必ず陰的**に扱う。

- **定常 (block-DPLUR, 主経路)** — *decoupled point-implicit*: 流れ 5 変数を block-DPLUR で解いた後、化学種を
  $$\Big[\Big(\tfrac{V}{\Delta t}+D^{\rm conv}_s\Big)\delta_{sk}-V\frac{\partial\dot\omega_s}{\partial(\rho Y_k)}\Big]\delta(\rho Y_k)=R_s+(\text{隣接 DPLUR 項})$$
  の $n_s\times n_s$ ブロックをセル毎に LU で解く (Candler et al. 2013 と同型)。既存の scalar-DPLUR sweep の対角 `src_jac_Y{s}` を
  ブロックに拡張した形。案C (`speciesImplicitCoupling=2`) の予測 $\delta(\rho Y)^*$ に反応寄与が入るので、流れ側は EOS-JVP 経由で反応熱を
  同一 block 解の中で見る。
- **非定常 (dual-time / 陽解法 RK, 副経路)** — Strang 分離: 対流拡散半ステップ → セル内 ODE $d(\rho Y)/dt=\dot\omega$ を backward Euler
  で sub-cycle (同じ Jacobian・LU) → 対流拡散半ステップ。sub-cycle 数は $\Delta t/\tau_c$ (Gershgorin 近似) から適応。
- **反応熱の陰的注入 (Phase 2a, 2026-09-04)**: 反応熱 $\dot Q^n$ を陽に流れブロックへ入れると、燃焼室組成のままのスロート付近
  (再結合熱 $\dot Q\sim7\times10^{11}$ W/m³, $\tau_c\sim4\times10^{-8}$ s) では擬似時間刻み $10^{-6}$ s で 1 step に内部エネルギーの
  数十 % が加わり数 step で発散する (case/46 run_0002–0007: 結合方式・`cfl_pseudo` 0.2 でも同じ)。案C の予測子 $\delta(\rho Y)^*$
  (種ブロック陰解の出力) を使って
  $$\dot Q^*=\dot Q^n+\sum_k\frac{\partial\dot Q}{\partial(\rho Y_k)}\,\delta(\rho Y_k)^*,\qquad
    \frac{\partial\dot Q}{\partial(\rho Y_k)}=-\sum_s c_s\,J^{\rm tot}_{sk}$$
  と線形化し、差分 $V(\dot Q^*-\dot Q^n)$ を EOS-JVP と同じ場所で `res_roe` に加える。これは離散恒等式
  $\Delta E=-\sum_s c_s\,\Delta(\rho Y_s)$ (ホストテストの反応器と同じ) の線形化で、化学種が陰的に消費した分だけの熱を流れが見るので
  過大な発熱を出さない。**反応流の定常陰解法は `speciesImplicitCoupling: 2` + `jacobianMode: 2` を必須**とする
  (coupling 0/1 では予測子が無く陽的注入のまま)。非定常陽解法 (RK, 小さい $\Delta t$) は陽的注入で問題ない (case/35 run_0049)。
- **正値性**: 消費項 (対角負) を陰、生成項を陽に置く Patankar 型で $\rho Y_s\ge0$ を保つ。残りは既存の `species_renormalize_d` でクリップ・再正規化。

### 5. 乱流‐化学相互作用 (TCI)

既定は **laminar finite-rate (No-TCI)**。`chemistry.tci: 1` で PaSR: $\dot\omega_s\to\kappa\,\dot\omega_s$ (Jacobian・$\dot Q$ も同じ κ),
$\kappa=\tau_c/(\tau_c+\tau_{\rm mix})$, $\tau_c$ は `tciTauChem` で選択 (0: $1/\max_s|\partial\dot\omega_s/\partial\rho Y_s|$ = 最速ラジカル時間 ~1e-8 s で κ≈0.03 となり Burrows–Kurkov で消炎 [run_0018]; 1 [既定]: 燃料/酸化剤 H₂・O₂ の消費時間 $\max(\rho Y_s/|\dot\omega_s|)$),
$\tau_{\rm mix}=C_{\rm mix}\sqrt{\nu/\varepsilon}$ (`tciMixModel: 0`) または $C_{\rm mix}k/\varepsilon$ (1), $\varepsilon=\beta^*k\omega$ (SST)。
κ の状態微分は Jacobian に入れない。実装は `chemistry_source_d` (ソース項経路のみ。Strang 経路は κ=1)。$C_{\rm mix}$ は未較正。**Cabra (case/48) で PaSR は C_mix 1/4 とも付着推移を変えない** — 自着火安定化火炎は化学律速で、混合速度リミッタでは着火遅れ統計 (最反応性混合分率・低散逸ポケット) を表現できない。EDC・flamelet は採らない。次段として **1st-order CMC を導入する** (`tci: 2`, 設計は [chemistry_cmc.md](chemistry_cmc.md), 計画 `plans/active/chemistry-cmc-tci.md`)。

### 参考文献

- Jachimowski, NASA TP-2791 (1988) — H₂/air 13 種 33 反応 (Table I)。
- Burke, Chaos, Ju, Dryer, Klippenstein, Int. J. Chem. Kinet. 44 (2012) — 高圧 H₂/O₂、Troe falloff。
- Bussing & Murman, AIAA J. 26 (1988); Eberhardt & Imlay, JTHT 6 (1992) — point-implicit ソース項。
- Candler, Subbareddy, Nompelis, AIAA J. (2013) — decoupled implicit。
- Niemeyer, Curtis, Sung, CPC 215 (2017) (pyJac); Perini, Galligani, Reitz, Energy & Fuels (2012) — 解析 Jacobian。
- 詳細は調査ノート §8。

---

## 実装

### 1. 入力配管 (実装済)

- **熱力学 DB**: `solver_density_cuda/tools/cea_thermo_to_species_db.py thermo.inp --species ... --out species_db.yaml`
  が NASA CEA `thermo.inp` から forge の `speciesDBFile` を生成する (NASA-9 2 温度域、LJ パラメータは Cantera h2o2/gri30 の
  transport 値の内蔵表)。内蔵 N2 係数との照合 `--check` で相対差 0。生成物の例:
  [`tools/mechanisms/species_db_h2air_cea.yaml`](../solver_density_cuda/tools/mechanisms/species_db_h2air_cea.yaml)
  (H2 O2 H O OH HO2 H2O2 H2O N2 N NO NO2 HNO AR)。種名キーは必ず引用符付き (`"NO"`, `"N"` は YAML 1.1 で真偽値に化ける)。
  Cantera との照合: $c_p/R$ は 4 桁一致、$h/RT$ は OH/HO₂ で生成エンタルピーの出典差 (CEA 2002 Ruscic vs GRI) により 2–5 % 違う
  (forge は CEA 側で統一)。
- **反応機構ファイル**: Cantera YAML 形式のサブセット (`units`, `phases[0].species`, `species[].thermo (NASA9)`, `reactions[]` の
  `equation` / `rate-constant {A,b,Ea}` / `type: three-body` / `efficiencies`、Phase 2 で `type: falloff` + `Troe`)。
  同梱: [`tools/mechanisms/h2air_jachimowski1988_13sp33r.yaml`](../solver_density_cuda/tools/mechanisms/h2air_jachimowski1988_13sp33r.yaml)
  (Table I 全 33 反応) と [`h2o2_jachimowski1988_9sp20r.yaml`](../solver_density_cuda/tools/mechanisms/h2o2_jachimowski1988_9sp20r.yaml)
  (反応 (1)–(20)、N₂ 不活性)。**同じファイルを Cantera が読める**ので、0-D/1-D の参照解は Cantera (CPU, `.venv-chem`) で作る。
- `solverConfig.yaml`: `physProp.chemistry: {enabled: 1, mechanismFile: ..., jacobianMode: 1, tMaxReaction: 6000, freezeBelowT: 0}`
  (キーの説明は [procedures/solver-settings.md](../procedures/solver-settings.md))。PaSR の `tci` は Phase 3。

### 2. カーネル構成 (Phase 1–2 実装済)

| 処理 | ファイル | 内容 | 状態 |
| --- | --- | --- | --- |
| 反応表 + 生成率評価 (host/device 共通) | `cuda_forge/chemistry_d.cuh` (`ReactionTable`, `chem_source`) | 固定長 POD の反応表。`chem_source(sp, rt, ρ, T, Y, ω, Q̇, jacMode, J, ∂ω/∂T)` が Arrhenius・三体・$K_c$ 逆反応・反応熱・解析 Jacobian (`jacMode` 0/1 対角/2 全 $n_s\times n_s$) を double で評価 | 済 |
| 機構ファイル読込 (host) | `cuda_forge/chemistry_mech_io.hpp` (`chem_io::loadMechanism`) | Cantera YAML サブセット → `ReactionTable`。単位換算 (cm³/mol/s, cal/mol → SI)、種名を `physProp.species` 順に対応付け。falloff (`high-P-rate-constant`/`low-P-rate-constant`/`Troe`, T2 省略可) も読む。Cantera `ck2yaml` 出力 (Li 2004, Troe 2 パラメータ形) をそのまま読めることを case/48 `run_0072` で確認 (Burke 2012 は未実行) | 済 |
| device 側ソース項 | `cuda_forge/chemistry_d.cu` (`chemistry_init`, `chemistry_source_d`, `chemistrySource_d_wrapper`) | セル毎に `res_roY{s} += Vω_s`, `res_roe += VQ̇`, `src_jac_Y{s} += max(0,−∂ω_s/∂ρY_s)`, 診断 `chemQdot`/`chemTau`。周期 node では部分体積 `volumePartial_d` を使う (seam 二重計上防止) | 済 |
| host インタフェース | `cuda_forge/chemistrySource_d.cuh` | `main.cpp`: `thermo_init_db` 直後に `chemistry_init`、`assembleResidual` の `speciesTransport_d_wrapper` 直後に `chemistrySource_d_wrapper` | 済 |
| datum 保持 | `thermo_d.{cuh,cu}` (`SpeciesThermo::h_datum`) | sensible datum で除いた $h^{abs}_s(T_{ref})$ [J/mol] を保持し、$\dot Q$ と $K_c$ (絶対 $H$) に使う。**構造体が変わったので full rebuild 必須** | 済 |
| 設定・変数 | `input/solverConfig.{hpp,cpp}` (`physProp.chemistry`), `variables.cpp` (`registerSpecies(n, chemistry)`) | キーは [solver-settings.md](../procedures/solver-settings.md) | 済 |
| node 周期の化学種残差 gather | `cuda_forge/periodicNode_d.cu` | `periodicNodeGather` が `res_roY{s}` (と凝縮モーメント残差) を合算し、`periodicMirrorNSState` が `roY{s}` をミラーする。**従来は流れ 5 変数と k/ω だけで、seam の化学種は部分体積分しか更新されなかった** (化学以前からの欠落。run_0049 で発覚・修正) | 済 |
| 種ブロック point-implicit (`jacobianMode: 2`) | `speciesTransport_d.cu` `species_dplur_block_sweep_d` (`speciesSweepOnce` が block/scalar を切替) | chemistry が `chem_jac` (温度結合込み $J^{\rm tot}$ の残差行列 $R=J^{\rm tot}+\mathrm{diag}(d)$, $d_s=\max(0,-J_{ss})$ は `src_jac_Y` へ) を書き、sweep はセル毎に $[(V/\Delta\tau+D^{\rm conv}+V\,\mathrm{src\_jac})I-VR]\delta=\mathrm{rhs}$ を部分ピボット LU (double) で解く。coupling 1 の解と coupling 2 の予測子の両方で使う | 済 (Phase 2a) |
| 反応熱の陰的注入 | `chemistry_heat_inject_d` (speciesTransport_d.cu, 案C 予測子直後) + `chem_cq` ($\partial\dot Q/\partial\rho Y_k$) | 理論 §4。`res_roe += V\sum_k(\partial\dot Q/\partial\rho Y_k)\delta(\rho Y_k)^*` | 済 (Phase 2a) |
| 反応熱のエネルギー対角 | `timeIntegration_d.cu` `implicit_defect_correction_block_d` (`src_jac_roe` 引数) | $(4,4)\ +=\ V\max(0,-\partial\dot Q/\partial\rho e)$ (precond 版 `_precond_d` は配線済: 2026-09-05 コード確認) | 済 (Phase 2a) |
| Jacobian 凍結 (`jacobianInterval: n`) | `chemistry_d.cu` (`frozenJac`, `chem_diag`) | 定常陰解法で n ステップに 1 回だけ Jacobian (`chem_jac`/`chem_cq`/`chem_jacroe`/対角 `chem_diag`) を再評価。凍結ステップは `chem_source(mode 0)` で ω・Q̇ だけ計算し、前回の対角を `src_jac_Y` に足す。κ (PaSR) も前回値。種ブロック LU (`species_block_factor_d`) は dt_local が変わるため毎ステップ分解 (sweep 間では再利用) | 済 (高速化 2) |
| Strang 分離 (`strang: 1`) | `chemistry_strang_d` / `chemistryStrangHalfStep_d_wrapper` (chemistry_d.cu), `advanceExplicitRK` (main.cpp) | 非定常 RK の前後で dt/2 ずつセル内 ODE ($d\rho Y/dt=\omega$, $\Delta\rho e=-\sum c_s\Delta\rho Y_s$) を backward Euler sub-cycle (v2: 刻みは「主要種の相対変化 <20 %・|ΔT|<30 K」で受理、超えたら半分にして再試行、受理が続けば倍に戻す。非反応セルは 1 sub-step・Newton 1 回。最大 64 sub-step、温度は $e$ から反転)。前半の後に N/M 始点を取り直す。有効時はソース項を RK に入れない | 済 (Phase 2c) |
| Falloff (Lindemann/Troe) | `chem_ln_kf_falloff` (chemistry_d.cuh), `chemistry_mech_io.hpp` (`high-P-rate-constant`/`low-P-rate-constant`/`Troe`) | $k=k_\infty\frac{P_r}{1+P_r}F$。$\partial\ln k/\partial T$, $\partial\ln k/\partial[M]$ は滑らかなスカラ関数の中心差分 (相対 1e-4)。GRI H₂/O₂ (`tools/mechanisms/h2o2_gri30_nasa9_10sp29r.yaml`) で Jacobian FD 2e-8, 着火遅れ 43.6 vs Cantera 44.0 µs | 済 (Phase 2b) |
| ホスト単体検証 | `tools/test_chemistry.cpp` + `tools/chem_reference_cantera.py` | Jacobian 有限差分照合 (相対 2e-8 PASS)、0-D 定積反応器 BDF1 が Cantera と一致 (着火遅れ 32.0 vs 32.2 µs, 平衡 T +1.3 K)、絶対/sensible datum の一致 | 済 |

ビルド: `g++ -O2 -std=c++17 -I solver_density_cuda solver_density_cuda/tools/test_chemistry.cpp -lyaml-cpp -o /tmp/tchem` →
`/tmp/tchem tools/mechanisms/h2o2_jachimowski1988_9sp20r.yaml tools/mechanisms/species_db_h2air_cea.yaml ref.csv 1200 1 5e-4`
(`ref.csv` は `.venv-chem/bin/python tools/chem_reference_cantera.py MECH.yaml ref.csv 1200 1 5e-4`)。

### 3. 精度方針

濃度・質量作用則の積・Jacobian 解は double、$\ln k_f$ と $\ln K_c$ の温度依存項は FP32 で評価して $\exp$ は 1 回
(thermophysics.md 実装 §6b の「桁落ちする量だけ double」方針)。

### 4. 検証

1. **0-D 定積着火 (済, 2026-09-04)**: ホストテスト (上記) と GPU 周期箱
   `case/35.uniform_periodic_box/run_0049_node_h2_ignition` (node 8³ 全面 periodic, 量論 H₂-air 1200 K 1 atm, explicit RK3
   dt=5e-9 s, jacobianMode=1)。全 729 ノードが一様に着火 (T 差 0.004 K, |u|<1e-4 m/s)、着火遅れ 32.0 µs (Cantera 32.2)、
   平衡温度 2948 K (Cantera 2944, +0.15 %)。図 `h2_ignition_vs_cantera.png`。
2. **Q1D ノズル再結合 (済)**: `case/46.nozzle_h2o2_kinetics run_0008` — $Y_{OH}$ が Cantera PFR (面積則) と一致、frozen 比で
   出口 $T$ +58 K・$M$ −1.1 %。反応熱を陽的に入れると定常陰解法はスロートで発散するため `speciesImplicitCoupling: 2` +
   `jacobianMode: 2` が必須 (plan §9 2026-09-04)。`check_convergence` は NOT CONVERGED (rms_roUy stalled) のまま。
3. **Burrows–Kurkov (一次結果)**: `case/47 run_0025` (メッシュ v2, 入口 BL δ 10 mm, 入口 k×3): 出口 $X_{H_2O}$ ピーク 0.505 @ 1.96 cm
   (実験 0.50 @ 2.0)、全温ピーク 1.08 vs 1.18、着火 4–6 cm (実験 18–25 cm, 未解決)。3 バグ修正後の再確認 `run_0028` で結論維持。
4. **Cabra H₂/N₂ (一次結果)**: `case/48 run_0069`–`run_0073` — 自着火 x/d≈20 → 上流伝播 → リップ付着 (実験 H/d≈10)。PaSR C_mix・
   機構 (Li 2004) とも非感応。文献どおり laminar-chemistry の構造的失敗であり、H/d 再現の実績が最も強いのは transported PDF / CMC (調査メモ参照)。
   準定常判定は `tools/check_quasisteady.py --quantity ignx,tmax,exit_massflux,exit_hflux,exit_y_<種>` (2026-09-05 追加)。
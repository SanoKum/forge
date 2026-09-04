# 有限速度化学 (H₂ 燃焼・ノズル化学非平衡) の導入

## メタ

- **area**: `thermophysics / chemistry`
- **status**: `in_progress` (Phase 0–2 完了、Phase 3 Burrows–Kurkov 一次結果あり 2026-09-04)
- **related_docs**:
  - `methods/chemistry.md` (現在仕様: 理論・実装方針)
  - `methods/thermophysics.md` (多成分 TP gas、sensible datum、種輸送・陰解法結合)
  - `notes/investigations/chemistry-finite-rate-h2-survey.md` (文献調査・方針の根拠)
- **related_plans**:
  - `plans/accepted/thermophysics-multicomponent-tpgas.md` (前提。§10 に反応は将来拡張と明記)
  - `plans/accepted/convection-multispecies-contact-pressure.md`
- **created**: `2026-09-04`
- **owner**: `CFD Dev`

## 1. 目的

forge の多成分 TP gas に化学反応ソース項を加え、(B) ノズル膨張中のラジカル再結合・凍結 (試験部の残留 OH/H/O/NO と $T$, $M$ への影響)、
(A) H₂ 燃焼加熱器 (低マッハ・RANS) を解けるようにする。完了時には frozen / equilibrium / finite-rate の三者比較が
ノズル設計チェーン (`methods/design/`) の確認 CFD メニューとして使える状態になる。

## 2. スコープ

- **やる**: 反応機構 (Cantera YAML サブセット) 読込、Arrhenius・三体・逆反応 ($K_c$)、反応熱の sensible datum 整合、解析 Jacobian、
  種ブロック point-implicit (定常) と Strang 分離 (非定常)、0-D / Q1D ノズル / Burrows–Kurkov 検証、No-TCI → PaSR。
- **やらない**: fully coupled $(5+n_s)$ block-DPLUR、EDC / flamelet、熱的非平衡 (多温度)、炭化水素機構 (CH₄ 加熱器は後続)、
  燃焼器の着火過渡 (定常は高温平衡 IC から始める)。

## 3. 関連 docs と前提

- 理論・実装方針: `methods/chemistry.md` (本 plan と同時に起草)。
- 既存資産: `speciesTransport_d` (種輸送・拡散・`src_jac_Y{s}` 配管・scalar-DPLUR sweep・案C EOS-JVP)、
  `condensationSource_d` / `ransSource_d` (point-implicit ソース項の雛型)、`thermo_d` (NASA-9, `THERMO_HD` host/device 共通)。
- 既知の制約: 多成分 TP 陰解法の安定 `cfl_pseudo` は 1–2 (thermophysics plan §10)。反応で剛性が上乗せされる。

## 4. 設計方針 (差分のみ。本文は `methods/chemistry.md`)

1. **datum**: `thermoHrefTemp: 298.15` を維持し、$\dot Q=-\sum_s c_s\dot\omega_s$ ($c_s=h^{\rm abs}_s(T_{\rm ref})$) をエネルギー式に陽注入。
   $K_c$ は焼き込み前の絶対 $a_7$ (`a7_abs`) で評価。
2. **剛性**: 既存 `src_jac_Y{s}` (対角) を $n_s\times n_s$ 密ブロック `src_jac_Y{s,k}` に拡張し `species_dplur_sweep_d` でセル毎 LU。
   流れ側は案C の EOS-JVP に反応寄与を載せ、$(5,5)$ 対角に $\partial\dot Q/\partial(\rho e_t)$。
3. **機構**: Jachimowski 1988 (9 種 20 反応 / 13 種 33 反応、falloff 無し) から。Burke 2012 (Troe) は Phase 2 で falloff 実装後。
4. **精度**: 濃度・Jacobian は double、$\ln k_f$/$\ln K_c$ は FP32。
5. **参照解**: Cantera (`.venv-chem`, CPU) と CEA (`.venv-cea`)。forge と同一 YAML・同一 NASA-9 を使う。

## 5. 実装ステップ

| Phase | 内容 | 主要ファイル | 状態 |
| --- | --- | --- | --- |
| **0** 前提整備 | CEA 平衡/凍結スクリーニング、`thermo.inp→species_db.yaml` ツール、Jachimowski YAML 化 + Cantera 検証、方針文書 | `tools/cea_thermo_to_species_db.py`, `tools/mechanisms/*.yaml`, `notes/investigations/scripts/cea_kinetics_screen/`, `methods/chemistry.md` | **done** (§9) |
| **0b** BC | 指定組成 `Yb` 配線: `inlet_Pressure` は M5/M7 で配線済 (case/44 `run_0091` が `Y0/Y1` を使用中、thermophysics plan §10 の記述は古い)。残りは `inlet_Pressure_dir` のみ (燃焼室出口組成は `inlet_Pressure` で与えられるので低優先) | `cuda_forge/boundaryCond_d.cu` | `inlet_Pressure` done / `_dir` todo |
| **1** ソース項 | `chemistry_d.{cuh,cu}` (機構読込・SI 換算・device 定数)、`chemistrySource_d` ($\dot\omega_s$, $\dot Q$, 解析 Jacobian)、`h_datum` 保持、config キー、ホスト 0-D テスト `tools/test_chemistry.cpp`、GPU 0-D 着火 (case/35 run_0049) | `cuda_forge/chemistry_d.*`, `chemistry_mech_io.hpp`, `thermo_d.*`, `input/solverConfig.*`, `main.cpp`, `periodicNode_d.cu` | **done** (§9) |
| **2** 陰解法結合 | 種ブロック point-implicit (LU)、案C 予測子経由の**反応熱の陰的注入**、$(4,4)$ 反応熱 Jacobian (precond 版含む)、falloff (Troe)、Strang 分離 (非定常 RK) | `speciesTransport_d.cu`, `chemistry_d.*`, `timeIntegration_d.cu`, `main.cpp`, `case/46`, `case/35` | **done** (§9) |
| **3** 燃焼器 RANS | Burrows–Kurkov (case/47): SST 壁関数 + No-TCI で自己着火・火炎を再現 (定量は混合不足)。PaSR 実装済 (未較正)。加熱器 (低マッハ) は未 | `case/47`, `chemistry_d` (tci) | **in_progress** (§9) |
| **4** 応用 | 加熱器 → ノズル → 試験部 end-to-end、設計チェーン確認メニューへ | `design/`, `methods/design/` | todo |

## 6. 検証

- **単体 / ビルド**: `tools/test_chemistry.cpp` (0-D 定積着火遅れ・定圧平衡到達 vs Cantera 同一 YAML、Jacobian 有限差分照合 $<10^{-6}$ 相対)。
- **検証ケース**:
  1. Q1D ノズル再結合: 新設 `case/46.nozzle_h2o2_kinetics` (H₂/O₂ 平衡入口、軸対称 node)。参照は Cantera PFR (面積則)。
     frozen / equilibrium / finite-rate の $T$, $Y_{OH}$, 凍結位置。
  2. Burrows–Kurkov (M2.44 vitiated air + H₂ 壁噴射、NASA GRC 検証アーカイブ)。
  3. 回帰: `case/28.cutler_coaxial_jet` (非反応多成分) と `case/44` (TP 単一擬似種) がビット不変 (`chemistry.enabled: 0`)。
- **判定基準**: 0-D 着火遅れ $\pm2$ %、平衡温度 $\pm$ 数 K; Q1D ノズルの $T$ 差 $<1$ %、$Y_{OH}$ 凍結値 $\pm10$ %;
  反応 ON で `cfl_pseudo` 上限が非反応と同等 (剛性が陰化で吸収されていること); 全残差 `check_convergence.py` PASS。

## 7. 影響範囲

- 新規: `cuda_forge/chemistry_d.*`, `tools/test_chemistry.cpp`, `tools/mechanisms/`, `tools/cea_thermo_to_species_db.py`。
- 変更: `speciesTransport_d.cu` (ブロック Jacobian), `thermo_d.*` (`a7_abs`), `species_eos_coupling_d.cuh`, `solverConfig.*`, `main.cpp`,
  `boundaryCond_d.cu` (`Yb`)。既定 (`chemistry.enabled: 0`) では全経路ビット不変。
- ドキュメント: `methods/chemistry.md`, `methods/index.md`, `methods/thermophysics.md` (datum 節に反応時の注記), `procedures/solver-settings.md` (config キー)。

## 8. 完了条件

- [x] `methods/chemistry.md` 起草 (理論・実装方針)
- [x] Phase 1 実装・検証 (0-D 着火)
- [x] Phase 2 実装・検証 (ノズル化学非平衡 case/46, falloff, Strang)
- [ ] Phase 3 実装・検証 (§6)
- [ ] `status: done` + §9 変更ログ
- [ ] `plans/active/` → `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-04` — 初稿。Phase 0 完了:
  - **CEA スクリーニング** (`notes/investigations/scripts/cea_kinetics_screen/run_screen.py`, 結果 `screen_table.md`, $p_c$ 11.39 bar,
    出口 $A/A^*$ 15.07 (M≈4.2) と 53.2 (M≈6)): 燃焼室温度 $T_c$ ごとの平衡 vs 凍結 (nfz=1 燃焼室凍結 / nfz=2 スロート凍結) の出口差。
    | $T_c$ [K] | 組成 | $X_{OH,c}$ | $X_{NO,c}$ | ΔM (M4.2, eq−fz1) | ΔT (M4.2) [K] | ΔM (M6) |
    |---|---|---|---|---|---|---|
    | 1161 (case/44 va3) | va3 | 0 | 1.3e-4 | 0.000 | 0.1 | 0 (凝縮を除く) |
    | 1595 | va3 | 9e-5 | 1.8e-3 | −0.003 (0.07 %) | +2 | −0.004 |
    | 1979 | va3 | 8.9e-4 | 6.7e-3 | −0.013 (0.3 %) | +10 | −0.019 |
    | 2422 | va3 | 5.1e-3 | 1.8e-2 | −0.043 (1.1 %) | +41 | −0.067 |
    | 2049 (H₂ 加熱器, H₂O 26 %) | B40 | 2.6e-3 | 7.0e-3 | −0.017 (0.4 %) | +14 | −0.025 |
    → **$T_c\le1600$ K では化学非平衡は無視可 (0.1 % 未満)、$T_c\approx2000$ K で 0.3–0.4 %、2400 K で 1 %** (設計ゲート 0.5 % Md と同桁)。
    有限速度解は eq と fz2 の間に入る。case/44 (1161 K) と case/45 (1600 K) は現状の frozen 設計で十分、
    本機能が効くのは $T_c\gtrsim2000$ K の加熱器案件。M6 出口の A1161/B10 eq 列は CEA が H₂O(cr) 凝縮を含めるため T が不連続 (forge は凝縮を別モデルで扱う)。
  - **熱力学 DB ツール** `tools/cea_thermo_to_species_db.py`: CEA `thermo.inp` → `species_db.yaml`。内蔵 N2 係数と相対差 0、
    Cantera と $c_p/R$ 4 桁一致 ($h/RT$ の OH/HO₂ 2–5 % 差は生成エンタルピーの出典差)。生成物 `tools/mechanisms/species_db_h2air_cea.yaml`。
  - **機構ファイル** `tools/mechanisms/h2air_jachimowski1988_13sp33r.yaml` / `h2o2_jachimowski1988_9sp20r.yaml` (TP-2791 Table I、熱力学は CEA NASA-9)。
    Cantera 定積着火遅れ (量論 H₂-air, `notes/investigations/scripts/chem_ignition_check.py`):
    | $T_0$ [K] / p [atm] | Jach13 | Jach9 | Cantera h2o2 (GRI 系) |
    |---|---|---|---|
    | 1000 / 1 | 152 µs | 152 | 305 |
    | 1200 / 1 | 32.2 | 32.2 | 44.2 |
    | 1500 / 1 | 10.1 | 10.1 | 12.8 |
    | 2000 / 1 | 3.5 | 3.5 | 3.6 |
    Jach13 = Jach9 (N 化学は着火に無関係) で YAML 化の整合を確認。Jachimowski は低温側で GRI 系より 1.4–2 倍速い (Slack データに合わせた 1988 校正、既知傾向)。
  - **YAML の罠**: 種名 `NO` / `N` は YAML 1.1 (PyYAML) で真偽値に化ける → DB 生成はキーを引用符付きにした。yaml-cpp (forge) 側も同様の扱いを Phase 1 で確認する。
- `2026-09-04` — **Phase 1 完了** (反応ソース項・反応熱・解析 Jacobian・対角 point-implicit)。
  - 実装: `cuda_forge/chemistry_d.cuh` (`ReactionTable`, `chem_source`: Arrhenius・三体・$K_c$・$\dot Q$・$n_s\times n_s$ 解析 Jacobian・$\partial\omega/\partial T$)、
    `chemistry_mech_io.hpp` (Cantera YAML サブセット読込)、`chemistry_d.cu` (device ソース項、`res_roY`/`res_roe`/`src_jac_Y` 加算、`chemQdot`/`chemTau` 出力)、
    `SpeciesThermo::h_datum` (sensible datum の除去分を保持 → $\dot Q=-\sum c_s\dot\omega_s$、$K_c$ は絶対 $H$)、`physProp.chemistry` キー、
    `registerSpecies(n, chemistry)`。既定 (`enabled: 0`) では全経路ビット不変 (構造体変更のため full rebuild)。
  - ホスト検証 (`tools/test_chemistry.cpp`, 参照 `tools/chem_reference_cantera.py`): Jacobian 有限差分照合 $\partial\omega/\partial\rho Y$ 1.8e-8・$\partial\omega/\partial T$ 1.3e-9 (PASS)、
    $\sum\dot\omega_s\sim10^{-13}$ (質量保存)。0-D 定積 BDF1 反応器 (量論 H₂-air 1200 K 1 atm) の着火遅れは刻みを締めると Cantera へ収束
    (29.4 → 31.1 → 32.0 µs; Cantera 32.2)、平衡 T 2945.2 vs 2943.9 K。13 種 33 反応 (1500 K 2 atm) も 4.7 vs 5.0 µs・3065 vs 3064 K。
    絶対 datum と sensible datum で同一の平衡 T (反応熱経路の検証)。
  - **GPU 検証** `case/35.uniform_periodic_box/run_0049_node_h2_ignition` (node 8³ 全面 periodic = 定積反応器, explicit RK3 dt 5e-9 s, 10000 step 94 s):
    着火遅れ 32.0 µs (Cantera 32.2)、平衡 T 2948.4 K (Cantera 2943.9, +0.15 %)、全ノード一様 (T 差 0.004 K, |u| 1e-4 m/s)。
    図 `h2_ignition_vs_cantera.png`。
  - **発見・修正したバグ (化学以前からの欠落)**: node 周期の `periodicNodeGather` が化学種残差 `res_roY{s}` を合算していなかった
    (流れ 5 変数と k/ω のみ)。seam ノードでは化学種の移流・反応が部分体積分 (面 1/2・辺 1/4・角 1/8) しか効かず、エネルギーだけ合算されて
    不整合 (run_0049 初回: seam が着火せず圧力波で場が乱れた)。`res_roY{s}`・凝縮モーメント残差の gather と `roY{s}` の root→member ミラーを追加。
    体積ソースは `volumePartial_d` を使う (bodyForce/ransSource と同じ規約)。**既往の node 周期×多成分 run (case/28 等は非周期なので影響なし) を要確認**。
  - 残: Phase 2 (種ブロック point-implicit、案C 結合、falloff、Strang)、`chemQdot` の陰解法 $(5,5)$ Jacobian、FP32 化。
- `2026-09-04` — **Phase 2a/2b 完了** (ブランチ `feature/chemistry-finite-rate`, worktree `forge-chem`)。
  - **種ブロック point-implicit** (`jacobianMode: 2`): chemistry が温度結合込み $J^{\rm tot}_{sk}=J_{sk}+\partial\omega_s/\partial T\cdot\partial T/\partial\rho Y_k$ を
    `chem_jac` に書き、`species_dplur_block_sweep_d` がセル毎 $n_s\times n_s$ LU (double) で δ(ρY) を解く (coupling 1 の解・coupling 2 の予測子の両方)。
    $(4,4)$ 対角に $\max(0,-\partial\dot Q/\partial\rho e)$。
  - **発散と根治**: `case/46.nozzle_h2o2_kinetics` (燃焼室 CEA 平衡 2788 K・11.39 bar, va3 ノズル, 15 種) で反応熱を陽的に入れる構成
    (run_0002–0007: jac 1/2, coupling 0/1/2, cfl_pseudo 2/0.2) は**全て数 step で NaN**。スロート付近の再結合熱 $\dot Q\approx7.5\times10^{11}$ W/m³
    (Cantera 評価, $\tau_c$ 4e-8 s) が擬似 Δt 1e-6 s で内部エネルギーの 30 %/step。案C 予測子 δ(ρY)* で線形化した
    $V\sum_k(\partial\dot Q/\partial\rho Y_k)\delta(\rho Y_k)^*$ を `res_roe` に注入 (`chemistry_heat_inject_d`) すると
    **run_0008 は cfl_pseudo 2 (frozen と同じ) で安定・12000 step 完走**、`rms_roe` 4.3 桁低下 (frozen 基準と同じプラトー水準)。
    → 定常陰解法の反応流は `speciesImplicitCoupling: 2` + `jacobianMode: 2` 必須 (config に警告を追加)。
  - **検証** (run_0008 vs Q1D Cantera 参照 `q1d_kinetics_ref.py`): 軸線 $Y_{OH}$ が forge スロート状態 (M=1.3) 始点の Q1D と全域一致
    (出口 6.39e-5 vs 6.44e-5)、NO は Zeldovich が遅く凍結 (3.07e-2, Q1D 同値; CEA 平衡スロートの 2.3e-2 は非現実)。
    frozen (run_0001, 燃焼室組成凍結) 比で出口 T +58 K, M −0.045 (−1.1 %); CEA 平衡スロート始点 Q1D の fr−fz は +40 K, −0.044。
    T/M の絶対値は 1D 面平均と軸値の差 (出口 M 3.85 vs 3.91) があり比較対象にしない。出口量は 4000→12000 step で T 0.5 K・$Y_{OH}$ 0.3 % の変化 (準定常)。
    `check_convergence`: run_0001/run_0008 とも STALLED プラトー (NOT CONVERGED; run_0091 系の既知 warm 床)。
  - **falloff (Lindemann/Troe)**: `chem_ln_kf_falloff` + 温度/[M] 微分は滑らかなスカラ関数の中心差分 (相対 1e-4)。
    Cantera h2o2.yaml (GRI H₂/O₂ 29 反応, Troe 2 本) を CEA NASA-9 熱力学に付け替えた `tools/mechanisms/h2o2_gri30_nasa9_10sp29r.yaml` で
    Jacobian FD 2.1e-8 PASS、着火遅れ 43.6 vs Cantera 44.0 µs (1200 K 1 atm)。Burke 2012 機構は同形式で読める。
  - 残: Strang 分離 (非定常)、`_precond_d` (低マッハ前処理) への $(4,4)$ 配線、`inlet_Pressure_dir` 組成、Phase 3 (燃焼器 RANS・PaSR)、
    反応 ON での `cfl_pseudo` 上限スイープ (frozen 同等の 2 は確認、4 以上は未)。
- `2026-09-04` — **Phase 2c (Strang 分離) + CFL 上限 + precond 配線**。
  - `chemistry.strang: 1`: 非定常 RK の前後 dt/2 でセル内 ODE を backward Euler sub-cycle ($h\le0.5\tau_c$, ≤64 回, 2 Newton, $\Delta\rho e=-\sum c_s\Delta\rho Y_s$)。
    有効時はソース項を RK に入れず、前半後に N/M 始点を取り直す (`advanceExplicitRK`)。
    case/35 `run_0052` (dt 5e-9) / `run_0053` (dt 5e-8, 10 倍): いずれも着火遅れ 32.0 µs・平衡 2945 K (Cantera 32.2 µs・2944 K) で dt に依らず一致、場は厳密一様。
    ソース項経路 (run_0049: 2948.4 K) より平衡温度が正確。**コストは重い** (run_0052 913 s vs run_0049 94 s; 729 ノードで GPU 未飽和 + local memory の 16×16 double 行列) → 要最適化 (FP32 化・レジスタ化)。
  - `cfl_pseudo` 上限 (case/46 run_0009/0010, 反応 ON 定常陰解): 4 は安定 (残差は 2 と同水準)、8 はリミットサイクル (残差 2 桁悪化) → 実用上限 ≈4。
  - `implicit_defect_correction_block_precond_d` (lowMachPrecond) にも反応熱の $(4,4)$ 対角 (`src_jac_roe`) を配線 (加熱器の低マッハ用。未検証)。
  - 残: Phase 3 (燃焼加熱器 RANS: 低マッハ precond + SST + No-TCI → PaSR、Burrows–Kurkov)、`inlet_Pressure_dir` 組成、Strang/ソースカーネルの性能最適化。
- `2026-09-04` — **Phase 3 一次結果 (Burrows–Kurkov, case/47)**。
  - 形状・条件は NASA TM X-2828 から (`papers/combustion/Burrows_Kurkov_1973_NASA-TM-X-2828.pdf`): 主流 M2.44/1270 K/1 atm (Y O₂ .258/N₂ .486/H₂O .256),
    H₂ スロット 0.4 cm (リップ 0.076 cm) M1/254 K, 燃焼器 35.6 cm (上壁 9.38→10.5 cm 直線拡大)。gmsh 平面構造 68809 ノード (品質 PASS)。
    実験の出口プロファイル (Fig.6 組成, Fig.13 全温) は図の目視読み取り (±0.03)。
  - **起動の罠 (node モード)**: ①`wall_isothermal` の壁温キーは `Ts` (`T` は無視され 0 K)。②IC の壁ノードは u=0・T_wall で作る (流速込み roe は T=0→NaN)。
    ③**入口と壁を共有する角ノード**で `inlet_uniformVelocity` と壁 Dirichlet が衝突し P 数十 bar → NaN (1 次では延命、2 次で発散)。
    入口 BL プロファイル (`inletProfile` + Crocco–Busemann ρ) でも解消せず → **未解決の制約**として記録 (別 plan 候補)。
    回避: 上流ダクト/スロット壁と燃焼器上壁を slip、no-slip 等温壁 (壁関数) はステップ面と燃焼器下壁のみ。段階起動 (1 次 cfl 0.5 → 2 次 cfl 2)。
  - **結果** (`run_0007` 混合 → `run_0008` 反応, jacobianMode 2 + coupling 2, cfl_pseudo 1, 20000 step, NaN なし):
    自己着火と壁噴流火炎を再現。出口 X_H₂O ピーク 0.51 (実験 0.50) だが位置 y 1.48 cm (実験 2.0 cm)、全温ピーク 1.10 (実験 1.18)、
    着火点 x 5–9 cm で上流へドリフト中 (実験 18–25 cm)。混合基準でも H₂ 層が実験より薄く、主流 M 2.55 (上壁 slip)。
    → 差の主因は**混合不足** (上流壁 slip で初期 BL/せん断層厚ゼロ、入口乱流量が仮置き) と着火の早さ (No-TCI)。
    残差は振動プラトー (`check_convergence` NOT CONVERGED、リップ後流)。
  - 次: 入口角ノードの解消 (solver 側) → 全壁 no-slip + 入口 BL、入口 k/ω の感度、PaSR (`tci: 1`) の C_mix 較正、cell モードとの比較。
- `2026-09-04` — **Phase 3 更新: 入口∩壁角の根治後 (plan `boundary-node-inlet-corner-wall`)**。全壁 no-slip 等温 + 入口 BL プロファイルで
  run_0014–0017 が安定。反応 run_0017: 出口組成が実験にほぼ重なり、X_H₂O ピーク 0.51 @ 1.78 cm (実験 0.50 @ 2.0)、全温ピーク 1.06 (実験 1.18)、
  着火点 x≈9 cm (定常化; 実験 18–25 cm)。残差はプラトー (NOT CONVERGED)。残る差: 主流 M 2.57 (実験 2.2, 壁関数 BL の blockage 不足)、
  全温ピーク不足と早い着火 (No-TCI/壁関数/入口 k-ω 仮置き)。次: PaSR `tci: 1` の C_mix 感度、入口 k/ω・δ 感度、y⁺ 確認。


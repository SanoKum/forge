# SU2-NEMO contact/interface thermo 取り扱い調査 (forge mixed-order face-state 比較)

<!-- ファイル名規約: <area>-<short-slug>.md -->

## メタ

- **area**: `convection / thermophysics`
- **status**: `draft` <!-- 調査専用。コード変更は未着手 -->
- **related_docs**:
  - `methods/thermophysics/` (多成分 TP gas)
  - `methods/convection/` (SLAU/Roe 再構成)
  - `methods/time_integration/` (block-DPLUR / general-EOS Jacobian)
- **related_plans**:
  - `thermophysics-multicomponent-tpgas.md`
  - `thermophysics-species-implicit-coupling.md`
  - `time_integration-general-eos-jacobian.md`
  - `architecture-axisym-axis-singularity.md`
- **created**: `2026-06-20`
- **owner**: `CFD Dev`

## 0. 調査対象リビジョン (記録)

forge ローカルの `.external/su2/` は**バイナリ配布のみ** (`SU2-v8.5.0-linux64-omp.zip` + `bin/`)
で**ソースは存在しない**。そのため同一バージョンのソースを read-only で取得して調査した。

- 取得先: `https://github.com/su2code/SU2.git` tag `v8.5.0`
- 配置: `.external/su2-src/` (`.gitignore` 済、commit 対象外)
- `git rev-parse HEAD`: `12eb826f049ef7f67df974dfcb44cf36ee07c0f8`
- `git describe --tags`: `v8.5.0`
- ローカルバイナリ (`SU2-v8.5.0-linux64-omp`) と同一タグ。

> 注: 設問の `git rev-parse HEAD` 等をローカル `.external/su2` で実行すると forge リポジトリの
> HEAD (`11162c2`, branch `feature/axisym-near-axis-cfl`) が返る (su2 ディレクトリは独立 git でないため)。
> SU2 側の真の対象 commit は上記 `12eb826` (v8.5.0)。

調査は本ディレクトリの実コードを根拠とする。以下、**[確認]=コードで確認** / **[推定]=構造からの推論** を区別する。

---

## 1. 調査サマリー (結論先出し)

**最重要の発見**: SU2-NEMO は forge の mixed-order face-state を**持たない**。
species と flow を**同一次数 (同一 MUSCL ループ・単一 limiter) で co-reconstruction** する。
ただし forge と別種の非整合がある — **face の `T,P,ρ,h` を EOS から再回復しない** (raw な
linear extrapolation をそのまま使う)。

| 問 | 回答 (要旨) |
| --- | --- |
| 1. forge と同じ mixed-order face-state があるか | **無い**。species も flow と同じ 2 次で再構成。forge 固有の「2次 P,ρ + 1次 Y_cell で T_f」構造は存在しない |
| 2. species と flow は同一 face state から flux を作るか | **Yes**。`Primitive_i/j` 全成分 (species 含む) を numerics に渡し、同一 face state から flux 構成 |
| 3. 多成分 contact の偽圧力振動を特別処理するか | **No**。double-flux / pressure-equilibrium 補正は無い。完全保存形式 → Abgrall/Karni 型の偽圧振動に**晒される** |
| 4. double-flux / quasi-conservative を使うか | **No** (tree 全 grep でヒット無し) |
| 5. 高 CFL contact limit-cycle を抑える仕組み | 専用機構は無い。間接的に (a) species/flow 同次 co-reconstruction (b) `CheckNonPhys`→該当 face を 1 次へ revert (c) `T∈[50,8e4]` clamp (d) full coupled Jacobian。CFL は residual-based adaptation |
| 6. implicit Jacobian の coupling 範囲 | **full `nVar×nVar` block** (species⇔flow を `dPdU` で連成)。近似 (frozen-reconstruction, limiter 微分なし)。chemistry/vib-relax/viscous の source Jacobian あり |
| 7. forge へ取り入れる価値 | (i) **species を flow と同次で再構成** (最優先), (ii) full-coupled `∂P/∂ρ_s` Jacobian (既に general-EOS plan で着手済), (iii) 非物理 face の 1 次 revert |
| 8. それでも double-flux が別途必要か | **可能性高い**。SU2-NEMO 自身も double-flux を持たず偽圧振動に晒される。完全保存形式である限り、γ が大きく変わる contact では同次再構成だけでは根治しない |

**forge の問題への含意**: forge の limit-cycle の核心は「2次 `P_f,ρ_f` と 1次 `Y_cell` の**次数不整合**が
`T_f = P_f/(ρ_f R(Y_cell))` を通じて face thermo を歪めること」。SU2-NEMO はこの**次数不整合を作らない**
(Y も 2 次) 点で forge より素直。ただし SU2 も contact 偽圧振動の根治機構 (double-flux) は持たないので、
「SU2 を真似て species を同次再構成する」のは**必要だが十分でない**可能性がある。

---

## 2. コードパス・関数・行番号一覧

ルート: `.external/su2-src/SU2_CFD/`。行番号は近傍値。

### 変数・primitive (CNEMOEulerVariable / CNEMOGas)
- 保存ベクトル `U=[ρ_1..ρ_Ns, ρu,ρv(,ρw), ρE, ρEve]`, `nVar=nSpecies+nDim+2` — `src/solvers/CNEMOEulerSolver.cpp:100,104`
- `Cons2PrimVar` (T,Tve,P,a,h 回復の本体) — `src/variables/CNEMOEulerVariable.cpp:168-277`
- primitive index 定義 (`RHOS/T/TVE/VEL/P/RHO/H/A/RHOCVTR/RHOCVVE_INDEX`) — `CNEMOEulerVariable.cpp:56-67`
- `ComputeTemperatures` (T 閉形式 + Tve Newton/二分, clamp `[50,8e4]`) — `src/fluid/CSU2TCLib.cpp:2109-2196`
- `ComputePressure` (Dalton 則) — `src/fluid/CNEMOGas.cpp:103-114`
- `ComputeSoundSpeed` (frozen `a²=(1+Ru·conc/(ρCv_tr))P/ρ`) — `CNEMOGas.cpp:83-101`
- `ComputedPdU` (`∂P/∂ρ_s` 含む) — `CNEMOGas.cpp:156-206` (species 項 `:203-205`)

### MUSCL 再構成 (Upwind_Residual)
- 再構成ループ `iVar=0..nPrimVarGrad` (全 primitive・species 含む, 単一 `lim_ij`) — `CNEMOEulerSolver.cpp:551-576`
- `nPrimVarGrad=nSpecies+nDim+8` — `CNEMOEulerSolver.cpp:107`
- `ComputeConsistentExtrapolation` (**eves と gamma のみ**再計算, T/P/h は再回復せず) — `CNEMOEulerSolver.cpp:654-678`
  - 自認コメント「`this doesnt compute Cvves/dPdU,etc.yet`」 — `CNEMOEulerSolver.cpp:654`
- `RecomputeConservativeVector` (Roe/MSW のみ, 再構成 primitive→U_f) — `CNEMOEulerSolver.cpp:683-718`
- `CheckNonPhys` (face 非物理→1 次 revert) — `CNEMOEulerSolver.cpp:720-750`, 適用 `:604-611`
- 再構成 state を numerics へ — `CNEMOEulerSolver.cpp:604-606`

### 数値流束
- Roe `CUpwRoe_NEMO::ComputeResidual` — `src/numerics/NEMO/convection/roe.cpp:92`
  - Roe sound speed (frozen, 再計算) `:132`; entropy fix `:146-163`; residual+Jacobian ループ `:170-195`
- AUSM 基底 `CUpwAUSM_SLAU_Base_NEMO::ComputeResidual` — `convection/ausm_slau.cpp:374`
  - 派生: `CUpwAUSM_NEMO:446`, `CUpwAUSMPLUSUP2_NEMO:485`, `CUpwAUSMPLUSM_NEMO:560`
  - species flux = `ρ_s × split mass flux` — `ausm_slau.cpp:421-436`
- MSW (FVS) `CUpwMSW_NEMO` — `convection/msw.cpp:110`
- Lax (centered) `CCentLax_NEMO` — `convection/lax.cpp:76`
- 物理流束 `GetInviscidProjFlux` — `src/numerics/NEMO/CNEMONumerics.cpp:88-135`
- factory dispatch (SLAU は NEMO 未対応, default error) — `src/drivers/CDriver.cpp:1894-1959`

### implicit Jacobian / source / dt
- `GetInviscidProjJac` (full block, species⇔flow `dPdU` 連成) — `CNEMONumerics.cpp:153-216` (species 結合 `:186-189`)
- chemistry source Jac `ComputeChemistry`/`ComputeNetProductionRates` — `src/numerics/NEMO/NEMO_sources.cpp:69-110`
- vib-relax source Jac `ComputeVibRelaxation` — `NEMO_sources.cpp:112-168` (τ 微分は frozen, コメント `CNEMOEulerSolver.cpp:823`)
- viscous Jac `GetViscousProjJacs` — `CNEMONumerics.cpp:309-352`, 呼出 `NEMO_diffusion.cpp:168-184`
- local dt `dt=CFL·Vol/Σ(|u·n|+c)` — `include/solvers/CFVMFlowSolverBase.hpp:506-514`, 固有値 `:421-437`
- CFL adaptation (residual-based) `AdaptCFLNumber` — `src/solvers/CSolver.cpp:1745-1934`

### axisymmetric
- `ComputeAxisymmetric` (1/r source, `yinv` 軸ガード) — `NEMO_sources.cpp:170-312` (軸ガード `:194-196`)
- 軸対称 source Jacobian (inviscid のみ) — `NEMO_sources.cpp:211-268`, 対角加算 `CNEMOEulerSolver.cpp:882-883`
- 軸=`SYMMETRY_PLANE`, `BC_Sym_Plane` — `CNEMOEulerSolver.cpp:1419+`; dt で axis/wall edge 除外 `CFVMFlowSolverBase.hpp:454-455`

---

## 3. face-state 計算フロー (実関数名付き)

SU2-NEMO (MUSCL on, Roe/AUSM):

```
cell U  ──Cons2PrimVar(CNEMOEulerVariable.cpp:168)──▶  cell primitive V
                                                        [ρ_s, T, Tve, u, P, ρ, h, a, ...]
   │
   ├─ SetPrimitive_Gradient / Limiter (base CSolver)        ← 全 primitive の勾配・limiter
   ▼
Upwind_Residual (CNEMOEulerSolver.cpp:530)
   │  iVar = 0 .. nPrimVarGrad-1   ← species と flow を区別せず全部
   │  Primitive_{i,j}[iVar] = V_{i,j}[iVar] ± 0.5·lim_ij·ProjGrad   (:573-576)
   │      → ρ_s,f / ρ_f / P_f / u_f / T_f / Tve_f / h_f / a_f すべて 2 次 (単一 limiter)
   │
   ├─ ComputeConsistentExtrapolation (:594)  →  eves, gamma のみ再計算
   │       ※ T_f, P_f, ρ_f, h_f は再構成値のまま (EOS 再回復なし)
   ├─ RecomputeConservativeVector (:599, Roe/MSW)  →  U_f = f(reconstructed prim)
   ├─ CheckNonPhys (:581) → 非物理なら face 全体を cell 値へ revert (1 次)
   ▼
numerics->SetPrimitive(Primitive_i, Primitive_j)  (:605)
numerics->ComputeResidual (roe.cpp:92 / ausm_slau.cpp:374)
   →  flux + Jacobian (同一 face state)
```

**forge との決定的差**:
- forge: `T_f = P_f^{(2)} / (ρ_f^{(2)} · R_mix(Y_cell^{(1)}))` ← **2次 P,ρ と 1次 Y の次数混在**
- SU2 : `T_f` は**直接 2 次再構成**された値 (Y_f も 2 次)。次数混在なし。
  ただし `T_f` は `P_f,ρ_f` から独立に再構成されるため、`P_f ≠ ρ_f R(Y_f) T_f` の**別種非整合**は残る。
  これを `CheckNonPhys`→1 次 revert と `T` clamp で抑える設計。

### 再構成表 (SU2-NEMO) [確認]

| 量 | セル値 | 1次 | 2次再構成 | limiter付き | faceで再計算 |
| --- | :-: | :-: | :-: | :-: | :-: |
| species density ρ_s | — | (revert時) | **✔** | ✔ (共有) | — |
| density ρ | — | (revert時) | **✔** | ✔ | — |
| pressure P | — | (revert時) | **✔** | ✔ | — |
| velocity u | — | (revert時) | **✔** | ✔ | — |
| temperature T | — | (revert時) | **✔** | ✔ | **✘ (EOS再回復せず)** |
| vibrational T_ve | — | (revert時) | **✔** | ✔ | ✘ |
| enthalpy h | — | (revert時) | **✔** | ✔ | ✘ |
| sound speed c | — | (revert時) | **✔** | ✔ | ✘ |
| eves (vib energies) | — | — | — | — | **✔ (再構成 ρ_s,Tve から)** |
| gamma | — | — | — | — | **✔** |

→ **pattern 判定**: 提示の A〜D では純粋な A/B/C いずれでもなく **「A 相当だが EOS 再回復が stub」(=D)**。
species と flow を同次 co-reconstruction (B 否定) するが、`T_f,h_f,c_f` を `P_f,ρ_f,Y_f` から作り直さない
(A の最終節を満たさない)。primitive を再構成 (C の conservative 再構成ではない)。

---

## 4. implicit Jacobian 構造

[確認] **full `nVar×nVar` block** (`nVar=nSpecies+nDim+2`)。

- inviscid: `GetInviscidProjJac` (`CNEMONumerics.cpp:153-216`) が species 行/列を含む全結合を構築。
  熱力学連成は `val_dPdU` (=`∂P/∂U`) 経由:
  - species→momentum `:186`: `+= dPdU[iSpecies]·n - proj_vel·u`
  - species→energy `:188`: `+= (dPdU[iSpecies]-H)·proj_vel`
  - species→vib `:189`: `+= -proj_vel·ρEve/ρ`
  - `∂P/∂ρ_s` は `CNEMOGas::ComputedPdU` `:203-205` で計算。
- `dTdU`,`dTvedU` は inviscid Jac には入らず、**source Jacobian** (chemistry/vib-relax) でのみ消費。
- **exact か近似か**: LHS は RHS と**同一の再構成 face state**で評価 (`SetPrimitive(Primitive_i/j)`)
  → 古典的 1 次 defect-correction では**ない**。ただし
  - limiter / 勾配の neighbor 微分 (`∂Primitive_i/∂U_neighbor`) は LHS に**入らない** (frozen-reconstruction)
  - Roe `|Λ|=P|Λ|P⁻¹` は U に対し定数扱い
  → **「近似 (frozen-reconstruction, limiter 微分なし) full-coupled Jacobian」**。
- source Jacobian:
  - chemistry `∂ω̇/∂U` (`dTdU`/`dTvedU` 使用) — `NEMO_sources.cpp:94-105`, 対角へ `SubtractBlock2Diag`
  - vib-relax (τ 微分は frozen) — `NEMO_sources.cpp:146-149`
  - viscous `GetViscousProjJacs` (NS) — `CNEMONumerics.cpp:309-352`
- axisymmetric source Jacobian: inviscid 1/r 項のみ LHS 対角へ (`NEMO_sources.cpp:211-268`)。
  viscous 1/r 項は residual のみ (LHS なし)。

**forge との対比**: forge の general-EOS Jacobian plan
(`time_integration-general-eos-jacobian.md`) が目指す「真の `∂F/∂Q` 整合 (`∂P/∂ρ_s` 含む固有系)」は
SU2-NEMO の `GetInviscidProjJac`+`dPdU` 連成と**同じ方向**。SU2 も exact ではなく近似である点は共通。

---

## 5. forge との比較表

| 項目 | forge 現状 | SU2-NEMO v8.5.0 |
| --- | --- | --- |
| flow reconstruction | 2次 | 2次 (primitive) |
| species reconstruction | **1次 (cell値)** | **2次** (flow と同一ループ・同一 limiter) |
| face thermo 用 species | **cell 値** | **再構成値** (face Y_f) |
| `T_f` の作り方 | `P_f/(ρ_f R(Y_cell))` (EOS, 次数混在) | **直接 2 次再構成** (EOS 再回復なし) |
| enthalpy | TP mixture | 直接再構成 (`h_f`), eves/γ のみ再計算 |
| species-flow Jacobian | segregated / partial (案B は無効と確定) | **full coupled** (`dPdU` で `∂P/∂ρ_s` 連成) |
| pressure-equilibrium correction | なし | **なし** |
| double-flux | なし | **なし** |
| contact 近傍低次化 | なし | 間接的 (`CheckNonPhys`→非物理 face を 1 次 revert) |
| positivity 処理 | **ρY_s≥0 クランプ + ΣY_s=1 再正規化あり** (`species_renormalize_d`)、realizability floor、DPLUR/block 更新で `max(·,0)`、T clamp `[50,6000]` | ρ_s floor `1e-20`; T/Tve clamp `[50,8e4]`; P floor `1e-20`; **ΣY=1 正規化はなし** |
| CFL 制御 | cfl_pseudo 固定 (実効) | residual-based adaptation (`CFL_AdaptParam`) |
| axis source Jacobian | scalar/ block 対角に hoop 源を陰化 | inviscid 1/r 源を LHS 対角へ (viscous は explicit) |
| 近軸 extra spectral radius | `λ_axis=β(|u_r|+c)A_planar` (本ブランチで追加) | **なし** (face 固有値 Σ(|u·n|+c) のみ; axis edge は除外) |
| sound speed (dt/flux) | (TP frozen) | frozen `a²=(1+Ru·conc/ρCv_tr)P/ρ` |
| species mass flux | `F_{ρY_s}=massflux·Y_{s,up}` (1次風上, segregated scalar; massflux=flow の SLAU/AUSM mdot)。**ΣF=mdot で総質量と整合** [確認] | `F_{ρ_s}=ρ_s·(split mass flux)`=`Y_s,up·F_ρ` (AUSM); Roe/MSW は固有系。ΣF_{ρ_s}=F_ρ [確認] |
| species 再構成の所在 | **segregated scalar 移流 (1次風上のみ, MUSCL なし)** + face thermo で cell 組成使用 → 2 重に 1 次 | flow と同一 MUSCL primitive ループ (2次) |
| 既存の contact 低次化機構 | **あり (診断, 既定 off)**: `FORGE_CONTACT_1ST` (組成センサで flow を hard 1次化), `g_contactBlend` (連続 2→1次ブレンド)。species 流束は不変 | `CheckNonPhys`→非物理 face を 1 次 revert (常時) |

---

## 6. 取り入れるべき設計の優先順位

1. **[最優先] species を flow と同次で再構成する** — forge の mixed-order (1次 Y) を廃し、
   `Y_s`(または `ρ_s`) を `P,ρ,u` と同じ 2 次・同一 limiter で face へ。これで
   `T_f` を `R_mix(Y_f)` で評価しても**次数混在が消える**。SU2 の `Upwind_Residual` ループ
   (`CNEMOEulerSolver.cpp:551-576`) が直接の参照実装。
   - **forge では 2 箇所の改修が要る** (調査で確定): (i) species 移流を 1 次風上から MUSCL 2 次へ
     (`scalarTransport_d.cu:25` の `phi_upwind=cell値` を再構成値に)、(ii) face thermo の
     `Rmix_cell[ic0/ic1]` を**再構成した face 組成 `Y_face`** から作る (`convectiveFlux_d.cu:436-439`)。
     現状は移流・face thermo の双方が cell 組成 (2 重に 1 次) なので、両方を揃えないと次数整合しない。
   - segregated scalar 構成 (`scalarTransport`) を保ったまま MUSCL 化は可能だが、face 組成を
     convectiveFlux 側 thermo に渡す配線が新規に必要 (現状 face 組成は thermo ローカルに無い)。
2. **[高] 非物理 face の 1 次 revert** (`CheckNonPhys` 相当) — 再構成 face で `ρ_s<0 / P<0 /
   T∉[Tmin,Tmax] / a<0` を検出したら当該 face を cell 値へ落とす。contact での過大再構成を局所抑制。
3. **[中] full-coupled `∂P/∂ρ_s` Jacobian** — forge は既に
   `time_integration-general-eos-jacobian.md` で着手。SU2 の `GetInviscidProjJac`+`ComputedPdU`
   を参照に `dPdU` 連成を確認・補完。
4. **[低/検討] residual-based CFL adaptation** — limit-cycle 域で CFL を自動で絞る安全弁。
5. **[要別途検討] double-flux / pressure-equilibrium 補正** — §7 参照。SU2 にも無いため
   「SU2 を真似る」だけでは得られない。forge 固有の根治レバーとして別 plan 化を検討。

---

## 7. SU2 を参考にしても double-flux が別途必要か

**必要になる可能性が高い [推定、ただしコード根拠あり]**。

- SU2-NEMO は完全保存形式で、tree 全 grep
  (`double.?flux|quasi.?conservative|pressure.?equilibrium|Abgrall|Karni|frozen.?gamma|...`)
  でも double-flux / pressure-equilibrium 補正は**皆無**。AUSM+up2/+M の `M_p` 低マッハ項は
  all-speed 散逸であって contact 偽圧振動の cure ではない。
- したがって SU2-NEMO 自身が、γ が大きく変わる多成分 contact では Abgrall/Karni 型の
  **偽圧力振動に晒される**。SU2 の安定性は「species/flow 同次再構成 + 非物理 1 次 revert +
  full-coupled 近似 Jacobian + residual-based CFL」の総合で得ており、**振動の根治ではない**。
- forge の症状が「高 CFL での limit-cycle」である以上、まず §6-1 (同次再構成) で
  forge 固有の次数混在を除去し、それでも残る場合は **double-flux / quasi-conservative**
  (例: face で γ_eff を左右セルから凍結し energy flux を L/R 別評価) を forge 独自に導入する
  必要がある。これは SU2 移植では賄えない。

---

## 8. 不明点・未確認事項

### 解決済 (forge 側 `convectiveFlux_d.cu` / `scalarTransport_d.cu` / `speciesTransport_d.cu` を調査)
- **species mass flux 構成 [確認]**: `Y_s,up·F_ρ` 型で総質量と整合。`scalarTransport_d.cu:23-26`
  が `mdot=massflux[ip]; phi_up=(mdot≥0)?Y[ic0]:Y[ic1]; flux=mdot·phi_up`。`massflux[ip]`
  は `convectiveFlux_d.cu:557` で flow の SLAU/AUSM mdot を一度だけ確定 → species は同一 mdot を
  共有する **segregated scalar 移流**。ΣY=1 再正規化が効くため Σ_s F = mdot で連続式と整合。
  独立 Riemann ではない。SU2-NEMO AUSM と同型。
- **species は 2 重に 1 次 [確認]**: (i) 移流自体が 1 次風上 (MUSCL なし, `scalarTransport_d.cu:25`)、
  (ii) face thermo も cell 組成 `Rmix_cell[ic0/ic1]` を使用 (`convectiveFlux_d.cu:420,436-439`)。
  T_face = `P_L/(ro_L·Rmix_cell[ic0])` がまさに forge の mixed-order 式。
- **positivity [確認]**: forge は SU2 より手厚い。`species_renormalize_d`
  (`speciesTransport_d.cu:98-118`) が ρY_s≥0 クランプ後 **ΣρY_s=ρ へ再正規化 (ΣY=1)**;
  scalarTransport realizability floor (`:169`); DPLUR 更新 `max(roYN+dq,0)` (`:307`);
  block-triangular `max(v,0)` (`:392`); face thermo の局所 Y クランプ+再正規化 (`:429-433`);
  T clamp `[50,6000]` (`DEPVAR_TMIN/TMAX`, `dependentVariables_d.cu:8-9`)。SU2 は ΣY=1 再正規化なし。
- **forge は既に contact 低次化の実験機構を保有 [確認]**: env ゲート (既定 off, ビット不変):
  `FORGE_CONTACT_1ST` (組成センサ `s_Y=max|ΔY|/(ΣY+ε)>thresh` で flow を hard 1 次化,
  `convectiveFlux_d.cu:336`)、`g_contactBlend` (連続 2→1 次ブレンド, `:391-398`)、
  `FORGE_CONTACT_LOG` (L/R 状態 printf, `:491-501`)。**いずれも species 流束は不変** = flow のみ低次化。
  SU2 の `CheckNonPhys`→1 次 revert と方向性は同じだが、forge は「組成センサ起動・flow のみ」で
  常時適用ではない点が異なる。

### 未確認 (残る)
- SU2 の `ComputeConsistentExtrapolation` が `T_f,P_f` を再回復しないことの**数値的影響**は
  コード読みだけでは定量化不可。SU2 実行による確認は未実施。
- forge の `λ_axis` 補正に対応する機構は SU2 に無い [確認] が、SU2 が近軸 CFL でどの程度
  安定かは未測定。
- forge の既存 `g_contactBlend` を limit-cycle ケースに当てた効果 (振幅減衰) は本調査では未測定
  (§9 の小テストで `max|p-p_0|/p_0` を測る対象)。

---

## 9. 小テスト案 (多成分 contact advection, コード改変なし)

目的: mixed-order face-state が偽圧振動を生むかを最小構成で定量化する。

- 構成: 1D (または擬 1D) 移流。**左右で `p` 一定・`u` 一定**、`ρ,T,Y` のみ不連続
  (例: 左 He / 右 N2, あるいは異なる組成比)。理想は contact が `p` を乱さず平行移流。
- 既存ケースの活用候補: `case/05.sod_shock_tube/` の TP 派生 run
  (`run_0003_slau_tp_n2_m2regr` 等) を雛型に、初期場を「p,u 一定・Y 不連続」へ置換した
  新規 `run_*` を作る (既存 run は使い回さない)。
- sweep:
  - 再構成次数: 1 次 (`limiter`/MUSCL off) vs 2 次
  - 時間積分: explicit vs implicit (block-DPLUR)
  - CFL: 例 0.5 / 1 / 2 / 5 / 10 を sweep
- 測定量: 各ステップで `max|p - p_0|/p_0` (p_0=初期一様圧) を記録。limit-cycle なら
  時間方向に減衰せず一定振幅で持続する。`residual_history.csv` 全列 +
  `tools/check_convergence.py` で収束/未収束を判定。
- 期待: 2 次 + 高 CFL で `max|p-p_0|/p_0` が下げ止まり (limit-cycle)、1 次では消える、なら
  forge の mixed-order が主因と切り分けられる。§6-1 の同次再構成適用前後で比較する。

> 本テストはコード改変なしで初期場・config だけで構成可能。投入時は run パスを
> リポジトリルート相対で明示し、case の README run 一覧表を同期すること (AGENTS.md)。

---

## 10. 結論 8 問への直接回答

1. **同じ mixed-order face-state があるか** → **No**。species も flow と同次 2 次。
2. **species と flow は同一 face state から flux か** → **Yes** (`SetPrimitive(Primitive_i/j)` 全成分)。
3. **多成分 contact 偽圧振動の特別処理** → **No** (専用補正なし)。
4. **double-flux / quasi-conservative** → **No** (全 grep でヒットなし)。
5. **高 CFL contact limit-cycle 抑制機構** → 専用なし。間接策 = 同次 co-reconstruction +
   非物理 face の 1 次 revert + T/ρ_s/P clamp + full-coupled 近似 Jacobian + residual-based CFL。
6. **implicit Jacobian の coupling** → full `nVar×nVar`、species⇔flow を `dPdU` で連成、
   chemistry/vib/viscous の source Jacobian あり、近似 (frozen-reconstruction)。
7. **forge へ最小取り入れ** → ①species 同次再構成 (最優先) ②非物理 face 1 次 revert
   ③`∂P/∂ρ_s` full-coupled Jacobian。
8. **それでも double-flux が必要か** → **Yes の可能性高い**。SU2 も double-flux 非保有で
   偽圧振動に晒される。forge 独自に検討すべき根治レバー。

## 変更ログ

- `2026-06-20`: 初版。SU2 v8.5.0 (`12eb826`) ソースを `.external/su2-src/` に取得し
  CNEMOEulerVariable / CNEMOEulerSolver / NEMO convection (roe/ausm_slau/msw/lax) /
  CNEMONumerics / NEMO_sources を read-only 調査。コード変更なし。

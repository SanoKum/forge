# KEEP スキーム復活 (modern API・cell/node) + node WALE — LES/ILES 構成

## メタ

- **area**: `convection / turbulence`
- **status**: `in-progress`
- **related_docs**: [`methods/convection/implementation.md`](../../methods/convection/implementation.md), [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md)
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 0. 現在地 (2026-07-19 整合メモ — 着手前に必読)

- KEEP 復活 (cell/node, modern API)・純化・pRef/advGauge 対応まで**完了**。陰解法化は §7 で
  「新規実装不要 (Roe-Jacobian LHS 流用で可・TGV 検証済)」と決着済み。
- **§1 の「散逸は WALE が担う」構成はその後の検証で更新された**: 64³ TGV では **WALE off の
  KEEP + ES 散逸 (matrix, `keepDissJump: 2`, σ=0.02) が最良** ([`turbulence-wale-fix`](../accepted/turbulence-wale-fix.md),
  [`convection-keep-diss-recon-jump`](../accepted/convection-keep-diss-recon-jump.md))。WALE は修正後も
  遷移流で ν_t を早く立てすぎる。未解像高 Re 乱流向けの明示 SGS は σ-model
  ([`turbulence-sigma-model`](turbulence-sigma-model.md)) が第一候補として検証中。
- したがって**本 plan を根拠に「KEEP+WALE」構成を新規計算に採用しないこと**。現在の推奨 LES 構成は
  memory [[keep-es-dissipation-status]] / [[wale-inactive-fix]] を参照。
- 残作業は「実ケース (backstep 等) での node KEEP LES 定量検証」(§5) のみ。σ-model 検証と合流して
  実施するのが自然で、それを別 plan で扱うなら本 plan はクローズ可。

## 1. 目的

node-centered (median-dual) で **KEEP 対流スキーム + WALE SGS** による LES/ILES を回せるようにする。
KEEP は運動エネルギー・エントロピー保存の低散逸中心流束で、LES の乱流エネルギースペクトルを
損なわない。散逸は SGS (WALE) が物理的に担う、という構成を狙う。

## 2. 出発点 (調査結果)

- **KEEP**: `cfg.solver` は `SLAU/SLAU2/HLLE/ROE` のみ有効。KEEP カーネルは
  `cuda_forge/convection/legacy/convectiveFlux_ausm_keep_d.inc.cuh` の `KEEP_FVS_d` に
  **旧フラット引数 API・`ip<nNormalPlanes` 限定**で退避され dispatch から到達不能。
  さらに**真の KEEP 中心流束 (Ctilde/Mtilde/Ktilde/Itilde) は `if(false)` で無効**、
  単純中心平均 + Roe 行列散逸のみ稼働していた。
- **WALE**: `WALE_d` は現役だが `dimGrid_cell` 起動で **cell 専用**。node では起動されず
  `vis_turb=0` (= LES OFF) になっていた。
- 対流 dispatch (`convectiveFlux_d_wrapper`) は `FaceGeom`+`convPlaneBound` で cell/node を
  吸収するので、KEEP を modern API で復活させれば node でも自動で効く。

## 3. 設計・実装

### 3.1 KEEP_d (新規 `convectiveFlux_keep_d.inc.cuh`)
- `KEEP_FVS_d` を modern bundled API (`FaceGeom/PrimState/ResidualOut/LimiterFields/GradFields`)
  へ移植。シグネチャは `ROE_d` に倣う (`conv_scheme, limit_scheme, ga, cnd, geom, st, reso, lim, grd`)。
- **KEEP 中心流束を有効化** (legacy の `if(false)` ブロックを常時実行に変更)。
  - 中心流束は隣接対 `(ic0,ic1)` の**生値**で構成 (KE/エントロピー保存)。
  - 散逸は MUSCL 再構成 L/R 状態の Roe 行列 $R|\Lambda|L$ を L/R で評価し $0.5(|A_L|Q_L-|A_R|Q_R)$。
- 周回面 = `geom.nLoopPlanes` (= `convPlaneBound`)。境界 ghost (`ic1>=nCells`) は `conv_scheme=-1` で 1 次化。
- `massflux[ip]` に散逸込み総質量流束を書く (スカラー輸送整合)。
- dispatch: `convectiveFlux_d.cu` に `else if (cfg.solver=="KEEP")` を追加。cell モードの
  `skipBoundaryFluxKernel` に KEEP を追加 (主ループが境界 ghost を処理するので専用境界カーネルは skip)。

### 3.2 node WALE
- `turbulent_viscosity_d_wrapper` の `WALE_d` 起動グリッドを `dimGrid_cell` → `dimGrid_normalcell`
  に変更 (SST と同じ)。WALE_d 本体は DOF ごとの勾配/体積/wall_dist を読むだけで discretization 非依存。

## 4. スコープ

- **やる**: KEEP_d (cell/node, KEEP 中心流束 + Roe 散逸) の復活と dispatch、node WALE 有効化、
  backstep での node+KEEP+WALE LES 計算。
- **やらない (おいおい対処)**: Roe 散逸の Ducros スケーリング/低マッハ補正による LES 用の散逸低減、
  KEEP_SLAU / AUSM 系の復活。KEEP の陰解法化 (現状 explicit RK3 のみ) は方針を §7 に記録、実装は別タスク。

## 5. 検証

- backstep 3D node + KEEP + WALE (LESorRANS=1, LESmodel=1) を非定常で回し、発散しないこと・
  剥離せん断層/再付着/spanwise 乱れが発達することを確認する。SLAU ILES (run_node3d_iles_unsteady) と
  剥離構造・再付着長を比較。
- KEEP 中心流束の KE 保存性は一様流 free-stream 保持 + 等方性乱流減衰の挙動で確認 (段階)。

## 6. 影響範囲

- 新規: `cuda_forge/convection/convectiveFlux_keep_d.inc.cuh`。
- 変更: `convectiveFlux_d.cu` (include + dispatch + skip リスト)、`turbulent_viscosity_d.cu` (node WALE グリッド)。
- 既存挙動: `solver != "KEEP"` かつ既存 LES/RANS は**不変** (KEEP は新規分岐、WALE は cell でも normalcell で同一網羅)。

## 7. KEEP の陰解法化 (DP-LUR) — flux Jacobian の方針

KEEP を block-DPLUR / dual-time で陰解法化するときの LHS flux Jacobian の設計判断。実装は別タスクだが、
方針として確定し記録する (誤って「KEEP を線形化したヤコビアンを LHS に使う」方向へ進まないため)。

### 7.1 結論

**KEEP の流束をそのまま線形化したヤコビアンを LHS に使ってはいけない。LHS には既存の散逸付き
Roe 分割ヤコビアン $A^\pm$ ([`block_dplur_jacobian_d.cuh`](../../solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh) の
`accumulate_split_jacobian_cf`) を Roe/SLAU と同じくそのまま流用する。** KEEP 専用の新規ヤコビアンは不要。

### 7.2 理由

- **KEEP 厳密ヤコビアンは DP-LUR を破綻させる**: KEEP は純中心・無散逸なので $\partial F/\partial Q$ の
  固有値が純虚 (散逸ゼロ)。対角優位にならず、対角優位前提の DP-LUR / LU-SGS 内部 sweep が収束しない。
  DP-LUR の対角ブロック $D_i=\tfrac{V}{\Delta\tau}I+\sum_f A^+_f S_f$ が安定するのは $|A|$ 散逸が対角を
  支配するためで、純中心ヤコビアンにはそれが無い。厳密 Newton をやるなら DP-LUR でなく Krylov (GMRES) が要る。
- **defect-correction 構成**: RHS = KEEP (中心), LHS = 一次近似の散逸ヤコビアン $A^\pm=R\Lambda^\pm L$ という
  非対称構成にする。中心系/歪対称系の陰解法化の標準手法。
  $$\Big[\tfrac{V}{\Delta\tau}I+\sum_f A^+_f S_f+\text{粘性}\Big]\Delta Q_i=-R_i^{\text{KEEP}}-\sum_f A^-_f S_f\,\Delta Q_{\text{nbr}}$$

### 7.3 KEEP の保存性は壊れない (最重要)

擬似時間 (定常) / dual-time 内部反復が収束すると $\Delta Q\to0$ → 方程式は $R_i^{\text{KEEP}}=0$ に帰着し、
解は RHS の KEEP 中心流束だけで決まる。LHS の $|A|$ 散逸は $\Delta Q$ に掛かるので**収束点では寄与ゼロ**。
→ 物理解に散逸を混入させず、KE/エントロピー保存を維持したまま陰解法の安定性だけを買える。

### 7.4 LES/非定常での要件 (落とし穴)

dual-time で時間精度を要する LES では:
- 各物理ステップで内部反復を**十分収束させれば** $\Delta Q\to0$ となり LHS 散逸は漏れない。
- 内部反復が**不十分** (nStepInner 過少 / cfl_pseudo 過大) だと $\Delta Q\neq0$ が残り、LHS の数値散逸が
  物理場に漏れて LES のエネルギー保存を汚染する。
- → 検証では **Taylor-Green の KE/エントロピー保存量を内部反復回数の関数として測り**、陽解法 KEEP と
  同等の保存性が出る内部反復数を確認する。定常用途ならこの心配は不要。
- **(2026-06-29 追記・要修正)** 実測の結果、この「内部反復不足→散逸漏れ」は副次的で、大 dt の KE 劣化は
  **BDF2 時間離散誤差が支配**だった。詳細と訂正は §7.6 を参照。

### 7.5 付随事項

- **EOS 整合**: `accumulate_split_jacobian_cf` の `thermallyPerfect` 分岐をそのまま使い (CPG=閉形式 /
  TP=$\chi_{\text{eos}}=c^2-\kappa h$ 補正)、RHS/LHS を同一セル・同一 thermo 評価で揃える。
- **低マッハ**: KEEP は純中心ゆえ odd-even 市松が出やすい ([[backstep-lowmach-checkerboard-precond2]])。
  `lowMachPrecond=2` ($\Gamma^{-1}A$ 前処理) は LHS 操作なので §7.3 の理由で保存性に影響せず、併用が有効。
- **scalar 版** (`blockDPLUR=0`, スペクトル半径対角) も「散逸付き近似ヤコビアン」の一種なので KEEP の
  LHS として原理上使える (安定 CFL は低い)。起動・フォールバック用。
- **実装の実体**: dispatch で KEEP に対しても Roe と同じ block-DPLUR Jacobian 蓄積を呼ぶだけで、
  新規ヤコビアンコードは不要。

### 7.6 検証結果 (2026-06-29, Taylor-Green 物理CFL掃引)

コードを読んだ結果、**実装は不要** (`solver=KEEP` + `dualTime=1` + `blockDPLUR=1`, `timeIntegration=11` の
config だけで動く: block-DPLUR の Jacobian 蓄積も RHS 残差も陰解法経路も `cfg.solver` に依存しないことを確認)。
TGV 32³ cell・非粘性で explicit(物理CFL≈0.05) vs implicit を物理CFLを上げて掃引した
([`case/09.Taylor-Green/`](../../case/09.Taylor-Green/) の `run_0011`〜`run_0021`、解析 `analyze_implicit_sweep.py`)。

- **安定性**: KEEP+block-DPLUR は **物理CFL≈16 (explicit の320倍) まで NaN なしで完走**。Roe-Jacobian LHS の
  流用で陰解法が問題なく動く (§7.1–7.2 を実証)。
- **§7.3 を実証**: 物理CFL を explicit と一致させた run は K/K0=1.0033 で **explicit と4桁一致**、エントロピー変化も
  同オーダ。**LHS の Roe 散逸は収束解を汚染しない**。CFL≤2 まで explicit と KE/エントロピー保存が一致。
- **§7.4 の機構は要修正 (重要)**: 大 dt で KE は単調減衰 (CFL16 で K/K0=0.954, エントロピー +1.8e-3) するが、
  これは当初 §7.4 で書いた「内部反復不足による LHS 散逸漏れ」**ではなく、外側 BDF2 dual-time の時間離散誤差が支配**
  であることが判明した。根拠: CFL2 で内部残差は 6.8e-5→3e-9 と完全収束し `nSubIterDualTime` 20→50 で結果が
  **完全一致**、CFL16 (内部残差は ~2e-7 で頭打ち) でも 20→80 で KE が回復しない (劣化は irreducible)。サブ反復を
  増やしても改善しない=LHS 散逸漏れではなく時間離散誤差。**したがって LHS の Jacobian 選択は解の質にほぼ無関係**で、
  大 dt の精度限界は時間積分 (BDF2) が決める。LES では物理 dt を時間精度が許す範囲に保つことが要件 (理由は
  「内部反復を十分回す」ことではなく「物理 dt を上げすぎない」こと)。
- 結論: KEEP の陰解法化に新規実装は不要。Roe-Jacobian を LHS 流用する方針 (§7.1) は妥当で、収束解の保存性も
  確認できた。残課題は LES 適用時の許容物理 dt (BDF2 時間精度) の見極めであり、Jacobian 設計の問題ではない。

## 8. 変更ログ

- 2026-06-28: KEEP_d 移植 + dispatch、node WALE 有効化、methods 更新。ビルド・backstep 計算で検証中。
- 2026-06-28: **KEEP_d を純粋 KEEP に簡素化** (user 依頼)。Roe 行列散逸・MUSCL 再構成・リミタ・Ducros・
  `keepDissipation` 切替を撤去し、引数を `KEEP_d(ga, geom, st, reso)` に縮約。`space.keepDissipation`
  config も廃止 (残っていても無視)。検証: Taylor-Green M0.4 で cell・node とも運動量~1e-7・KE0.4%・
  エントロピー~1e-5 保存 (`case/09`)、回帰 smoke PASS。**注意**: 低マッハ checkerboard 抑制の散逸が無くなったため、
  homogeneous 方向の市松は `lowMachPrecond`/SLAU 併用で対処する (旧 keepDissipation=1 経路は廃止)。

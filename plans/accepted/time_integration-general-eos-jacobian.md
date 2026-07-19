# block-DPLUR の一般EOSヤコビアン化 (TP gas の陰解法律速の根治)

## メタ

- **area**: `time_integration` / `thermophysics`
- **status**: `done` (TP cfl 上限 2→100・残差 9.6e-8→4e-11、CPG ビット不変。§9)
- **related_docs**:
  - `methods/time_integration/theory.md` (block-DPLUR A± 分割の一般EOS化)
  - `methods/time_integration/implementation.md`
  - `methods/thermophysics.md` (EOS 微分 κ, χ)
- **related_plans**:
  - `thermophysics-multicomponent-tpgas.md` (§10 / γ[ic] 化の上位修正)
  - `thermophysics-species-implicit-coupling.md` (要因2。本件はその上位律速)
- **created**: `2026-06-20`
- **owner**: `CFD Dev`

## 1. 目的

`thermalMethod==2` (TP/NASA-9) の陰解法 (block-DPLUR) が `cfl_pseudo`≈2-3 で頭打ち・残差プラトーになる
**真因**を根治する。真因は、5×5 flux Jacobian の閉形式 FVS 固有系が **CPG 専用の EOS 関係 `h=c²/(γ−1)`・
`∂P/∂ρ|_e=0` を前提**にしており、TP では別方程式を線形化していること (FD 検証で TP 誤差 469% = 単なる近似精度
ではない)。**一般EOS の固有系**を実 `H_t`・`κ=∂P/∂(ρe)|_ρ`・`χ=∂P/∂ρ|_{ρe}` で構築し、LHS を真の ∂F/∂Q に
整合させる。狙い = TP で CPG 並みの cfl 上限・収束率を回復する (CPG は同一ノズル+SST で cfl100 まで rms 1e-12)。

## 2. スコープ

- **やる**: block-DPLUR の 5×5 ブロック Jacobian (A± 分割) を一般EOS固有系に置換 (TP 経路のみ)。各セルの
  `H_t, c², κ, χ` を同一セル状態・同一 thermo 評価から作り、まず **数値 L=R⁻¹ (部分ピボット double 5×5 LU)** で
  解く (Method A)。デバイス単体検証 (LR=I / FD / split / 行列作用) と CPG 回帰・ノズル回帰。
- **やらない**: CPG 経路は閉形式のまま保持 (回帰基準・ビット不変)。閉形式の一般EOS左固有ベクトル化 (Method B) は
  正しさ確認後の最適化として後続。残差床の EOS 反転誤差側の根治 (別トラックで切り分けのみ)。低マッハ前処理
  (`lowMachPrecond=2`) の一般EOS化は後続 (まず非前処理経路)。

## 3. 関連 docs と前提

- 真因の確定: `case/16.nozzle_wys/README.md`「ノズル CFL 上限」(run_0121〜0155, CPG/TP 対照) と
  `/tmp/eos_eig_verify.py` の FD 検証 (一般EOS固有系が TP の真の ∂F/∂Q に機械精度一致, 現行 CPG 形は 469% 誤差)。
- 既存の閉形式 FVS: `cuda_forge/timeIntegration_d.cu` の `accumulate_split_jacobian_cf` (lowMachPrecond 0/1 経路) と
  `build_jacobian_split` (legacy/precond)・`implicit_defect_correction_block_precond_d` (lowMachPrecond 2)。
- 既存 thermo: `thermo_d.cuh` の `thermo_cph_mix`(cp,h 融合)・`thermo_R_mix`・`thermo_gamma_mix`。
- 既に per-cell `gamma[ic]=γ_mix(T)` 化済 (κ=γ−1 は正しい)。**唯一不足はエンタルピー項と χ 項**。

## 4. 設計方針

詳細は [methods/time_integration/theory.md](../../methods/time_integration/theory.md)。要点と**実装条件 (レビュー指摘)**:

### 4.1 一般EOS固有系 (FD 検証済)

法線 `n`、接線 `t1,t2`、`U_n=u·n`。固有値 `U_n−c, U_n, U_n, U_n, U_n+c`。右固有ベクトル:
- 音響: `r∓ = [1, u∓c·n, H_t ∓ c·U_n]`
- 接触: `r_c = [1, u, H_t − c²/κ]`   ← **CPG で `H_t−c²/κ=K` に簡約** (χ=0)。一般EOS補正の本体。
- せん断: `r_s1=[0, t1, u·t1]`, `r_s2=[0, t2, u·t2]`

`H_t = h(T)+½|u|²` は**実 NASA エンタルピー**。`κ=γ(T)−1=R/c_v`。`χ=c²−κ·h` (CPG=0)。
恒等式 `c²=χ+κh`、`H_t−c²/κ = K−χ/κ`。

### 4.2 実装条件 (plan 明記必須)

1. **CPG 経路は閉形式のまま分岐保持** (`thermalMethod==0 → 現行閉形式 / ==2 → 一般EOS+数値逆行列`)。
   数値 LU は CPG 閉形式と演算順序が違うため**ビット一致は期待しない**; CPG 経路を残すことが回帰基準になる。
2. **R 構築・5×5 LU・L=R⁻¹・RΛL 積算は (試作段階で) double**。`H_t−c²/κ=K−χ/κ` は大きな二数の差になり得るので、
   float で組むと EOS 不整合とは別の条件数由来誤差が入る。GPU 負荷評価は正しさ確認後。
3. **部分ピボット付き LU** (`PA=LU`)。ピボットなし Gauss-Jordan は不可。デバッグ時に最小ピボット・`‖LR−I‖`・
   NaN/Inf・R の簡易条件数指標を監視。エッジ状態 (M→0, U_n→0, κ 小, 高温で c_p 大) を試験に含める。
4. **同一セル・同一 thermo 評価で統一**: `c²=γ(T)RT, κ=γ(T)−1, h=e(T)+RT` を**同じ T・同じ組成・同じ NASA データ・
   同じエネルギー datum**から評価する。`sonic[ic]` は別経路、`Ht[ic]` は保存量から、`gamma[ic]` は別関数、という混在を避け、
   1 つの thermo helper が `{p,h,cp,cv,kappa,chi,a2(=c²)}` をまとめて返す形にする (`ThermoDerivatives`)。
5. **セル状態と面状態を混ぜない**: DPLUR のセル対角・非対角ブロックを凍結するなら `R_i,L_i,Λ_i` はすべて同じセル i の
   `(ρ_i,u_i,H_{t,i},c_i,κ_i,χ_i)` で作る。RHS 数値流束が面/Roe 平均でも可 (LHS は近似ヤコビアン) だが、`H_t` だけ面・
   `c,κ` はセル・`U_n` は面平均、という混成は今回直す不整合を再導入するので不可。

## 5. 実装ステップ

1. **thermo helper** (`thermo_d.cuh`): `ThermoDerivatives{p,h,cp,cv,R,kappa,chi,a2}` を 1 セル状態 (ρ,e,Y) から返す
   `thermo_derivatives_mix(...)` を追加 (`thermo_cph_mix`+`thermo_R_mix` を 1 スイープ流用、内部 double)。
2. **固有系+LU** (`timeIntegration_d.cu` or 新 `cuda_forge/eos_jacobian_d.cuh`): 一般EOS R 構築 + 部分ピボット
   double 5×5 LU で `L=R⁻¹`、`A± = R Λ± L` を返す device 関数。`accumulate_split_jacobian_cf` を
   `thermalMethod` で分岐 (CPG=現行閉形式 / TP=一般EOS)。
3. **デバイス単体検証** (host コンパイル可能な test, `tools/test_eos_jacobian.cpp`): §6 の Level1/2。
4. **block-DPLUR 配線** + Level3 (行列作用) + CPG ビット回帰。
5. **ノズル回帰** (§6 Level4) + 残差床 別トラック診断。
6. (後続) 速度問題があれば閉形式一般EOS L へ (Method B)。lowMachPrecond=2 経路の一般EOS化。

## 6. 検証 (レベル別)

- **Level1 固有系単体**: 各状態で `‖LR−I‖`・`‖RΛL−A_analytic‖/‖A_analytic‖`・`‖RΛL−A_FD‖/‖A_FD‖`。
  状態: CPG低温 / TP 250K / TP 1000K / TP 高温 NASA 区間 / M=0 / 亜音速 / 音速近傍 / 超音速 / 一般方向 n (軸非平行)。
- **Level2 split**: `A+ + A− = A`、法線反転 `A_{-n}^+ = −A_n^-`・`A_{-n}^- = −A_n^+` (ic0/ic1 法線問題の再発防止)。
- **Level3 DPLUR 行列作用**: 小メッシュで任意 v に対し `J_LHS·v` が組み立てた対角・非対角ブロックと一致 (流束 FD との
  完全一致でなく、実装した DPLUR 演算子と行列構築の一致)。
- **Level4 ノズル回帰**: `cfl_pseudo` 2/5/20/50/100 で 発散ステップ / 初期残差減少率 / 最終残差床 / 出口 Mach /
  質量流量 / 全エンタルピー保存 / Newton 反転回数 を CPG (run_0145) と TP 新旧で比較。駆動は
  `case/16.nozzle_wys` の TP-N2 hot 等エントロピー (run_0151 系)。
- **CPG 回帰**: `thermalMethod==0` が閉形式経路でビット不変 (分岐により保証)。
- **残差床は別トラック (混同しない)**: 一般EOS化で cfl 上限・収束率は改善するはず。最終床は `T→e→T` 往復誤差・
  NASA polynomial の float 誤差・Newton 停止判定・保存量→primitive 丸め・IC とソルバ EOS 差・面とセルの演算順序で
  決まり得る。**機械ゼロまでは保証しない**。TP 床を論じるときは `rms_ro` だけでなく
  `max_i |e(T_i)−e_target,i|/max(|e_target|,e_ref)` も並べ、EOS 反転床との対応を見る。
  (float CPG が 1e-12 まで落ちるのは double reduction や完全定常 IC との差分を見ている可能性に留意。)

## 7. 影響範囲

- `cuda_forge/thermo_d.cuh` (ThermoDerivatives)、`cuda_forge/timeIntegration_d.cu` (分岐+一般EOS固有系) または
  新規 `cuda_forge/eos_jacobian_d.cuh`、`tools/test_eos_jacobian.cpp` (新規)。
- CPG 経路・既存ケースは不変 (分岐で保護)。lowMachPrecond=2 は当面据え置き。

## 8. 完了条件

- [x] docs theory/impl 更新
- [x] ThermoDerivatives + 一般EOS固有系 + (検証用) 部分ピボット double LU + **閉形式 (本実装)**
- [x] Level1/2 (host) + CPG ビット回帰 合格、Level4 ノズル実証
- [x] TP ノズルで cfl 上限が CPG 並みに改善 (cfl≥5 で即発散しない) を実証
- [x] 残差床トラックの切り分け (EOS 反転床 ~4e-11 が残ると確認)
- [x] README / plan 同期、`status: done`
- [x] Level3 行列作用の自動テスト (`tools/test_eos_jacobian.cpp`、実 `accumulate_split_jacobian_cf` を共有ヘッダ化して照合)
- [x] lowMachPrecond=2 経路の一般EOS化 (構築レベル。3 箇所の CPG 仮定を修正・CPG ビット不変・TP↔CPG 挙動一致)
- [ ] (後続) precond=2 TP の **end-to-end 収束検証** (低マッハ TP ケースが必要。超音速 wys は CPG でも precond=2 発散で不可)・閉形式 L の更なる検証

## 9. 変更ログ

- `2026-06-20` — 初稿。CPG/TP 対照 (case/16 run_0121〜0155) で真因=TP の LHS 固有系不整合と確定、
  `/tmp/eos_eig_verify.py` の FD 検証で一般EOS固有系 (音響 H_t±cU_n・接触 H_t−c²/κ) が TP の真の ∂F/∂Q に
  機械精度一致・現行 CPG 形は 469% 誤差を確認。Method A (数値 L=R⁻¹) で着手。実装条件 (§4.2)・検証レベル (§6) を明記。
- `2026-06-20` — **実装・検証完了 (TP 律速を根治)**。
  - **閉形式で実装** (数値 LU は検証参照のみ)。一般EOS音響左固有ベクトルを解析導出
    (`l± = [χ+κK∓cU_n, −κu±cn, κ]`, `l±·r±=2c²`) → 既存 rank-2 構造への**3項改変**:
    `accumulate_split_jacobian_cf` の音響右エネルギー成分を実 `Ht` に、左密度成分に `χ·is` (χ=c²−κh) を追加。
    `thermallyPerfect` フラグ (`thermalMethod==2`) で分岐、CPG はビット不変。`Ht[ic]` は既存 kernel 引数を流用。
  - **検証** (`/tmp` → `tools/test_eos_jacobian.cpp`): Level1/2 — 閉形式 A± が数値 LU・FD と**機械精度一致**
    (CPG/TP, M0–3, T250–2500K; `‖clsd−num‖~1e-16`, `‖RΛL−A_FD‖~1e-8`, `‖LR−I‖~1e-12`, split/法線反転)。
  - **Level4 ノズル** (`case/16` TP-N2 hot 等エントロピー, `run_0161〜0166`): **TP cfl 上限 2→≥100**
    (旧: cfl2 で 9.6e-8 プラトー・cfl5 で step14 発散 / 新: cfl 2/5/20/50/100 すべて 3000step 完走・
    rms_ro ~4e-11)。**残差プラトーも 9.6e-8→~4e-11 (3.6 桁改善)**。CPG 回帰 (`run_0166` vs 旧 `run_0143`)
    は ~atomicAdd 一致 (thermallyPerfect=0 は旧コードと同一)。物理健全 (T∈[296,505]K, 超音速, NaN 無)。
  - **残差床は別トラックと確認**: 新 TP は ~4e-11 で下げ止まり (CPG は 1.4e-12)。これは Jacobian でなく
    NASA Newton 反転・単精度 EOS 評価の床 (レビュー予測どおり「機械ゼロまでは保証しない」)。後続課題。
- `2026-06-20` — **③ Level3 自動テスト + ② precond=2 経路の一般EOS化**。
  - **③ (完了)**: `accumulate_split_jacobian_cf` を共有ヘッダ `cuda_forge/block_dplur_jacobian_d.cuh` へ移設
    (host/device 両対応, `max`→`blkdplur_max` 三元 max でビット不変)。`tools/test_eos_jacobian.cpp` に Level3
    (実カーネル関数の `diag += face_area·A⁺`・`nbr += (−A⁻)·sdq` を検証済 `eos_split_jacobian_general_closed` と照合)
    を追加 → **machine 精度一致 (~5e-16) で PASS**。標準経路 (precond=0) は再 run で TP step1 全桁一致 (ビット不変)。
  - **② (構築レベル完了・end-to-end 未検証)**: lowMachPrecond=2 の precond カーネルが持つ **3 箇所の CPG 仮定**を特定・修正:
    (i) `build_jacobian_split` の R[4] エンタルピー再構成 → 検証済 `eos_Mg_general` で a_plus/k_off を構築 (TP 分岐),
    (ii) Γ_c の g ベクトル `Htot=c²/(γ−1)+ek` → 実 `Ht`,
    (iii) ∂p/∂Q の r ベクトル `rvec[0]=κek` → `χ_eos+κek` (χ_eos=c²−κh)。いずれも `thermallyPerfect` 分岐で
    **CPG はビット不変** (再 run で確認)。**TP precond=2 は CPG precond=2 と同一挙動** (両者 supersonic wys で step6 発散)
    = EOS 不整合は解消 (旧 CPG 形なら TP は step3 で早期発散)。**ただし precond=2 の収束を正で示す低マッハ TP ケースが無く、
    超音速 wys は CPG でも precond=2 発散のため、end-to-end の収束検証は未** (後続課題)。標準経路 (precond=0/1) が TP の
    検証済・推奨経路。
- `2026-06-20` — **precond 経路の固有系を統一 (Option A, クリーンアップ)**。`build_jacobian_split` の CPG 専用べた打ち
  R/L (RL≠I の近似) と TP 分岐を廃止し、検証済み `eos_split_jacobian_general_closed` (eos_jacobian_d.cuh) を呼ぶ
  薄いラッパに統一 (`a_plus=A⁺`, `k_off=−A⁻`)。H=実 Ht・χ_eos=c²−κh を常に使い、CPG は χ_eos≈0 で簡約。Γ_c の
  g/r ベクトルも実 Ht・χ_eos に統一し `thermallyPerfect` 分岐を precond カーネルから撤去 (標準経路 accumulate は
  CPG ビット不変のため分岐維持)。**CPG precond=2 は legacy 近似→厳密 Jacobian に変わる**が収束先不変。
  **回帰スポット確認** (`case/23.axi_nozzle` CPG inviscid precond=2, eps0.15): step4000 rms_ro 1.60e-5 (旧 legacy 1.38e-5)・
  NaN 無で同一収束域 = 回帰なし。precond=2 が標準経路と同じ検証済み固有系を共有するようになり保守性向上。
- `2026-06-20` — **リネーム整備 (ビット不変)**: `accumulate_split_jacobian_cf` のローカル名を明確化
  (`chi`→`kappa_over_c` (=(γ−1)/c, EOS の χ ではない)・`inv_chi`→`c_over_kappa`・`s2`→`inv_sqrt2`・`is`→`inv_sonic`、
  `kappa=γ−1`/`chi_eos=c²−κh` を明示)、フラグ `generalEOS`→`thermallyPerfect` (κ=γ−1 仮定=TP ideal gas 専用と明記)、
  `#ifdef FORGE_JACOBIAN_SANITY` のデバッグ検査 (c≤0/非有限/κ<κ_min) を追加。再 run で TP step1 が全桁一致を確認。

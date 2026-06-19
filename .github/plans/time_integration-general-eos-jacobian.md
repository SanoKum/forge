# block-DPLUR の一般EOSヤコビアン化 (TP gas の陰解法律速の根治)

## メタ

- **area**: `time_integration` / `thermophysics`
- **status**: `in_progress`
- **related_docs**:
  - `docs/time_integration/theory.md` (block-DPLUR A± 分割の一般EOS化)
  - `docs/time_integration/implementation.md`
  - `docs/thermophysics/theory.md` (EOS 微分 κ, χ)
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

詳細は [docs/time_integration/theory.md](../../docs/time_integration/theory.md)。要点と**実装条件 (レビュー指摘)**:

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

- [ ] docs theory/impl 更新
- [ ] ThermoDerivatives + 一般EOS固有系 + 部分ピボット double LU 実装
- [ ] Level1〜4 + CPG ビット回帰 合格
- [ ] TP ノズルで cfl 上限が CPG 並みに改善 (cfl≥5 で即発散しない) を実証
- [ ] 残差床トラックの切り分け結果を記録
- [ ] README / plan 同期、`status: done`

## 9. 変更ログ

- `2026-06-20` — 初稿。CPG/TP 対照 (case/16 run_0121〜0155) で真因=TP の LHS 固有系不整合と確定、`/tmp/eos_eig_verify.py`
  の FD 検証で一般EOS固有系 (音響 H_t±cU_n・接触 H_t−c²/κ) が TP の真の ∂F/∂Q に機械精度一致・現行 CPG 形は 469%
  誤差を確認。Method A (数値 L=R⁻¹) で着手。レビュー指摘の実装条件 (§4.2) と検証レベル (§6) を明記。

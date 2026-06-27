# 新セッション引き継ぎプロンプト — 軸対称 近軸の陰解法を混合精度で root-fix

> 新しいセッションの冒頭に、この内容を貼って開始してください。

---

## タスク

forge (CUDA 密度ベースソルバ) の **軸対称 RANS で出る「軸中心 k スパイク」を根治**する。
根本原因は特定済み: **float32 の陰解法 (block-DPLUR) が近軸第一セルの平均速度 `u_r` を収束させきれない**
こと (乱流モデルは無関係)。今回は **`papers/precision/itref_phys.pdf` (Baboulin et al. 2008, mixed-precision
iterative refinement) に沿った混合精度** —「**状態 U と残差 R は double、DPLUR 線形 solve は float**」—
で根治する。

**最初に必ず読むこと**:
1. `.github/plans/precision-mixed-axisym.md` (本タスクの実装計画・落とし穴・検証基準)
2. `.github/plans/architecture-axisym-axis-singularity.md` (根本原因の確定記録)
3. `papers/precision/itref_phys.pdf` (混合精度反復改良。Algorithm 1 が指針。
   `pdftotext papers/precision/itref_phys.pdf -` で読める)
4. `AGENTS.md` の開発フロー・収束確認ルール・`.github/forge-su2-cross-check.md`

## 確定している事実 (再調査不要・これらは潰し済み)

- 真因 = float32 陰解法が近軸で `ΔU=D⁻¹R` を潰す。**explicit float も global double も正しい**
  (laminar conical 第一セル `Uy`: float陰解=**−0.64 固着**、explicit=−14.9、double=−15.1 が正解)。
- スパイクは偽の `Uy` → 偽ひずみ → SST 生産。**フープ項・Kato–Launder・dilatationCorrection は無関係**(マスク)。
- **以下はすべて実測で無効**: 5×5 solve 倍精度 / 軸対称 source-Jac 除去 / FVS 面 Jac 倍精度 /
  DPLUR 反復まるごと倍精度 / **残差の倍精度蓄積 (`doubleResidual`)** / cfl 上下 (0.1–8) /
  `nStepInner` 増 (線形系は既に収束=FGMRES でも不可) / scalar DPLUR (cfl0.1 で収束するが固着)。
- **効くのは状態 U の倍精度 (global double) のみ**。桁落ちは「残差の atomic sum」ではなく
  「**float 状態から計算する per-face 流束の値**」に由来 → 残差を正しく出すには状態 U の double 評価が必須。
- SU2 (= 後退Euler + FGMRES + ILU0 + **double**) が無傷なのは **double 精度**ゆえで手法ではない。

## やること (実装方針)

`itref_phys.pdf` Algorithm 1 の対応付け: `x`=保存状態 `U`、`r`=非線形残差 `R(U)`、`Az=r`=DPLUR 線形 solve。
- **状態 `U` (ro,roUx..roe) と残差 `res_*` を double 配列で保持・評価**。フラックス (SLAU)・軸対称ソース・
  粘性を **double で計算** → double 残差。
- **DPLUR 線形 sweep・Jacobian・前処理・勾配・限界子・乱流・出力は float のまま** (混合精度)。
- 更新 (commit) は double で。

段階:
1. **正しさ優先で native double** の混合精度を実装し、まず効くことを確認 (下記検証)。
2. **速度が問題なら double-float (compensated, float ペア+TwoSum/TwoProduct)** へ。消費者GPU(RTX3060)は
   FP64=1/32 なので、`itref_phys.pdf` の思想どおり高精度は残差・更新に限定し、桁落ちする近軸圧力差まわりに局所化。

## 検証 (AGENTS.md の計算・収束ルールに従う)

- `case/29.bell_vs_conical/run_axis_lam_slau` と**同条件** (laminar conical, `blockDPLUR=1`,
  `cfl_pseudo=2`, `nStepInner=20`, viscMethod=1) で第一セル `Uy`(x=40mm) が **−15 近傍**に収束するか。
  - 比較対象: `run_axis_lam_slau`(float固着 −0.64)、`run_axis_lam_slau_double`(global double 正 −15.1)、
    `run_axis_lam_slau_expl`(explicit 正 −14.9)。
  - 確認スクリプト例:
    ```bash
    cd case/29.bell_vs_conical
    python3 -c "import h5py,numpy as np; f=h5py.File('<run>/res_20000.h5'); g=h5py.File('<run>/nozzle.h5'); cc=g['/CELLS/centCoords'][:].reshape(-1,3); m=np.abs(cc[:,0]-0.04)<4e-4; o=np.argsort(cc[m,1]); print('first-cell Uy=', float(f['/VALUE/Uy'][m][o][0]))"
    ```
- 効けば SST (`run_su2cmp_forge_sst` 相当) で軸中心 k が **SU2 同様「軸で最小」** (スパイク消失) を確認。
  SU2 比較は `.github/forge-su2-cross-check.md` と `case/29.bell_vs_conical/compare_forge_su2.py`。
- 回帰: case 26 flat_plate / 27 CEA / 29 推力 (mdot·λ) 悪化なし。**収束は rms_ro だけで判断しない**
  (近軸の `rms_roUy`/`rms_roK` も見る — AGENTS.md 収束ルール)。
- 速度計測は native (`.github/forge-development-environment.md`)。

## 環境・参照物

- ビルド: `solver_density_cuda/build-verify`(float, arch86 に設定済)。global double 参照ビルドは
  `solver_density_cuda/build-double`(arch86, `FORGE_CUDA_BLOCKSIZE=64` で実行。**SST は発散したので注意**)。
- `flow_float` は `solver_density_cuda/flowFormat.hpp`(float)。混合精度は全面 double 化ではなく
  状態・残差のみ double にする (global double より軽い狙い)。
- 主要ソース: `cuda_forge/convectiveFlux_d.cu`(SLAU_d), `axisymmetricSource_d.cu`,
  `viscousFlux_d.cu`, `timeIntegration_d.cu`(block-DPLUR `implicit_defect_correction_block_d`),
  `variables.{hpp,cpp}`(配列確保)。
- AGENTS.md 開発フロー: 実装前に docs/plan を更新、検証済みで commit/push (feature ブランチ
  `feature/nozzle-opt-survey`、main 直接禁止)。`case/` の run 成果物は commit しない。

## 注意

- 「軸に近いところだけ double」「ハイブリッド陽的」はユーザ却下/不可。混合精度 (状態+残差 double・
  線形 solve float) を全セル一様に適用する方針。
- このスパイクは**推力・平均流の妥当性は損なわない** (case29 で検証済) ので、速度が許容外なら
  「割り切り (近軸 k は信用しない)」も選択肢。費用対効果を見ながら進める。

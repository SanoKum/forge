# 混合精度 (iterative refinement) で軸対称 近軸の陰解法を root-fix する

## メタ

- **area**: `precision` (time_integration / architecture 横断)
- **status**: `draft`  <!-- 設計のみ。実装は新セッションで -->
- **related_docs**:
  - `docs/time_integration/` (block DPLUR)
- **related_plans**:
  - `architecture-axisym-axis-singularity.md` (**根本原因の確定記録 — 必読**)
  - `time_integration-implicit-stable-cfl.md` (block DPLUR)
- **reference**: `papers/precision/itref_phys.pdf`
  (Baboulin, Buttari, Dongarra et al., *"Accelerating Scientific Computations with
  Mixed Precision Algorithms"*, 2008 — iterative refinement の枠組み)
- **created**: `2026-06-13`
- **owner**: `CFD Dev`

## 1. 背景 (確定済みの根本原因)

`architecture-axisym-axis-singularity.md` で確定: case 29 軸対称 SST の「軸中心 k スパイク」は
**float32 の陰解法 (block-DPLUR) が近軸第一セルの平均速度 `u_r` を収束させきれない**ことが真因。
乱流モデル (フープ/Kato–Launder) は無関係。**explicit float と global double はいずれも正しい。**

**精度の切り分け (すべて実測)**:
| 倍精度化した箇所 | 結果 |
| --- | --- |
| 5×5 solve / source-Jac / FVS-Jac / DPLUR反復全体 | 無効 |
| 残差の倍精度**蓄積** (`doubleResidual`) | 無効 |
| cfl 上下 (0.1–8) / `nStepInner` 増 / scalar DPLUR | 無効 (scalar は cfl0.1 で収束するが固着) |
| **状態 (P,ρ,roU) も倍精度 = global double** | **有効** |

→ 桁落ちは「atomic sum (残差蓄積)」ではなく、**float の状態から計算する per-face 流束の値**に由来。
よって**残差 R を正しく出すには状態 U を高精度で評価する必要がある**。SU2 (double, FGMRES+ILU) が
無傷なのは double 精度ゆえで、手法 (FGMRES) ではない (`nStepInner` テストで線形系は既に収束と判明)。

## 2. 方針: iterative refinement (混合精度)

`itref_phys.pdf` Algorithm 1 (Newton 等価の反復改良):
- **残差 `r_k = b − A x_{k-1}` を高精度 (εd)** で計算。
- **解の更新 `x_k = x_{k-1} + z_k` を高精度 (εd)**。
- **分解・補正解 (`A z = r`) は低精度 (εs)**。

forge への対応付け (時間進行 = Newton/defect-correction の外側反復):
- `x` = 保存状態 `U` (ro, roUx..roe)。**double で保持・更新**。
- `r` = 非線形残差 `R(U)` (フラックス+ソース)。**double で評価** (= 状態 U が double なら正しく出る)。
- `A z = r` の補正解 = **block-DPLUR 線形 sweep。float で可** (Jacobian・前処理は低精度で良いと実測・文献とも一致)。

すなわち **「状態 U と残差 R は double、DPLUR 線形 solve は float」** の混合精度。
これが効くことは global double で実証済 (状態 double が効いた)。iterative refinement の主眼は
**最も重い部分 (forge では多数回の DPLUR sweep) を float に保ったまま double 精度の解を得る**こと。

### 消費者GPU (FP64 が遅い) への配慮 — 本 plan の肝
RTX3060 は FP64 = FP32 の 1/32。よって「残差・状態の double 化」を native FP64 で全面実装すると遅い。
段階的に検証する:
1. **まず正しさ優先で native double**: 状態 U と残差評価を double 配列で実装し、DPLUR solve は
   float のまま (既存カーネル流用) → 近軸が直るか・速度低下が許容かを測る。
2. **速度が問題なら double-float (compensated)**: 状態・流束の double を「float ペア (hi,lo) +
   誤差なし変換 (TwoSum/TwoProduct)」で FP32 ハードのまま ~double 精度に。`itref_phys.pdf` の
   思想 (高精度は残差・更新のみ) に沿い、**桁落ちする近軸の圧力差まわりに限定**してコストを抑える。

## 3. 実装ステップ (新セッション)

1. `architecture-axisym-axis-singularity.md` と `itref_phys.pdf` を読み、方針確認。
2. **state+residual を double に**: `flow_float` の全面 double 化 (global double) ではなく、
   保存量 `U` と残差 `res_*` を double 配列で持ち、フラックス/ソースを double で評価。
   勾配・限界子・乱流・出力・DPLUR sweep は float のまま (混合精度)。
   - 最小構成: SLAU 経路 + 軸対称ソース + 粘性 を double 評価 → double res、DPLUR rhs は double res、
     DPLUR solve は float、更新後の U は double で commit。
3. **検証**: `case/29.bell_vs_conical/run_axis_lam_slau` と同条件 (laminar conical, blockDPLUR=1,
   cfl_pseudo=2) で第一セル `Uy` が **−15 近傍**に収束するか (現状 float は −0.64 固着、double は −15.1)。
4. 効けば SST でも軸中心 k が SU2 同様「軸で最小」になるか確認、case 26/27/29 回帰。
5. 速度計測 (native 環境)。許容外なら double-float に置換。

## 4. 検証基準

- laminar conical 第一セル `Uy` ≈ −15 (float 固着 −0.64 から脱出)。
- SST conical 軸中心 k が SU2 同様に軸で最小 (スパイク消失)。
- 全軸対称回帰 (case 26 flat_plate / 27 CEA / 29 推力 mdot·λ) 悪化なし。
- 速度低下が許容範囲 (global double 比で改善していること)。

## 5. 落とし穴 / 既知事項 (新セッションが踏まないように)

- **残差の蓄積だけ double にしても無効** (`doubleResidual` で実証)。状態 U の double 評価が必須。
- **DPLUR solve / Jacobian の double 化だけでは無効** (B で実証)。
- **FGMRES 等への線形ソルバ変更は不要・無効** (`nStepInner` テストで線形系は収束済)。
- **cfl・scalar DPLUR・ハイブリッド陽的は不可/却下**。
- global double ビルドは `solver_density_cuda/build-double` (arch86・`FORGE_CUDA_BLOCKSIZE=64`)。
  double ビルドでの **RANS (SST) は本調査で発散**したので、SST 検証時は安定性に注意。
- 参照 run: `case/29.bell_vs_conical/run_axis_lam_slau` (float 固着), `_double` (global double 正),
  `_expl` (explicit 正), `run_su2cmp_*` (SU2 比較)。

## 6. 変更ログ

- `2026-06-13` — plan 起票 (draft)。根本原因確定 (別 plan) を受け、`itref_phys.pdf` の iterative
  refinement に沿った混合精度 (状態+残差 double / 線形 solve float) を新セッションで実装する方針。

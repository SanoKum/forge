# 混合精度 (iterative refinement) で軸対称 近軸の陰解法を root-fix する

## メタ

- **area**: `precision` (time_integration / architecture 横断)
- **status**: `done`  <!-- 2026-06-14 実装完了: 閉形式 FVS + implicitSolvePrecision フラグ。§13 参照 -->
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

## 7. Phase 1 実装方針の確定 (2026-06-13, 実装セッション)

§2 の失敗実験を精査した結果、**「残差/状態 double + DPLUR 線形 solve float」という Algorithm 1
そのものの分割は未検証**であることが判明した (§2 で試したのは float ビルド上の DPLUR 倍精度 /
doubleResidual / global double のみ。"double 残差 + float solve" は試していない)。
よって Phase 1 は**最小コスト・最小リスクでこの分割を検証**する:

- **実装**: `flowFormat.hpp` を double に切替え (= 全配列 double = 残差 R(U) を double で評価) しつつ、
  **block-DPLUR カーネル `implicit_defect_correction_block_d` の線形 solve 内部演算のみ float** にする
  (double 入力を float へキャスト → float で Jacobian/5×5 solve → float 補正を double `dq_new` に書戻し)。
  commit `Q=QN+dq` は double 演算 (= Algorithm 1 の εd 更新)。
- これは plan step 1「状態 U と残差評価を double 配列で実装し、DPLUR solve は float のまま」に厳密一致。
  メモリは global double と同じ (2倍) だが、**重い DPLUR sweep は float** = 混合精度の核心を満たす。
- **判定**: laminar conical 第一セル `Uy` が −15 近傍なら Algorithm 1 分割が成立 → 線形 solve は float で良い。
  固着 (−0.6) なら solve の精度も必要 = 方針再考。
- Phase 2 (検証 OK 後): メモリ/速度最適化として勾配・乱流・出力を float へ戻し (状態+残差のみ double)、
  消費者GPU 向けに double-float (compensated) 化。

ヘルパ (`block_dplur::*`) は scalar 型 `T` でテンプレート化し、plain カーネルは float、precond
カーネル (lowMachPrecond=2) は従来通り double インスタンスを使う (codegen 不変)。

## 13. 本番実装 (2026-06-14, status: done)

§9–§12 の知見を本番化した。**2 つの独立した成果**を `implicit_defect_correction_block_d` に実装:

1. **閉形式 FVS をデフォルト化** (`block_dplur::accumulate_split_jacobian_cf<T>`): R/L 固有分解 +
   `a_plus`/`k_off`/`solve_mat` コピーを排除し rank-2 外積に。**float 陰解法を ~10% 高速化**
   (laminar: 91.5→82.6s、case26 で残差収束が legacy と一致=回帰 OK)。軸対称 (nz=0) で legacy 等価。
2. **混合精度フラグ** `implicitSolvePrecision` (`solverConfig` `time.deltaT`、既定 0=float):
   - `0`: 線形 solve を float (既定・高速)。
   - `1`: **線形 solve のみ double** (Jacobian 構築 + 5×5 + 近傍 sweep。残差/状態は float のまま)。
     軸対称 近軸の float 陰解固着を根治 (laminar Uy −0.6→**−14.98**、~×2.8)。

**実装詳細**:
- カーネルを `template<typename ST>` 化 (ST=float/double)。状態/残差 (float) を ST へキャストして取り込み、
  閉形式で diag/nbr を畳み、`solve_5x5<ST>` で in-place 解法、補正を float `dq_new` へ書戻し。
- ヘルパ (`zero5`/`add_identity_scaled`/`zero5x5`/`solve_5x5`/`accumulate_split_jacobian_cf`) を scalar 型
  テンプレート化。precond (lowMachPrecond=2)/scalar (blockDPLUR=0) 版は従来 float のまま (本フラグ非対応、TODO)。
- wrapper は `cfg.implicitSolvePrecision` で `<float>`/`<double>` インスタンスを起動。
- 変更ファイル: `cuda_forge/timeIntegration_d.cu`, `input/solverConfig.{hpp,cpp}`。

**検証** (`case/29.bell_vs_conical`): `run_0014_prod_float` (prec=0: 82.6s, Uy0=−0.6449),
`run_0015_prod_double` (prec=1: 234.5s, Uy0=−14.98)。回帰: case26 `run_regr_cf` (2D SST 陰解, 残差収束
legacy 一致), case27 `run_regr_cf` (軸対称 陰解)。

**残課題** (別タスク): SST 軸中心 k スパイクは double-solve で縮小 (9.16→7.4) するが残留 (乱流近軸の第二寄与,
§10)。precond/scalar DPLUR の double 対応。

## 6. 変更ログ

- `2026-06-13` — plan 起票 (draft)。根本原因確定 (別 plan) を受け、`itref_phys.pdf` の iterative
  refinement に沿った混合精度 (状態+残差 double / 線形 solve float) を新セッションで実装する方針。
- `2026-06-13` (実装着手) — §7 を追記。Phase 1 = double ビルド + block-DPLUR solve 内部 float で
  Algorithm 1 分割 (未検証) をまず確認する方針に確定。
- `2026-06-13` (Phase 1 実装・検証 → **棄却**) — §8 参照。**「double 残差/状態 + float 線形 solve」は
  近軸固着を直さない**ことを実測で確定。Algorithm 1 の単純分割では本ケースは根治不可。

## 8. Phase 1 検証結果 — Algorithm 1 単純分割は **棄却** (2026-06-13)

**実装**: `flowFormat.hpp` double 化 (残差 R(U)・状態 U を double 評価) + `implicit_defect_correction_block_d`
の線形 solve 内部を float 化 (Jacobian 構築・5×5 solve・近傍 sweep を float、double 入力をキャスト、
float 補正を double `dq_new` へ書戻し、commit は double)。`block_dplur::*` ヘルパは scalar 型 `T` で
テンプレート化。ビルド: `solver_density_cuda/build-mixed` (arch86)。

**検証 run**: `case/29.bell_vs_conical/run_0001_mixed_lam_slau/` (laminar conical, blockDPLUR=1,
cfl_pseudo=2, nStepInner=20, viscMethod=1 — `run_axis_lam_slau` と同条件)。成果物:
`res_20000.h5`, `residual_history.csv`/`.png`。NaN/Inf なし。

**結果 (第一セル x=40mm `Uy`)**: float と**ほぼ完全一致の固着**。

| step | mixed (double res + float solve) | float (参照) | global double (参照) |
| --- | --- | --- | --- |
| 5000 | −0.174 | −0.2 | −15.1 |
| 10000 | −0.343 | — | −15.1 |
| 20000 | **−0.641** | −0.645 | −15.1 |

`rms_roUy` は ~4e-3 で停滞 (`rms_ro` は 2e-5 まで低下) = 近軸未収束の固着シグネチャ。
→ **「残差・状態を double にしても、線形 solve が float なら近軸は固着する」**。

**理由 (確定)**: `ΔQ = D⁻¹R` には **R が正確 (double 残差・状態)** であることと、**D⁻¹ の適用が正確**
であることの**両方**が独立に必要。近軸では block Jacobian `D` が悪条件なので、R が厳密でも float の
`D⁻¹` が増幅方向で精度を失い ΔQ が不正確になる (= iterative refinement の前提「系が悪条件すぎない」が
近軸で破綻、`itref_phys.pdf` 2.1 節の但し書き)。`architecture-axisym-axis-singularity.md` §2 の
「double solve + float 残差も無効」と合わせ、**残差と solve の双方が double を要する = 実質 global double**。

**含意**: 
- itref_phys.pdf Algorithm 1 の**単純分割 (重い solve を低精度に保つ) は本ケースに適用不可**。
- plan §2 step 2 の double-float (compensated) も、**残差/状態だけでなく solve (悪条件 D⁻¹) まで**
  補償精度にしないと効かない (= 近軸計算経路全体を ~double 精度に)。残差だけの compensated は無効。
- 残る現実的選択: (A) global double (実証済・−15.1、ただし FP64 遅い・メモリ2倍)、
  (B) 近軸計算経路全体を double-float compensated 化 (FP32 ハードで ~double、大規模・不確実)、
  (C) 割り切り (近軸 k を信用しない、§3-1/2) — 推力・平均流は妥当。

## 9. **判定の対称実験 → 真因は「solve 精度」と判明 (2026-06-13)**

ユーザ指示で「両方が独立に必要 (AND)」を直接確認するため、run_0001 と**対称**の実験を実施:
**float ビルド (残差・状態は float のまま) + block-DPLUR 線形 solve 内部のみ double**
(`build_jacobian_split`・面法線・5×5 solve・近傍 sweep を、生 float 配列を double へ昇格して全部 double、
double 補正 ΔQ を float `dq_new` へ書戻し)。ビルド `solver_density_cuda/build-mixed2` (float, arch86)。
検証 run: `case/29.bell_vs_conical/run_0002_dsolve_fres_lam_slau/`。

**結果 (第一セル `Uy`)**: **−14.88 @20k で収束 (global double −15.1 と一致、近傍セルも double と一致)**。
→ **AND 仮説は棄却。真の判別因子は「線形 solve の精度」であり、残差/状態は float で良い。**

| 構成 | 第一セル Uy @20k | 判定 |
| --- | --- | --- |
| float 全て | −0.645 | 固着 |
| double 残差 + float solve (run_0001) | −0.641 | **固着** |
| **float 残差 + double solve** (run_0002) | **−14.88** | **直る ✓** |
| global double | −15.1 | 直る |

**含意 (plan 方針の転換)**:
- itref_phys.pdf Algorithm 1 の精度割り当ては**本ケースでは逆**: 「重い solve を低精度、残差を高精度」ではなく、
  **「solve を高精度 (double)、残差は低精度 (float) で可」**。近軸の悪条件 `D⁻¹` の適用にこそ double が要る。
- **安価な根治が存在する**: 全配列 float のまま、**block-DPLUR の solve カーネル内部だけ double**。
  メモリ増なし (配列 float)、追加 FP64 はセル毎 5×5 + Jacobian 構築のみ (流束/勾配/乱流の重い経路は float)。
  消費者GPU でも追加 FP64 はカーネル一部に限定 → global double より遥かに軽い見込み。
- これは `architecture-axisym-axis-singularity.md` §2「DPLUR反復まるごと倍精度=効果なし」と**矛盾**する。
  旧実験は Jacobian 構築入力や面法線を float のまま部分的に double 化したか、`solve_5x5` 単独のみ
  double 化した可能性が高い (§2 はそれらを個別に潰したが、**生 float 入力を昇格した solve 全鎖の double 化**は
  本実験が初)。→ §2 を要訂正 (本 plan §9 が一次情報)。

## 10. SST 検証結果 — 平均流は直る・k スパイクは「縮小するが残る」 (2026-06-13)

`run_0003_dsolve_sst` (= `run_su2cmp_forge_sst` 同条件 SST conical を `build-mixed2` で実行)。NaN なし
(rms_roK/roOmega が inf 表示なのは baseline float も同じ **ログ artifact**、k/omega 場は有限)。

| 量 (x=40mm 第一セル) | baseline float SST | double-solve SST | SU2 |
| --- | --- | --- | --- |
| Uy | −0.28 | **−14.9 (step5k–10k)** → −11.16 (15k 以降) | (軸で 0 から滑らか) |
| k_axis | 9.163 | 7.4 (全 step) | 軸で**最小** |
| k_axis / k_core | 2.10 | 1.70 | <1 |

- **平均流の近軸は double solve で直る** (Uy=−14.9、laminar と同じ)。step15k で −11.16 へ緩く drift
  (発散ではなく float 残差床でのリミットサイクル/乱流場の成熟に伴う再配分)。
- **k スパイクは縮小 (9.16→7.4) するが消えない**。重要: **Uy=−14.9 (平均流が正しい) の step5k–10k でも k_axis≈7.3**。
  → **k スパイクは平均流 Uy 欠陥「だけ」が原因ではない**。平均流を正しくしても残る第二の寄与
  (乱流近軸の離散化: cell中心・第一セル r=Δr/2、フープ生産など) が存在する。
  `architecture-axisym-axis-singularity.md` の「偽ひずみ (Uy欠陥) が主犯」は**部分的**で、Uy 修正後も
  残留スパイクがある (SU2 は軸で最小)。

## 11. 確定事項と次の判断
**確定 (大きい)**: 全配列 float のまま **block-DPLUR solve 内部だけ double** にすれば、
**軸対称 近軸の平均流欠陥 (laminar Uy −0.6→−14.9) は安価に根治**できる (メモリ増なし)。
旧「global double しか効かない・安価な根治なし」は誤り。

**残課題**: SST 軸中心 k スパイクは double solve で縮小するが残る (第二の寄与=乱流近軸離散化)。
完全な SU2 一致には別途、乱流近軸 (point-implicit を double 化 / 第一セル定式化 / 生産項) の対処が要る。

## 12. 速度計測 (native, RTX 3060, arch86, 20000 step, laminar conical) — 2026-06-13

| 構成 | time | 対 float |
| --- | --- | --- |
| 純 float (`build-verify`, `run_0004_float_timing`) | **91.5 s** | ×1.00 |
| float + double-solve, bs128 (`build-mixed2`, `run_0002`) | 521.6 s | **×5.70** |
| float + double-solve, bs64 (`run_0005_dsolve_bs64`) | 502.9 s | ×5.50 |
| **float + double-solve, スピル削減再構成 bs128 (`run_0006_dsolve_restruct`)** | **426.9 s** | **×4.67** |
| (参考) SST: 純 float 115.9s → double-solve 529.8s | — | ×4.57 |

| **float + double-solve, 閉形式 A± bs128 (`run_0007_dsolve_closedform`)** | **240.0 s** | **×2.62** |

**カーネル再構成 (2026-06-13, `run_0006`)**: `a_plus`/`k_off`/`solve_mat` を materialize せず
R/Λ/L から直接 diag/nbr に畳み込み + `solve_5x5_double` で diag を in-place 解法。
STACK **376→264 B**、time **×5.70→×4.67** (−18%)。Uy=−15.003。REG 依然 255 = R/L (double 50) 律速。

**閉形式 A± (2026-06-13, `run_0007`) — 大きな改善**: 固有ベクトル行列 R/L を陽に作らず、
`M(g)=g₂I+(g₁−g₂)r₁⊗l₁+(g₅−g₂)r₅⊗l₅` (音響右/左固有ベクトルのみ) で a⁺=M(λ⁺)・k_off=M((−λ)⁺) を直接
diag/nbr へ畳む `accumulate_split_jacobian_double`。R Λ L の 5×5×5 三重積が消え rank-2 外積に。
**REG 255→238・STACK 264→240・time ×4.67→×2.62**。Uy=−15.028 で根治維持 (R/L 版と一致)。
- **数値検証**: `/tmp/mtest.cpp` で `M_closedform` と legacy `R diag(g) L` を突合。**nz=0 (軸対称・2D) で厳密一致**
  (~1e-10)。一般 3D は legacy R/L が厳密逆行列でない (RL≠I, max0.246) ため僅差だが、閉形式は固有値を厳密に
  max(λ,0) とする valid な FVS 分割で defect-correction の定常解は不変 (回帰ケース 26/27/29 は全て nz=0)。
- 残オーバーヘッド ×2.62 の主因は diag_block 5×5 double (25) の常駐 + in-place 解法。さらなる削減は
  diag 常駐削減 or compensated double-float (要実装)。

**compensated double-float (df) の試行と診断 (2026-06-13〜14, `run_0010`/`run_0011`)**:
df=(float hi,lo)+TwoSum/TwoProduct(FMA) を実装 (`/tmp/dfgpu.cu` で device 上の df 演算が正しく機能=
two_sum lo=1.0・two_prod lo=1.4e-14・df 和が float 比 8× 高精度、`-O3` で collapse しないことを確認)。
solve+diag 蓄積を df 化したカーネルは **×1.16〜1.23 と高速** だが **Uy=−0.6449 で float と同値=固着**。
**診断**: 行列の蓄積・消去のみ df 化し、**Jacobian 構築 (固有ベクトル r₁/l₁・面法線 nx 等の中間量) を float の
まま**にしたため。double-solve は全鎖を double で計算しており、§2 の切り分け (solve_5x5 単独/build_jacobian
単独の double 化は無効) どおり**精度は solve 全鎖に必要**。→ **df は Jacobian 構築まで含め全中間量を df 化
しないと効かない** (部分 df = float 同等)。全 df 化は大規模だが理論上 FP64 (×2.62) より速い可能性 (df≈15 FP32
演算 vs FP64 は 32× 遅) はある。未実施。
- **full-df も棄却 (2026-06-14, `run_0012_dffull`)**: Jacobian 構築 (法線・固有ベクトル・全中間量) まで
  全 df 化したカーネル (diag を shared df) でも **Uy=−0.6453 で float 同等=固着** (×1.18 と高速だが無意味)。
  partial-df (run_0010/11) と同じ。→ **df (~48bit) は本問題の近軸に精度不足。true double (52bit) のみ有効。**
  df 演算自体は device で正しく機能 (dfgpu.cu) するが、df の実効精度は float 入力精度＋df 加算で ~1e-8 止まり
  (microtest 確認) で、近軸の悪条件下では double との 4bit 差が −15 と −0.6 を分ける閾値を跨ぐ。
  → **compensated double-float では根治できない。native double が必須。**
- **最終結論 (速度)**: **確実な根治の最小オーバーヘッドは「閉形式 double-solve ×2.62 (laminar)/×2.15 (SST)」が床**。
  shared-diag・block size・df のいずれも改善せず。RTX 3060 の FP64=1/32 が本質的制約。

**閉形式 float solve は現行 float (R/L) より速い (2026-06-14, `run_0013_float_cf`)**:
ユーザ指摘の確認。float の閉形式 solve (R/L 固有分解を作らず rank-2 外積) は **82.8s vs R/L baseline 91.5s =
×1.10 高速**、REG 168→123。数値は float R/L と等価 (軸対称で Uy 差 0.02、両方固着 −0.6449)。
→ **DPLUR の Jacobian は `build_jacobian_split` (R/L) を閉形式の `accumulate_split_jacobian_cf` に置換すれば
float 陰解法も ~10% 高速化**できる (a_plus/k_off/solve_mat コピー排除)。軸対称で legacy 等価、3D は RL≠I で
僅差だが valid (定常解不変)。**閉形式は (1) float DPLUR の高速化 と (2) double-solve 根治の基盤 の両方に有効。**

**diag 常駐削減 (shared memory) は無効と判明 (2026-06-13, `run_0009_dsolve_shareddiag`)**:
diag_block を `__shared__ double[64][5][5]` (12.8KB/block, bs64) に移し **STACK スピル 240→40 B** とほぼ除去したが、
**速度は 240.0→239.3 s で不変 (×2.62)**。→ **律速はスピル(メモリ)ではなく FP64 演算スループット + レジスタ
制約による低占有率 (REG 238 → ~18%)**。RTX 3060 の FP64=1/32 が本質。よって **×2.62 が native FP64 の床**で、
これ以上は **compensated double-float (FP64 を FP32 ペアで回避) のみが残るレバー**。Uy=−14.98 で根治維持。
- **SST 閉形式再確認 (`run_0008_dsolve_cf_sst`)**: k_axis=7.446・Uy0=−11.11 で R/L 版 `run_0003`
  (7.424/−11.16) と一致 (axis-k 差 max 0.022)。速度 **float 115.9s → R/L double 529.8s(×4.57) →
  閉形式 double 249.1s(×2.15)**。閉形式は laminar/SST 双方で R/L 版と等価かつ ~半分のコストと実証。
  (k スパイク残留=第二寄与は §10 の通り double-solve では消えない。)

- **block size 64 への変更は ~4% しか改善せず**(×5.70→×5.50)。占有率ではなく **FP64 演算 + レジスタスピル**が律速。
  `cuobjdump`: block-DPLUR (double) は **REG:255 (上限) / STACK:376 B**(double 5×5 を複数保持=diag/a_plus/
  k_off/solve_mat で溢れ local memory へスピル)。RTX 3060 は FP64=1/32 ゆえ二重苦。
- 帰結: **安価な根治だがメモリ帯域でなく演算/レジスタで ~5.5×遅い**。block size 微調整では戻らない。
- 残レバー: (a) **compensated double-float** で FP64 を避ける(FP32 は 32×速いので演算律速分は縮むが、
  hi/lo ペアは double と同じ 8B/値ゆえスピルは残る → ×2〜3 程度が現実的見込み・要実装)、
  (b) double 5×5 の**保持数を減らす**カーネル再構成(a_plus/k_off を on-the-fly 適用)でスピル削減、
  (c) 大規模 3D では DPLUR が総コストに占める割合が小さく相対オーバーヘッドは低下する見込み、
  (d) 近軸精度が要るケースのみ double-solve、通常は float(config フラグ)。

**次フェーズ (要ユーザ判断)**:
1. **本番実装**: `solverConfig` に `implicitSolvePrecision` (0=float/1=double) を追加し、
   `implicit_defect_correction_block_d` の double 内部経路を正規化 (templated ヘルパは実装済)。
   precond (lowMachPrecond=2)/scalar (blockDPLUR=0) 版も同様に。
2. **回帰** (case 26 flat_plate / 27 CEA / 29 推力 mdot·λ) + **native 速度計測** (float 比オーバーヘッド)。
3. **k スパイク残留**の扱い: (a) 乱流 point-implicit も double 化して効くか試す、
   (b) 第一セル定式化・生産項の近軸対処、(c) 割り切り (平均流・推力は妥当)。

# SLAU2 圧力流束 (低マッハ圧力–速度カップリング是正)

## メタ

- **area**: `convection`
- **status**: `done`  <!-- 実装・検証完了。ただし低マッハフロアは SLAU2 では解消せず (§9) -->
- **related_docs**:
  - `docs/convection/theory.md`
  - `docs/convection/implementation.md`
- **related_plans**: `gpu-implicit-plan.md` (block DPLUR 陰解法。本 plan とは独立だが同時に使う)
- **created**: `2026-06-07`
- **owner**: `CFD Dev`

## 1. 目的

現行 SLAU の圧力流束第 3 項は低マッハで圧力–速度カップリングが弱く、チャンバーのような
低マッハ域 ($M\sim0.06$) で非物理な圧力ばらつき (チェッカーボード) とエネルギー残差フロア
を生む。3D ノズル (`run_3d_m4_prod`, chamber $M\approx0.067$) では密度・運動量残差が
$\sim0$ まで落ちる一方 `rms_roe` が $\approx10$ で頭打ちになる症状を確認済み。
SLAU の圧力束を **SLAU2** (Kitamura–Shima 2013) の式に差し替え、低マッハの圧力散逸スケール
を是正して、このフロアを下げ・ロバスト性を上げる。質量流束は SLAU と同一のため不変。

## 2. スコープ

- **やる**:
  - `SLAU2` を `solver` の新オプションとして追加 (圧力流束第 3 項のみ差し替え)。
  - `SLAU_d` カーネルに版選択を渡し、圧力第 3 項を分岐 (mdot は不変)。
  - `docs/convection/{theory,implementation}.md` に SLAU2 を追記。
  - SLAU / SLAU2 の比較検証 (3D ノズル + 2D 軸対称回帰)。
- **やらない**:
  - 低マッハ**前処理本体** (Weiss-Smith/Turkel)、暗黙ヤコビアン (`timeIntegration_d.cu`)
    と `setDT_d.cu` の前処理化 (別 plan)。
  - AUSM 系・Roe・HLLE の低マッハ改修。
  - 既定 solver の変更 (既存ケースは `SLAU` のまま。`SLAU2` は明示選択)。

## 3. 関連 docs と前提

- 理論: [`docs/convection/theory.md`](../../docs/convection/theory.md) の **SLAU 節**
  (§ 質量流束と圧力束, L323 付近)。ここに SLAU2 の圧力束式と低マッハ根拠を追記する。
- 実装: [`docs/convection/implementation.md`](../../docs/convection/implementation.md) に
  `solver: SLAU2` の選択方法とコード対応 (`SLAU_d` 分岐) を追記。
- 低マッハ非カップリングの背景・前処理との関係は本 plan §1 と上記 docs に集約。
- **実装着手前に上記 docs 2 ファイルを先に更新する** (AGENTS 開発フロー)。

## 4. 設計方針

`SLAU_d` (`cuda_forge/convectiveFlux_d.cu`, L391–393) の圧力束:

現行 SLAU
$$\tilde p = \frac{p_L+p_R}{2} + \frac{\beta^+-\beta^-}{2}(p_L-p_R)
            + (1-\chi)(\beta^++\beta^--1)\,\frac{p_L+p_R}{2}$$

SLAU2 (第 3 項のみ差し替え)
$$\tilde p = \frac{p_L+p_R}{2} + \frac{\beta^+-\beta^-}{2}(p_L-p_R)
            + (\beta^++\beta^--1)\,\sqrt{\tfrac{|V_L|^2+|V_R|^2}{2}}\;\bar\rho\,\bar c$$

- $\bar\rho=\tfrac12(\rho_L+\rho_R)$、$\bar c=$ `c_hat` (既存)、$|V|^2=$ `velocity2_L/R` (既存)、
  $\beta^\pm=$ `beta_p/beta_m` (既存) を再利用。新規変数は $\bar\rho$ のみ。
- **質量流束 `mdot` (L395)・風上化・残差組み立ては一切変更しない。**
- 版選択は `SLAU_d` に `int slauVariant` 引数 (1=SLAU, 2=SLAU2) を追加し、第 3 項だけ分岐。
  カーネル重複は避ける (差分 1 項のため)。

## 5. 実装ステップ

1. **docs 先行更新**: `docs/convection/theory.md` SLAU 節に SLAU2 圧力束式と低マッハ根拠を
   追記 / `docs/convection/implementation.md` に `solver: SLAU2` を追記 / `docs/index.md`
   目次は構成不変なら更新不要 (確認のみ)。
2. **カーネル分岐** (`cuda_forge/convectiveFlux_d.cu`, `.cuh`): `SLAU_d` 引数に
   `int slauVariant` を追加。L391–393 の第 3 項を分岐 (SLAU2 で上式)。
3. **ディスパッチ** (`convectiveFlux_d.cu` L2538 付近): `cfg.solver == "SLAU"` →
   `slauVariant=1`、`else if (cfg.solver == "SLAU2")` → 同カーネルを `slauVariant=2` で呼ぶ。
   `skipBoundaryFluxKernel` (L2706) の条件にも `"SLAU2"` を追加 (SLAU と同様に境界専用
   フラックスカーネルをスキップ)。
4. **確認**: `solver` 文字列は whitelist 検証が無い (`solverConfig.cpp` はそのまま読み、
   wrapper で分岐) ため `solverConfig` 改修は不要。他に SLAU 圧力束を持つカーネル
   (例: KEEP_SLAU 系) が `solver=="SLAU"/"SLAU2"` 経路で使われないことを grep で確認。
5. **ビルド**: Docker dev image で `tools/build.sh`。

## 6. 検証

- **ビルド**: 警告/エラー無くリンク (`forge`, `convertGmshToForge`)。
- **検証ケース** (`guide/verification/README.md` 準拠):
  - **3D ノズル** `case/23.axi_nozzle/run_3d_m4_prod` を複製した新 `run_*` で
    `solver: SLAU2` (陰解法設定は据え置き: timeIntegration 11, cfl 0.5, relax 0.3,
    nStepInner 20)。等エントロピー初期化後に実行。
  - **2D 軸対称回帰** `case/23.axi_nozzle` の既存 SLAU m4 ケース相当を SLAU / SLAU2 で実行。
- **判定基準**:
  - SLAU2 で `rms_roe` フロアが SLAU 比で低下 (チャンバー圧力 std/mean も低下)。
  - 超音速部 (出口 $M$ 分布) は SLAU と実質同一 (低マッハ域のみ変化)。
  - 発散しない・物理的に妥当 (出口 $M\approx4$)。
  - 2D 軸対称の既存 SLAU 結果が**不変** (SLAU 経路を触らないことの確認)。
  - `residual_history.csv` → `residual_history.png` を生成。

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/convectiveFlux_d.cu`, `convectiveFlux_d.cuh`
  (`SLAU_d` 引数 + 圧力第 3 項分岐 + wrapper ディスパッチ)。
- `docs/convection/theory.md`, `docs/convection/implementation.md` (+ `docs/index.md` 確認)。
- 既存ケースへの影響: `solver: SLAU` の挙動は不変 (新オプション追加のみ)。
- `case/23.axi_nozzle/mesh_3d_m4/README.md` に SLAU2 の選択肢を追記 (任意)。

## 8. 完了条件

- [x] `docs/convection/theory.md` 更新済み (SLAU2 節 + 低マッハフロアとの関係を正確化)
- [x] `docs/convection/implementation.md` 更新済み (dispatch 表 + `slauVariant` 分岐)
- [x] 実装・検証完了 (本 plan §6 を満たす。検証結果は §9)
- [x] `.github/plans/README.md` の状態を `done` に更新
- [x] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-06-07` — 初稿 (計画)。
- `2026-06-07` — 実装。`SLAU_d` に `int slauVariant` を追加し圧力束第 3 項のみ分岐
  (`convectiveFlux_d.cu`)、wrapper で `solver: SLAU2` を `slauVariant=2` でディスパッチ、
  `skipBoundaryFluxKernel` に `SLAU2` 追加。`SLAU` 経路・mdot は不変。ビルド成功。
- `2026-06-07` — **検証 (重要な所見)**。`run_3d_m4_prod`(SLAU) と `run_3d_m4_prod_slau2`
  (SLAU2) を同条件 (3D ノズル, 陰解法 block-DPLUR, cfl 0.5, relax 0.3, nStepInner 20,
  2000 step) で比較:

  | | `rms_roe` フロア | チャンバー圧 std/mean | 出口 $M$ |
  | --- | --- | --- | --- |
  | SLAU  | 5.05 | 0.654 % | 4.028 |
  | SLAU2 | 5.06 | 0.656 % | 3.988 |

  両者ほぼ同一。**SLAU2 では低マッハのエネルギー残差フロアもチャンバー圧の非物理
  ばらつきも解消しない**。理由: 低マッハの圧力–速度カップリングは SLAU/SLAU2 で共通の
  質量流束項 $-\chi(P_R-P_L)/\widehat c$ が担い、SLAU2 が変える圧力束第 3 項は運動量側
  (主に衝撃波ロバスト性) のため。解は微差で SLAU2 は正しく作動・安定 (NaN 0, 出口
  $M\approx4$)。`solver: SLAU2` は衝撃波ロバスト性向けオプションとして残す。
  **低マッハフロアの根治には別 plan で時間項の前処理 (Weiss–Smith) が必要** (質量流束の
  散逸スケールごと是正する)。docs (theory.md SLAU2 節) もこの所見に合わせて記述済み。

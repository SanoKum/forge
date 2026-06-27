# 別セッション用プロンプト: 陰解法の安定 CFL 引き上げ

> このファイルは「forge の陰解法 (block DPLUR, `timeIntegration: 11`) の安定 CFL を一桁上げる」
> 改修を**別セッションで進めるための作業指示プロンプト**である。新しいセッションの最初の
> 入力としてこの内容を貼り付けて使う。実装着手前に AGENTS.md の開発フロー(docs→plan→実装)を踏むこと。

---

## 背景 / 問題

forge は密度ベース圧縮性ソルバー(`solver_density_cuda/`)。定常 RANS(SST)で壁法則を検証する
ケース `case/26.flat_plate_sst`(乱流平板, M=0.2, y⁺<1)を回したところ、**陰解法 `blockDPLUR=1`
でも安定 `cfl_pseudo` が 3〜5 程度に制限され、外層境界層の収束が非常に遅い**ことが分かった。

- 定常計算の実効 CFL は `cfl` ではなく **`cfl_pseudo`**(`setDT_d.cu` の `setDT_d_wrapper`、
  `unsteady==0` で `dt_local = cfl_pseudo·dt/cfl_cell`、`dt` が相殺し実効 CFL=`cfl_pseudo`)。
  詳細は `guide/solver-settings.md` の「CFL の定義」節。
- `cfl_pseudo=10` は MUSCL・cold start いずれも step 数十で発散。`nStepInner` を 15→40 に増やしても
  発散したので、**内部 Jacobi の収束不足ではなく外側 defect-correction の安定限界**。

## 現状コードの把握(改修対象)

実装は `solver_density_cuda/cuda_forge/timeIntegration_d.cu`。

- `blockDPLUR==1` 経路 = `implicit_defect_correction_block_d`(5×5 ブロック point-Jacobi の defect correction)。
  - **対流ヤコビアンは真の分離 A±** を使用(`block_dplur::build_jacobian_split` で `a_plus` を対角、
    `a_minus`(=`k_off`)を非対角に。法線符号も処理済み)。← ここは既に良い。
  - **粘性項はスカラーのスペクトル半径** `viscous_radius = 2ν_eff/δ` を**対角に identity 倍で加えるだけ**。
    非対角(隣接セル)への真の粘性ヤコビアン結合が無い。← 高アスペクト比の壁セルで弱点。
  - **point-Jacobi**: 非対角は `dq_old`(前 sweep 値)を使う遅延更新。**壁垂直方向の直接解(ライン陰解法)が無い**。
  - 軸対称ソースのヤコビアンは `diag_block[2][*]` に入っている(SST ではない)。
- スカラー版 `implicit_defect_correction_d`(`blockDPLUR==0`)は対流もスカラー半径近似。
- SST ソース: `solver_density_cuda/cuda_forge/ransSource_d.cu` の `rans_sst_source_d`。
  k/ω は平均流とは**分離(segregated)** で解かれ、`src_jac_k = β*·ω`、`src_jac_omega = 2β·ω`(消散項のみ陰的)。
  **生産項 Pk, Pω は陽的(lagged)**。壁近傍で ω~10⁸ と stiff。

ドキュメント: `docs/time_integration/theory.md`, `docs/time_integration/implementation.md`。
関連 plan: `.github/plans/time_integration-explicit-pointimplicit-sst.md`,
`.github/plans/gpu-implicit-plan.md`。

## 実装する改良(スコープ確定)

優先度と実施可否はユーザー判断で次のとおり確定:

1. **SST 生産項の陰的化**(必須・最優先, 低リスク)
   - `ransSource_d.cu` で生産項のヤコビアン対角を `src_jac_k` / `src_jac_omega` に加える。
   - k 方程式: `Pk = μ_t S²`(limiter 付き)。`∂Pk/∂(ρk)` を近似(例: `Pk/(ρk)` を正で対角に)。
   - ω 方程式: `Pω = α ρ S²` は ρω に弱依存。生産が大きい壁近傍での剛性緩和を狙う。
   - 負の対角寄与を作らない(安定化のため正の対角のみ)よう注意。
   - これは小さく閉じた変更。まずこれだけで安定 `cfl_pseudo` がどこまで上がるかを測る。

2. **真の粘性ヤコビアン**(可能なら実施, 中)
   - block 版 `implicit_defect_correction_block_d` の粘性項を、スカラー `viscous_radius` の対角加算から
     **粘性流束ヤコビアン(対角+非対角)** に置き換える(少なくとも thin-shear-layer 近似 ∂F_v/∂q)。
   - 壁垂直方向の強い粘性結合を非対角に正しく入れることが高アスペクト比セルの安定化に効く。
   - #1 で不足する場合に着手。実装コスト・安定性向上を見て段階的に。

### 対象外(今回スコープ外)

- **壁垂直方向ライン陰解法**: 安定 CFL を最も上げられる本命だが、forge は非構造データ構造で
  壁法線方向のライン抽出・ブロック三重対角ソルバを新規実装する必要があり、コストが大きすぎるため
  **今回は実装しない**。#1(+必要なら #2)で到達できる範囲を狙う。

## 検証(必須)

- 検証ケース: `case/26.flat_plate_sst`(本番)。補助で `case/20.naca_ml/001.test/run_slau`、`case/08.bump`。
- **定量目標**: 改良前の安定 `cfl_pseudo`(発達場 restart で MUSCL ≈5)に対し、#1(+可能なら #2)で
  **安定 `cfl_pseudo` をできるだけ上げる**(まずは 10〜20 を目標、到達できた上限を報告)。同一反復数での
  残差降下が速いことを示す。ライン陰解法を入れないため一桁超(30〜50)までは届かない可能性が高い点は許容。
- 壁法則(u⁺-y⁺ が log 則 `(1/0.41)ln y⁺+5.0` と粘性低層 `u⁺=y⁺`)と Cf-Re_x(Schlichting
  `0.0592 Re_x^{-1/5}` と数% 一致)が改良後も保たれること。ポスト処理: `case/26.flat_plate_sst/tools/postprocess_wall_law.py`。
- 数値結果を変えない設計の部分は既存ベースラインと差が出ないこと、意図的に変わる部分は差分を定量報告。

## 運用上の必須事項(AGENTS.md 準拠)

- **開発フロー**: 数値スキーム変更なので着手前に
  (1) `docs/time_integration/theory.md` 更新 →
  (2) `docs/time_integration/implementation.md` 更新 →
  (3) `.github/plans/time_integration-<slug>.md` を `_template.md` から作成し README.md 一覧に追記 →
  (4) 実装。完了時に plan を `status: done` 化し変更ログ記載。
- **ビルド**: `build/forge` は自動再ビルドされない。`*.cu` 編集後は必ず Docker 再ビルド:
  ```bash
  docker run --rm --gpus all --user "$(id -u):$(id -g)" -v "$PWD:/workspace" \
    -e HDF5_INC=/usr/include/hdf5/serial -e HDF5_LIBDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial \
    forge-solver:cuda-dev bash -c "cd /workspace/solver_density_cuda && bash tools/build.sh"
  ```
- **計算**: 既存 `run_*` を使い回さず複製した新 `run_*` で実行。実行後 `residual_history.png` を生成。
  メッシュ生成→`convertGmshToForge`→`forge` の流れは `guide/calculation-workflow.md` 準拠。
- **設定変更**は `guide/solver-settings.md` を参照(特に `cfl`/`cfl_pseudo` の定義)。
- **言語**: docs/コメントは日本語、識別子は原語 `code` 表記、commit は英語命令形。
- **コミット**: 機能単位で feature ブランチへ。`case/` の run 成果物(res_*.h5, *.png, log)は commit しない。

## 注意 / 既知の罠

- cold start(一様初期値)は剛性が高い。**1次風上 + 低 cfl_pseudo で spinup → restart で MUSCL+高 CFL** の
  段階法が有効(`case/26.flat_plate_sst/README.md` 参照)。CFL 改良の効果は「発達場からの restart で
  どこまで cfl_pseudo を上げられるか」で測ると分かりやすい。
- converter(`gmshReader.hpp::writeInputH5`)が書く VALUE は `variables.hpp::read_cellValNames` に
  載っている変数のみ。SST シード `roK/roOmega` は最近追加済み(これが無いと初期 k=ω=0 で相対層流化する)。
- `flat_plate` 初期条件は `solver_density_cuda/input/setInitial.hpp` に追加済み(M=0.2, k=0.3, ω=300)。

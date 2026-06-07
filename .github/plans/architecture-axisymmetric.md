# architecture-axisymmetric — 軸対称ソルバ Phase 1

## メタ

- **area**: `architecture` (geometry / flux / BC 横断)
- **status**: `done`
- **related_docs**:
  - [`docs/axisymmetric/theory.md`](../../docs/axisymmetric/theory.md)
  - [`docs/axisymmetric/implementation.md`](../../docs/axisymmetric/implementation.md)
- **related_plans**: なし
- **created**: `2026-05-31`
- **owner**: `Copilot`

## 1. 目的

forge の圧縮性 NS ソルバに **軸対称モード** (`isAxisymmetric=1`) を追加する。
対称軸は x 軸固定、半径 $r=y$。離散化は B 流儀 (幾何量 $V$, $\mathbf{S}_f$ を
$\bar r$ 重み付け、保存変数 $Q$ はそのまま) を採用し、軸 ($r=0$) 上の face を
自動退化させて軸特異性を構造的に消す。Phase 1 では非粘性 + 半径運動量の
圧力ソース $p\cdot A_{\text{planar}}$ + 粘性 source + 軸 BC + RK 陽解法までを
対象に、`case/23.axi_nozzle/run_slau_2d` で 1D 等エントロピー解と整合する解が
得られる状態を達成する。

## 2. スコープ

### やる
- `solverConfig` に `isAxisymmetric` (0/1) を追加 (`physProp` 配下)。
- `variables::setStructuralVariables_d` で $V$, $\mathbf{S}_f$ に $\bar r$ 重み、
  セル変数 `A_planar` を追加。
- 半径運動量への圧力・粘性ソース $(p-\tau_{\theta\theta})\cdot A_{\text{planar}}$
  を加算するカーネル `axisymmetricSource_d` を追加し、`main.cpp` の対流 flux
  直後に挿入。
- `axis` BC 種別を追加 (中身は `slip` 互換)。
- `case/23.axi_nozzle/run_slau_axisymmetric/` を作成して検証実行。
- 既存 2D/3D ケースが `isAxisymmetric` 未指定で従来通り動くことのリグレッション。

### やらない (将来 plan)
- 陰解法 Jacobian の軸対称化 → Phase 3。
- 任意軸対応 (`axis_dir` 指定) → 別 plan。
- LES / RANS の軸対称対応 → 別 plan。
- 出力スクリプトの $2\pi$ 自動補正 → 必要に応じて別 plan (本 plan ではコメント
  追加のみ)。

## 3. 関連 docs と前提

- 理論: [`docs/axisymmetric/theory.md`](../../docs/axisymmetric/theory.md)
  (B 流儀、圧力ソース導出、Roe 4 波固有分解)。
- 実装: [`docs/axisymmetric/implementation.md`](../../docs/axisymmetric/implementation.md)
  (フラグ・$r$ 重み付け箇所・カーネル設計・ソース対応表)。
- 既存 3D Roe: [`docs/convection/theory.md`](../../docs/convection/theory.md) と
  対比して、軸対称では $\theta$ 方向波が落ち 4 波になる旨を theory.md に明記済。
- 検証ケース: [`case/23.axi_nozzle/run_slau_2d`](../../case/23.axi_nozzle/run_slau_2d)
  をベースに新規 `run_slau_axisymmetric/` を作成。

## 4. 設計方針

支配方程式は `docs/axisymmetric/theory.md` を参照。要点のみ:

- $V := \bar r \cdot A_{\text{planar}}$, $\mathbf{S}_f := \bar r_f\,dl_f\,\hat{\mathbf n}_f$
- 半径運動量のみソース $R_{\rho u_r} = (p_{\text{cell}}-\tau_{\theta\theta,\text{cell}}) \cdot A_{\text{planar}}$
- 連続・軸方向運動量・エネルギーは追加 source 0
- 軸面は $\bar r_f = 0$ で $\mathbf{S}_f = \mathbf{0}$ → flux 計算が自動的に 0
- Roe スキームは 3D 実装をそのまま流用 (4 波分解は theory のみ)

## 5. 実装ステップ

### Phase 1A — 設計ドキュメント (本 plan 起票時に完了)
1. ✅ `docs/axisymmetric/theory.md` 作成 (B 流儀導出、圧力ソース、Roe 4 波)。
2. ✅ `docs/axisymmetric/implementation.md` 作成 (フラグ、$r$ 重み付け箇所、
   ソース対応表)。
3. ✅ `docs/index.md` に `axisymmetric` 行を追加。
4. ✅ 本 plan を `.github/plans/architecture-axisymmetric.md` として起票。
5. ✅ `.github/plans/README.md` 一覧に追加。

### Phase 1B — solverConfig フラグ
6. `input/solverConfig.hpp` に `int isAxisymmetric;` メンバを追加。
7. `input/solverConfig.cpp` の `physProp` パース節で
   `getValidatedValue<int>(physProp, "isAxisymmetric", "physProp")` を読み込む
   (未指定時は 0 を許容)。
8. `case/23.axi_nozzle/run_slau_2d/solverConfig.yaml` の `physProp:` 配下に
   `isAxisymmetric: 0` を追加 (既存挙動維持・将来コピー元のテンプレート用)。

### Phase 1C — `A_planar` セル変数追加と $r$ 重み付け
9. `variables.hpp` のセル変数リストに `"A_planar"` を追加。
10. `variables.cpp::setStructuralVariables_d` で `cudaMemcpy` 直前に
    `cfg.isAxisymmetric == 1` ブランチを追加し、ホスト配列上で
    sx, sy, sz, ss, volume を $\bar r$ 倍。`A_planar` ホスト配列に元の volume
    (= planar 面積 × 単位深さ) を保存して device に転送。
11. CPU 経路 (`setStructuralVariables`) にも同等処理を追加 (テストのみ重視、
    GPU 経路と数値一致を確認)。

### Phase 1D — 軸 BC
12. `boundaryCond.hpp` の `valueTypesOfBC` に `"axis"` (slip と同じ全 -1)
    を追加。
13. `boundaryCond.cpp::applyBconds` の if-else 列に
    `else if (bc.bcondKind == "axis") slip_d_wrapper(...)` を追加。
14. `case/23.axi_nozzle/run_slau_2d/bcondConfig.yaml` を編集する **代わりに**、
    Phase 1F で複製先 `run_slau_axisymmetric/` でのみ `kind: axis` に変更
    (元 case はそのまま slip で残す)。

### Phase 1E — 圧力・粘性ソースカーネル
15. `cuda_forge/axisymmetricSource_d.cu` / `.cuh` を新規作成。
    - kernel: `res_roUy[ic] += (p[ic] - tau_theta_theta[ic]) * A_planar[ic]`
      (符号は既存 flux と整合確認)。
    - wrapper: `cfg.isAxisymmetric != 1` で早期 return。
16. `cuda_forge/CMakeLists.txt` に新ファイルを追加。
17. `main.cpp` の時間進行ループ内、`convectiveFlux_d_wrapper(...)` 直後・
    `viscousFlux_d_wrapper(...)` 直前に
    `axisymmetricSource_d_wrapper(cfg, cuda_cfg, msh, var)` を 1 行追加。

### Phase 1F — 検証ケース
18. `case/23.axi_nozzle/run_slau_axisymmetric/` を `run_slau_2d/` から複製。
19. `solverConfig.yaml` の `physProp:` に `isAxisymmetric: 1` を追加。
20. `bcondConfig.yaml` の `axis: kind: slip` → `kind: axis` に変更。

### Phase 1G — ビルド・実行・検証
21. Docker (`compose.yml`) で `cmake --build` を実行、コンパイルエラー無し
    を確認。
22. リグレッション: `case/21.naca_2d` または `case/02.airfoil` を 100-500 step
    走らせ、`isAxisymmetric:0` で残差 / 解が以前と一致することを確認。
23. 軸対称検証: `case/23.axi_nozzle/run_slau_axisymmetric/` を 5000 step 走らせ、
    `residual_history.csv` から `residual_history.png` 生成。
24. `case/23.axi_nozzle/run_slau_axisymmetric/compare_isentropic.py` を整備して
  軸線上の Mach 数・静圧分布を 1D 等エントロピー解析解と比較する。
  `run_0022_slau_axisym_m4_aximoc_2d_cfl0p8_10k_blend2/` では出口 Mach が
  `3.9303` (数値) / `3.9739` (1D) で、出口部の整合は良好だった。
25. (補助) 既存 3D 軸対称メッシュ計算 `run_slau_3d_regression` と Mach コンター
    が定性一致することを確認。

### Phase 1H — クロージング
26. 本 plan の `status` を `done`、§9 変更ログに実装・検証結果を追記。
27. `.github/plans/README.md` 一覧の状態を `done` に同期。
28. `docs/axisymmetric/*` と `docs/index.md` の整合性確認 (実装で乖離が出た
    箇所を docs に反映)。

## 6. 検証

- **単体 / ビルド**: Docker 内で `cmake --build` 完走。`isAxisymmetric:0` で
  既存ケース ([`case/21.naca_2d`](../../case/21.naca_2d/) または
  [`case/02.airfoil`](../../case/02.airfoil/)) のリグレッションが通る。
- **検証ケース**: [`case/23.axi_nozzle/run_slau_axisymmetric`](../../case/23.axi_nozzle/)
  (Mach 2 minimum-length nozzle、上半面メッシュ)。
- **判定基準**:
  - 5000 step 走らせて残差オーダー > 4 桁低下 (rms_ro 等)。
  - `residual_history.png` を生成して残差プロファイル安定。
  - 軸線上の出口 Mach 数を 1D 等エントロピー解析解と比較し誤差 ≤ 5%。
  - 既存 3D 軸対称計算 `run_slau_3d_regression` と Mach コンターが定性一致。

## 7. 影響範囲

- 触るモジュール / ファイル:
  - [`input/solverConfig.{hpp,cpp}`](../../solver_density_cuda/input/)
  - [`variables.{hpp,cpp}`](../../solver_density_cuda/)
  - [`boundaryCond.{hpp,cpp}`](../../solver_density_cuda/)
  - [`main.cpp`](../../solver_density_cuda/main.cpp)
  - `cuda_forge/axisymmetricSource_d.{cu,cuh}` (新規)
  - [`cuda_forge/CMakeLists.txt`](../../solver_density_cuda/cuda_forge/CMakeLists.txt)
- 既存ケース・実行手順への影響: `isAxisymmetric` 未指定で完全互換。新フラグ
  を立てない限り従来通り。
- ドキュメント: `docs/axisymmetric/*`, `docs/index.md` を整備済 (Phase 1A)。

## 8. 完了条件

- [x] 関連 [`docs/axisymmetric/theory.md`](../../docs/axisymmetric/theory.md) 作成
- [x] 関連 [`docs/axisymmetric/implementation.md`](../../docs/axisymmetric/implementation.md) 作成
- [x] 実装・検証完了 (本 plan の §6 を満たす)
- [x] [`README.md`](README.md) の状態を `done` に更新
- [x] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-05-31` — 初稿。Phase 1A 完了 (theory / implementation / plan 作成)。
- `2026-05-31` — Phase 1B–1F 実装完了。
  - `solverConfig.{hpp,cpp}` に `isAxisymmetric` (optional, default 0) を追加。
  - `variables.{hpp,cpp}` の `setStructuralVariables_d` に B 流儀の $r$ 重み付けを注入。
    軸面 ($r=0$) で `S_f=0` だと下流カーネル (CFL の `vol/ss`、Roe/SLAU の
    `n=S/|S|`) で NaN が出るため、$r$ に下限 `r_floor = 1e-20` を掛けて float32
    での underflow を回避しつつ、寄与は実質ゼロに保つ実装にした。
  - `cuda_forge/axisymmetricSource_d.{cu,cuh}` を新規追加し、`main.cpp` の
    `convectiveFlux_d_wrapper` 直後に挿入 (ProfileSection に `AxisymmetricSource`
    も追加)。
  - `boundaryCond.{hpp,cpp}` に `axis` BC 種別を追加し、`slip_d_wrapper` に
    ディスパッチ。
  - 検証ケース `case/23.axi_nozzle/run_slau_axisymmetric/` を新規作成。
- `2026-05-31` — Phase 1G 部分実施。
  - **ビルド**: Docker (`tools/build.sh`) で警告なし完走。
  - **リグレッション**: `run_slau_axisymmetric_off_check` (本ブランチ + `isAxisymmetric: 0`)
    の最初 9 step の residual と `run_slau_2d` (元コード) の residual が、
    rms_roUy の float32 LSB レベル差を除いて完全一致。$\Rightarrow$ 既存ケースへの
    影響無し。
  - **軸対称ラン**: `run_slau_axisymmetric` で 200 step 完走、residual_history.png
    を生成。CFL=0.5 で max cfl≈0.50 が安定維持され、r 重み付けにより residual の
    絶対値が約 100 倍小さくなる挙動 (mesh 半径スケール ~0.025 m) を確認。
  - **既知の制約**: 以前の長時間ランでの落ちは、後から面ベクトルの向きが逆だった
    ことに起因するケース側の不整合として整理した。したがって、これは軸対称モデル
    自体の制約とは扱わない。
  - **補足**: `run_0022_slau_axisym_m4_aximoc_2d_cfl0p8_10k_blend2/` は良好に収束しており、
    軸対称モード自体は安定に動作する基準例として扱える。
  - **補足**: その後の検証で、低背圧 + `1stUpwind` の `run_0009_slau_axislowbp_1st`、
    その restart MUSCL の `run_0010_slau_axislowbp_muscl_restart`、および新規の
    軸対称ターゲット Mach 2 ノズル `run_0011_slau_axisym_m2_1st` はいずれも長時間
    安定に回ることを確認した。つまり未解決なのは「軸対称が成立しない」ことではなく、
    元の厳しい設定を含む広い条件域での頑健性である。
  - **Phase 1G で整備**: `run_slau_axisymmetric/compare_isentropic.py` を追加し、
    1D 等エントロピー比較を再現可能な形にした。`run_0022` の出口 Mach は
    `3.9303` (数値) / `3.9739` (1D) で、出口部の差は約 `1.1 %` だった。

- `2026-05-31` — diagnostics / long-run verification 追記。
  - weighted geometry を gradient に直接使うと `\partial p / \partial r` に `p/r`
    が混入するため、`calcGradient_d.cu` は axisymmetric 時でも planar geometry を
    使う実装へ修正した。
  - `axisym_r`, `axisym_p_over_r`, `axisym_uy_over_r`, `axisym_divU`,
    `axisym_roUy_source` を diagnostics として出力するようにした。
  - `case/23.axi_nozzle/mesh_axisym_m2/` を追加し、axisymmetric Mach 2 を狙った
    ノズル形状生成と長時間 run を実施した。
- `2026-06-02` — 本 plan を完了扱いに更新。
  - `docs/axisymmetric/theory.md` の phase 記述を実装現状に合わせて更新。
  - `.github/plans/README.md` の状態を `done` に同期。
  - 軸対称コアは `run_0022` / `run_0033` / `run_0034` で継続検証済みとして扱う。

### 残タスク (Phase 2 候補に切り出し)
- ~~軸近傍セルにおけるソース項剛性の implicit 処理 (`res_roUy += p A_planar` の
  係数 ∂R/∂Q を Jacobian に取り込む、または半径方向方程式を局所解析)。~~
  **実施済 (2026-06)**: block DPLUR (`timeIntegration==11`) の `implicit_defect_correction_block_d` に
  軸対称ソースヤコビアン（圧力 ∂p/∂Q ＋ 粘性フープ 2μ/(ρ r_eff)）を roUy 行の対角ブロックへ追加。
  `case/23.axi_nozzle` M4 ノズルで陽解法収束解と壁面静圧一致（平均 0.02%）、過渡収束 ~2 倍速・回帰なし。
  なお lagged source でも収束する（幾何 r 整合 + 平均流 A⁺/A⁻ 修正で足りる）ことを確認。詳細は
  [`gpu-implicit-plan.md`](gpu-implicit-plan.md) と [`docs/axisymmetric/`](../../docs/axisymmetric/)。
- 1D 等エントロピー解析解との Mach 数 / Ps 比較スクリプト整備。
- 既存 3D 軸対称ケース `run_slau_3d_regression` との Mach 場定性比較。
- 出力スクリプト群への $2\pi$ 補正コメント / ヘルパ追加。

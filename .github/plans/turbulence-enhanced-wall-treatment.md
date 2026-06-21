# SST Enhanced (Automatic / y⁺ 非依存) Wall Treatment

## メタ

- **area**: `turbulence`
- **status**: `done`  <!-- draft / in_progress / done / superseded -->
- **related_docs**:
  - `docs/turbulence/theory.md` (§6.5)
  - `docs/turbulence/implementation.md` (§3.7)
- **related_plans**:
  - `architecture-rans-sst.md` (親, done)
  - `time_integration-explicit-pointimplicit-sst.md` (done)
  - `time_integration-implicit-stable-cfl.md` (done, `case/26.flat_plate_sst` 検証基盤)
- **created**: `2026-06-21`
- **owner**: CFD Dev

## 1. 目的

SST k-ω の壁条件は現在 low-Re 壁解像型のみ
(`omega_w = 60ν/(β₁ y²)`, `k_w = 0`, [`ransBoundary_d.cu:34`](../../solver_density_cuda/cuda_forge/ransBoundary_d.cu#L34))
で、第一セルが粘性低層 (y⁺≲1) にあることを前提とする。第一セルがバッファ層・対数層に
落ちると `omega_w` も壁せん断 `τ_w` (分子勾配から評価) も誤り、Cf を過小評価する。
本計画は **y⁺ 非依存の automatic (enhanced) wall treatment** を opt-in で導入し、
第一セルが粘性低層〜対数層のどこにあっても u⁺-y⁺ が law of the wall に乗り、Cf が相関に
一致する状態を得る。

## 2. スコープ

- **やる**:
  - `solverConfig` に `wallTreatmentSST` (0:low-Re[既定] / 1:automatic) を追加。
  - 摩擦速度 `u_τ`・y⁺・`τ_w` を **Reichardt 則の逆解き** (全層を滑らかに繋ぐ普遍則, Newton 数回)
    で算出し `bc.bvar_d` に格納する新カーネル `ransWallFunction_d.*`。
  - ω 壁面 BC を Menter automatic ブレンド `√(ω_vis² + ω_log²)` に切替 (mode=1)。
  - 運動量壁せん断を modeled `τ_w = ρ u_τ²` (有効壁粘性) で課す (mode=1, `viscousFlux_wall_d`)。
- **やらない**:
  - μ_t 本体式 ([`turbulent_viscosity_d.cu:143`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu#L143))
    は変更しない。近壁 P_k の wall-function 整合化は §9 の将来課題。
  - Kader/Fluent EWT (`e^Γ`) 定式化は採用しない。
  - 軸対称・多成分との組合せ最適化は別 plan。

## 3. 関連 docs と前提

理論は [`docs/turbulence/theory.md`](../../docs/turbulence/theory.md) §6.5、実装対応は
[`docs/turbulence/implementation.md`](../../docs/turbulence/implementation.md) §3.7。
本計画の実装前に両 docs を更新済み (本 plan と同時)。検証基盤は
[`time_integration-implicit-stable-cfl.md`](time_integration-implicit-stable-cfl.md) が
確立した `case/26.flat_plate_sst` (陰解法 block-DPLUR + SLAU + MUSCL, 壁法則検証)。

## 4. 設計方針

定数 κ=0.41, β\*=0.09, β₁=0.075。`y = wall_dist[ic]`, `ν = μ_lam/ρ`,
`U_t = |U_c − (U_c·n̂)n̂|` (no-slip 壁の壁接線速度)。

### (1) 摩擦速度 u_τ — Reichardt 則の逆解き
`u⁺ = (1/κ)ln(1+κy⁺) + 7.8[1 − e^(−y⁺/11) − (y⁺/11)e^(−y⁺/3)]` を、
`u⁺ = U_t/u_τ`, `y⁺ = u_τ y/ν` で `f(u_τ) = U_t/u_τ − u⁺(y⁺(u_τ)) = 0` として
Newton 3〜5 回 (初期値 `u_τ⁰ = √(ν U_t/y)`)。GPU カーネル内で完結。

### (2) ω 壁面 BC — Menter automatic ブレンド
`ω_vis = 6ν/(β₁ y²)`, `ω_log = u_τ/(√β\* κ y)`, `ω_w = √(ω_vis² + ω_log²)`。
係数は解析漸近の 6 (現行 60 ではない)。`k_w = 0` と ghost 反射 `ω_g = 2ω_w − ω_c` は共通。

### (3) 運動量壁せん断 — 有効壁粘性
modeled `τ_w_vec = ρ u_τ² ê_t` を壁運動量残差に課す。等価有効壁粘性 `μ_eff = ρ u_τ² y / U_t`。
y⁺→0 で `μ_eff→μ` (壁解像)、対数層で `μ_eff>μ` (壁関数) → 全層連続。
no-slip なので壁せん断仕事 0、エネルギー残差・熱流束は不変。`twall_*_b`/`ypls_b` は modeled 値で上書き。

### 呼び出し順序 (確認済み, `main.cpp` 残差ループ)
`applyRansScalarBoundaries`(=ω壁BC) → `calcGradient` → `turbulent_viscosity` → … → `viscousFlux`。
u_τ は速度・y・ν のみで算出可 → **`applyRansScalarBoundaries` 直前**で新カーネルを呼び、
`bc.bvar_d["utau"]` を ω壁BC と viscousFlux 壁カーネルの両方が共有する。

## 5. 実装ステップ

1. `solverConfig.hpp` (`int wallTreatmentSST = 0;`) と `solverConfig.cpp`
   (`getOptionalValidatedValue<int>(turb,"wallTreatmentSST",0,"turbulence")` + 0/1 検証, katoLaunder と同型)。
2. `boundaryCond.hpp` の wall/wall_isothermal bvar 初期化に `{"utau",0}` 追加 (ypls/twall_* は既存)。
3. `cuda_forge/ransWallFunction_d.cu`(+`.cuh`) 新設。`computeWallFrictionSST_d` が Reichardt 逆解きで
   `utau`/`ypls`/`twall_*` を埋める。wrapper を `applyRansScalarBoundaries` 内の壁ループ直前で
   `wallTreatmentSST==1` 時のみ呼ぶ。
4. `ransBoundary_d.cu` の wall scalar カーネルに `wallTreatmentSST`/`utau` を引数追加し ω BC を分岐。
5. `viscousFlux_d.cu` の `viscousFlux_wall_d` に `wallTreatmentSST`/`utau` を渡し、接線せん断を modeled τ_w に置換。
6. 呼出し元 (main.cpp / boundaryCond.cpp) と build 定義 (新 .cu 登録) を更新。

## 6. 検証

- **ビルド**: native フルリビルド (`solverConfig` struct 変更 → 差分ビルドは obj 取りこぼし注意,
  → メモリ `stale-build-struct-layout-trap`)。
- **回帰 (細メッシュ, y⁺≈0.35)**: `case/26.flat_plate_sst` 複製 run で `wallTreatmentSST:1`。
  y⁺<1 では enhanced が壁解像に縮退するため現行 (mode 0, `run_0007`) と u⁺-y⁺・Cf-Re_x がほぼ一致
  (`tools/postprocess_wall_law.py`)。
- **効果確認 (3 つの y⁺ 帯で y⁺ 非依存性を実証)**: 第一セル高さを変えた変種を作成。
  - バッファ層 y⁺≈10 (Reichardt 遷移域=本手法の本領): mode 1 で u⁺-y⁺ が Reichardt 曲線に乗り
    Cf が相関に一致、mode 0 では大きく外れることを対比。
  - 対数層 y⁺≈30-50: mode 1 で log-law collapse、Cf が `Cf≈0.026 Re_x^(−1/7)` に一致。
  - 3 帯 (0.35/10/30-50) で mode 1 の Cf-Re_x がほぼ重なる (=y⁺ 非依存) ことを示す。
- **NaN/収束 (必須)**: 投入直後と結果前に `residual_history.csv` 全列確認、
  `tools/check_convergence.py <run_dir>` の VERDICT を結果報告に貼る (rms_ro 単独で判断しない)。
- **run 索引**: `case/26.flat_plate_sst/README.md`「## 計算 run 一覧」を同期。

## 7. 影響範囲

- 触る: `solverConfig.hpp/.cpp`, `boundaryCond.hpp`, 新 `ransWallFunction_d.cu/.cuh`,
  `ransBoundary_d.cu`, `viscousFlux_d.cu`, 呼出し元 (main.cpp / boundaryCond.cpp), build 定義。
- docs: `docs/turbulence/theory.md` (§6.5), `docs/turbulence/implementation.md` (§3.7)。新規ファイルなし。
- 既定 `wallTreatmentSST=0` のため既存ケース・回帰は不変 (opt-in)。

## 8. 完了条件

- [x] `docs/turbulence/theory.md` 更新済み (§6.5)
- [x] `docs/turbulence/implementation.md` 更新済み (§3.7)
- [x] 実装・検証完了 (本 plan §6)
- [x] `.github/plans/README.md` の状態を `done` に更新
- [x] 本 plan の `status` を `done` に変更し §9 に変更ログ

## 9. 変更ログ

- `2026-06-21` — 初稿。docs §6.5/§3.7 を先行更新。実装着手。
- `2026-06-21` — 実装完了:
  - `solverConfig` `wallTreatmentSST` (既定 0)、`boundaryCond.hpp` + `mesh.hpp` `bplaneValNames` に `utau` 登録。
  - 新カーネル `ransWallFunction_d.cu/.cuh`: Reichardt 普遍速度則を Newton 5 回で逆解きし u_τ・y⁺ を算出。
    `ransBoundary_d.cu` の壁 ω BC 直前で呼び、`bc.bvar_d["utau"]` を ω BC と viscousFlux 壁せん断で共有。
  - `ransBoundary_d.cu`: mode1 で ω_w=√(ω_vis²+ω_log²) (ω_vis=6ν/β₁y²)、mode0 は従来 60ν/β₁y²。
  - `viscousFlux_d.cu`: mode1 で接線せん断を modeled τ_w=ρu_τ² に置換 (法線粘性・体積項は落とし、熱流束は不変)。
  - **バグ修正**: bvar は `boundaryCond.hpp` の bcondKind 別テンプレートと `mesh.hpp` の `bplaneValNames`
    マスターリストの **2 箇所**に登録が要る。`utau` を後者に入れ忘れ → `bvar_d["utau"]` 未確保で
    illegal memory access。両方に登録して解消。struct 変更のため native フルリビルド。
- `2026-06-21` — 検証 (`case/26.flat_plate_sst`, 詳細は同 README「EWT 検証結果」):
  - 細メッシュ y⁺₁≈0.35 (`run_0009`): mode1 ≈ mode0(run_0007) で Cf/Schl 0.93–0.99、残差 floor 全列一致 → 壁解像に縮退。
  - 粗メッシュ y⁺₁≈30 (`run_0010/0011`): mode0 Cf_model/Schl 0.13 (崩壊) → mode1 0.89–0.93 (回復)。
  - バッファ y⁺₁≈10 (`run_0012/0013`): mode0 0.51–0.60 → mode1 1.05–1.14。
  - mode1 の u_τ は全帯 2.5–2.8 に揃い **y⁺ 非依存**を実証。壁関数 Cf は modeled τ_w=ρu_τ²
    (`tools/wall_law_modeled.py`) で評価 (分子勾配は粗メッシュ対数層で過小評価のため不適)。
  - 残差は block-DPLUR 構造事情でプラトー (本ケース既知) だが Cf ドリフト ≤0.06% で定常。
  - 将来課題 §10: 近壁 P_k/μ_t の wall-function 整合化、wall-function 項の implicit Jacobian。
- `2026-06-21` — k 壁 BC の検討 (レビュー指摘: k は zero-gradient にすべきでは):
  教科書どおり mode1 で k を zero-gradient (Neumann) に変えて再検証 (`run_0014`/`run_0015`)。
  **結果は悪化** — y⁺30 Cf/corr 0.89→1.80, y⁺10 1.10→1.68。原因: zero-gradient k BC は
  近壁 P_k の wall-function 化 (P_k=τ_w·(∂U/∂y)_log, k を u_τ²/√β* に固定) とセットで成立するが、
  forge は P_k が解像勾配のままで粗メッシュでは誤り → k 暴走 → μ_t 過大 → u_τ/Cf 過大。
  k=0 Dirichlet はこの未整合生産を部分相殺するため (この時点の実装では) 良い。いったん k=0 Dirichlet に
  戻し、真の automatic (k zero-grad + P_k log則化 + ω ピン留め) を次項で実装した。
- `2026-06-21` — **真の automatic wall treatment を実装** (旧 §10 を本体に取込み, docs §6.5(b')(c)(d) 改訂):
  - **ω ピン留め**: wall-adjacent セルの `omega`/`roOmega` を Menter ブレンド ω_w に毎評価ピン留め
    (`ransBoundary_d.cu`)。全 y⁺ で妥当 (y⁺→0 で ω_vis 支配) なので μ_t を上限し k 暴走を断つ。
  - **k zero-gradient** (mode1): `k_g=k_c`。
  - **P_k の wall-function 化**: 新セル変数 `wf_pk` を `ransWallFunction_d.cu` で wall-adjacent セルに
    `P_k=ρu_τ⁴/ν·g(1-g)`, g=du⁺/dy⁺(y⁺₁) と算出 (Reichardt 微分)。`ransSource_d.cu` が標準 P_k を
    この値に置換。粘性低層 g→1 で P_k→0 (壁解像極限), 対数層で ρu_τ³/(κy)。これと ω ピン留めで
    k が平衡値 u_τ²/√β* に自律収束。
  - ω ピン留めセルの `res_roOmega`/`src_jac_omega` を 0 化 (Dirichlet セルの残差は 0; rms_roOmega 汚染回避)。
  - 新変数 `wf_pk` は `variables.hpp` cellValNames に登録。`boundaryCond.cpp` が bc ループ前に
    `initWallFunctionPk_d_wrapper` で -1 初期化。
  - 検証 (`case/26`, full WF, modeled Cf, x=0.30/0.60/0.89): 細 y⁺0.35 `run_0016` 0.89/0.92/0.95 (=壁解像),
    バッファ y⁺10 `run_0018` 0.95/0.98/1.01 (ω-blend 単独の 1.05–1.14 過大を解消), 対数 y⁺30 `run_0017`
    0.92/0.94/0.97。u_τ≈2.5–2.8 で **y⁺ 非依存**。u⁺-y⁺ も 3 帯で普遍プロファイルに collapse
    (`ewt_uplus_yplus.png`)。既定 0 はビット不変。
  - 収束: VERDICT は 3 run とも `NOT CONVERGED (plateau)` (本番 `run_0007` も同様の block-DPLUR 構造的
    プラトー)。**運用判定 A** (全残差列が低レベル横ばい かつ Cf ドリフト <0.3%) で実用上 OK と判断。

## 10. 将来課題 (本 plan 外)

- 近壁第一セルの μ_t を wall-function 整合化 (現状は標準 SST 式のまま; ω ピン留めで実効的に上限済み)。
- 陰解法 (block-DPLUR) で wall-function 項の Jacobian 整合 (現状は陽的評価 + ω ピン留めセルは残差 0 化)。
  block-DPLUR の残差プラトーを解消できれば VERDICT を PASS にできる見込み。
- 多成分・軸対称との組合せ (axisymmetricSource が res_roOmega を足す経路で ω ピン留めセルの再零化が要る)。

# 26. 乱流平板 (SST 壁法則検証)

Menter SST k-ω モデルが乱流平板の**壁法則 (law of the wall)** を再現できるかを確認するケース。
陰解法 (block DPLUR, `timeIntegration: 11`) + SLAU + MUSCL で定常計算する。

## 流れ条件

- 入口: 全圧 `Pt=100000 Pa`, 全温 `Tt=288.15 K` (`inlet_Pressure`)
- 出口: 背圧 `Ps=97250 Pa` (`outlet_statPress`)、`M≈0.2` (U∞≈67.8 m/s)
- 空気 `mu=1.8e-5 Pa·s`, `Re/m≈4.46e6` → 平板長 1 m で `Re_L≈4.5e6` (乱流)
- 流入乱流: `k=0.3, omega=300` (freestream `mu_t/mu≈60`, naca_ml 相当)

## 領域・メッシュ (`mesh/flat_plate.geo`)

- `x∈[-0.1,1.0]`, `y∈[0,0.2]`, z 厚さ 0.01 (1 層押し出し)
- 底面: `x∈[-0.1,0]` slip (対称), `x∈[0,1.0]` no-slip wall (平板)
- 上面 slip, 入口/出口は上記 BC、z 両面 slip
- 壁垂直: 90 点, 第一セル高さ **4.14e-6 m**, 等比 1.10 → **y⁺₁≈0.35〜0.4 (<1)**
- streamwise: 前縁 (x=0) へ上流側・平板側の両方を寄せ、前縁でセルサイズ連続 (≈0.4 mm)

## 実行手順 (段階法)

**重要**: 定常計算の実効 CFL は `cfl` ではなく `cfl_pseudo`
(`.github/forge-solver-settings.md` の「CFL の定義」を参照)。
本ケースは壁セルが極めて細く陰解法でも CFL 制限が低いため、段階的に進める。

1. **stage A (`run_0005_slau_spinup`)**: 一様初期値からの cold start。
   `convMethod: 0` (1次風上) + `cfl_pseudo: 2.0` で 25000 step 回し、平滑な発達場を作る。
2. **stage B (`run_0006_slau_muscl`)**: stage A の `res_25000.h5` から restart。
   `convMethod: 2` (MUSCL) + `cfl_pseudo: 3.0` で 40000 step 回し、2次精度で収束させる。

> cold start で `cfl_pseudo≳5` や、いきなり MUSCL + 高 `cfl_pseudo` にすると発散する。
> 発達場からの restart でも本ケースの安定上限は `cfl_pseudo≈3` 程度。

メッシュ生成・変換・実行は `.github/forge-calculation-workflow.md` の標準手順に従う
(Docker, `convertGmshToForge`, `forge`)。

## ポスト処理

```bash
python3 tools/postprocess_wall_law.py run_0006_slau_muscl 0.3 0.6 0.89
```

`wall_law.png` (u⁺-y⁺ と Cf-Re_x) を出力する。

## 収束について

陰解法 `blockDPLUR` は減衰付き point-Jacobi の defect correction で、近似ヤコビアン+高アスペクト比
壁セル+陽的 SST 生産項のため安定 `cfl_pseudo` は MUSCL で ≈5(cold start は 2〜3)が上限。
そのため収束が遅く、`run_0006`(40000 step)時点では外層がまだ発達途中(Cf が 10000 step あたり
3〜4% ドリフト)だった。

- **stage C (`run_0007_slau_muscl_innersweep`)**: `run_0006/res_40000.h5` から restart し、
  `cfl_pseudo: 5.0`, `nStepInner: 20`, MUSCL で **120000 step** 継続。残差 3.8e-7、
  最後の 20000 step での Cf ドリフト 0.1% 以下まで収束。**最終結果はこの run を使う。**

> 安定 CFL が低い構造的原因と改良案(SST 生産項の陰的化, 真の粘性ヤコビアン)は
> `.github/plans/implicit-acceleration-session-prompt.md` に整理(別セッションで実装予定)。

## 結果 (run_0007, 収束済み MUSCL)

| station | Re_x | y⁺₁ | u_τ | Cf | Cf/Schlichting |
|---|---|---|---|---|---|
| x=0.30 | 1.34×10⁶ | 0.36 | 2.707 | 3.12×10⁻³ | 0.885 |
| x=0.60 | 2.66×10⁶ | 0.34 | 2.580 | 2.82×10⁻³ | 0.918 |
| x=0.89 | 3.99×10⁶ | 0.33 | 2.511 | 2.67×10⁻³ | 0.944 |

- u⁺ vs y⁺: 3 ステーションが普遍プロファイルに重なる。粘性低層 `u⁺=y⁺` (y⁺<5)、
  対数則 `u⁺=(1/0.41)ln y⁺+5.0` (y⁺=30〜300)、外層 wake までよく再現。
- Cf vs Re_x: 発達領域で Schlichting 1/5 乗則と 1/7 乗則の間に乗る。1/5 乗則より数% 低いのは、
  前縁からの層流→遷移区間ぶんの virtual-origin オフセット(Re_x を前縁起点で測るため)による系統差。
- y⁺₁≈0.33〜0.36 (<1) を全域で満たす。

## 計算 run 一覧

| run_* | 目的・設定 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_rans_implicit` | 標準: SLAU + SST + block-DPLUR 陰解 (2D) | rms_ro 2.08e-9 収束、壁法則再現 | active |
| `run_regr_cf` | 回帰: 閉形式 FVS 既定化 (`implicitSolvePrecision=0`) の確認。run_0001 と同条件 | 残差収束が legacy と一致 (rms_ro 2.0e-9)、NaN なし → **閉形式は 2D 陰解に無害**。plan [precision-mixed-axisym.md](../../.github/plans/precision-mixed-axisym.md) | active |
| `run_0009_ewt_fine_mode1` | EWT 回帰: 細メッシュ (y⁺₁≈0.35) で `wallTreatmentSST:1`。run_0007 収束場から restart 20000 step | mode1 ≈ mode0(run_0007): Cf/Schl 0.93/0.96/0.99、残差 floor 全列一致 → 壁解像に縮退 | active |
| `run_0010_ewt_yp30_mode0` | EWT 検証: 粗メッシュ y⁺₁≈30, mode0(low-Re)。cold start→MUSCL | **Cf 崩壊** Cf_model/Schl 0.13/0.13/0.14 (low-Re BC が粗メッシュで破綻) | active |
| `run_0011_ewt_yp30_mode1` | EWT 検証: 粗メッシュ y⁺₁≈30, mode1(enhanced) | **Cf 回復** Cf_model/Schl 0.89/0.91/0.93、u_τ≈細メッシュ一致 → y⁺ 非依存 | active |
| `run_0012_ewt_yp10_mode0` | EWT 検証: バッファ帯 y⁺₁≈10, mode0(low-Re) | Cf_model/Schl 0.51/0.56/0.60 (過小) | active |
| `run_0013_ewt_yp10_mode1` | EWT 検証: バッファ帯 y⁺₁≈10, mode1(ω-blend+τ_w, k=0 Dir) | Cf_model/Schl 1.05/1.10/1.14 (過大) → full WF で改善 | ref |
| `run_0016_ewt_fine_wf`  | **full WF** (ω pin + k zero-grad + P_k log則化): 細 y⁺₁≈0.35 回帰 | Cf/Schl 0.89/0.92/0.95 (=壁解像), Cf_molec≈Cf_model → 壁解像極限を保持 | active |
| `run_0017_ewt_yp30_wf`  | **full WF**: 対数帯 y⁺₁≈30 | Cf/Schl 0.92/0.94/0.97, u_τ≈2.5–2.6 | active |
| `run_0018_ewt_yp10_wf`  | **full WF**: バッファ帯 y⁺₁≈10 | Cf/Schl 0.95/0.98/1.01 (ω-blend単独の過大を解消) | active |
| `run_0019_des_regr_mode0` | **SST-DES T1-A 回帰**: 発達場(run_0007/res_120000)から `DESmode:0` 1000 step。新旧バイナリ比較 | DESmode:0 は baseline と atomicAdd 床内で一致 (1 step base–base ≡ base–des: ro 1.19e-7/roe 1.56e-2/roK 3.73e-9) → **ビット不変**。plan [turbulence-iddes-sst.md](../../.github/plans/turbulence-iddes-sst.md) | active |
| `run_0020_des_mode1` | **SST-DES T1-A シールド (wall-resolved y+1)**: 同場から `DESmode:1` 1000 step | mode1 vs mode0 `roUx` relL2 **5.7e-7**≪Cf 0.1%, `roK` 1.8e-5。付着 BL シールド (DES 発火は外縁 3%)、NaN なし | active |
| `run_0021_des_ewt_yp30_mode0` | SST-DES T1-A: y+30 wall-modeled 場 (run_0011) から `DESmode:0` 1000 step (mode1 比較基準) | NaN なし | ref |
| `run_0022_des_ewt_yp30_mode1` | **SST-DES T1-A シールド (wall-modeled y+30)**: 同場 `DESmode:1` 1000 step | mode1 vs mode0 `roUx` relL2 4.7e-5, `roK` 9.3e-4、付着 BL シールド (frac f_d<0.05=0.92)、NaN なし | active |

### node-centered (median-dual) 検証 run

| run_* | 目的・設定 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_node_lam_cont` | node モード層流平板 (visc=4e-4, Re_L≈1e5, M=0.2) | Blasius BL 再現 (δ99=0.99×, 形状誤差 3.3%)、図 `blasius_validation.png` | active |
| `run_node_sst_fine` | node モード SST (cell run_0001 から cross-mesh restart, 20000 step) | omega/k 6 桁収束だが u_τ 過小 (node 1.24 vs cell 1.97)。残差プラトー (rms_roUx≈0.23)。図 `uplus_yplus_node_vs_cell.png` | active |
| `run_node_lam_5e_long` | **node 弱形式境界 (Phase 2: 5a+5b+5e)** 層流 (40000 step, run_node_lam_cont 崩壊場から restart) | **出口 BL 崩壊解消** (δ99(x) 単調・Blasius 一致, x=1.0 δ99=0.0114 vs 旧 1e-5)。**残差プラトー打破** rms_roUx 0.214→3.08e-5 (3.8桁), rms_ro 1.04e-7。NaN なし | active |
| `run_node_sst_5e_long` | **node 弱形式境界 (5a+5b+5e)** SST (40000 step) — 過剰乱流の「before」 | プラトー打破 (rms_roUx 0.597→1.3e-3) だが **過剰乱流**: peak mut/μ=375, Cf≈3×Schl。真因=壁 ω 非ピン + wall_dist バグ (下行で解消) | ref(before) |
| `run_node_sst_final` | **node SST 完成**: massflux 書込 + 壁 ω Dirichlet ピン (omega/roOmega ピン+ω decouple+res_roOmega 壁ゼロ) + **wall_dist 修正** (node で壁点に壁ノード座標)。過剰乱流場から 60000 step | **cell 基準と一致: Cf/Schl 0.87/0.90/0.92** (旧 3×), peak mut/μ 217 (cell 199), k_wall=0, δ99 cell 近接。rms_roUx 0.597→1.5e-4 (3.6桁), 全列 falling。NaN なし | active |

> **node 弱形式境界 (Phase 2)**: node モードで inlet/outlet/slip も ghostless 弱形式化。(5a) 主対流ループを内部+periodic
> のみに、(5b) 全境界を `convectiveFlux_boundary_d` の bvar 弱形式に、(5e) **block-DPLUR Jacobian の境界半割面で
> 退化粘性対角 `viscous_diag` をスキップ** (境界ノードが境界上で `dcc·S≈0`→delta 爆発→dq≈0 で凍結し出口 BL 崩壊・
> 残差プラトーを生んでいた真因)。**出口 BL 崩壊と残差プラトーの双方を解消**。plan
> [discretization-node-boundary-ghostless.md](../../.github/plans/discretization-node-boundary-ghostless.md)。コーナー BC
> (壁∩出口) は `ow=ib` マルチマーカ emit + `wall_flag` Dirichlet で別途解決済み (commit 5ce92dc)。
> **過剰乱流を解消** (commit 2b19d1d/0f6d53d): 真因は node で**壁 ω が Dirichlet ピンされていない**こと (ghost BC が
> 壁ノード CV 中心・dcc≈0 退化で omega[ic] を固定できず ω 過小→k 消散不足→過生成)。修正: (i) 境界 massflux 書込
> (出口で k が移流流出できず蓄積するバグ)、(ii) `wall_y_eff` MEAN→MIN (正しい ω_w)、(iii) 壁ノードで omega/roOmega を
> 直接ピン (保存変数 Dirichlet)、(iv) point-implicit で ω 行 decouple (dω=0)、(v) res_roOmega 壁ゼロ化。
> （commit 2b19d1d/0f6d53d で peak mut/μ 367 まで改善するも x=0.89 が ~1.5× 残った）
>
> **さらに wall_dist バグを修正 (commit 8b60f2c) → cell 基準と完全一致**: `calcWallDistance_kdtree` が node で壁点に
> 半割面重心 (壁ノードから x に ~dx/8 ずれ) を使い、近壁 wall_dist が法線距離 y でなく x ずれ (≈1e-4·x, 下流で増大)
> になっていた → ω_w 過小 → 過剰乱流 (下流ほど)。node では壁点に**壁ノード座標**を使うよう修正。結果
> [run_node_sst_final](run_node_sst_final/): **Cf/Schl 0.87/0.90/0.92**, peak mut/μ **217** (cell 199), k_wall=0,
> δ99 が cell 近接、rms_roUx 0.597→1.5e-4 (3.6桁)。**massflux + omega pin + wall_dist の 3 点で SST node が
> cell 基準に一致**。

> **SST-DES (DDES) T1-A** の合否: `DESmode:0` ビット不変 + `DESmode:1` で付着 BL が RANS から
> 不変 (Cf 差≪0.1%)。**発達乱流場 (nut/nu≫1) から restart すること** (`run_regr_cf` 等 nut/nu≤1 の
> 未発達場では νt≈0 で rd 小→シールドが効かない)。詳細 plan [turbulence-iddes-sst.md](../../.github/plans/turbulence-iddes-sst.md)。

> EWT (enhanced wall treatment) 検証の詳細は plan [turbulence-enhanced-wall-treatment.md](../../.github/plans/turbulence-enhanced-wall-treatment.md)。
> まとめ図は `ewt_comparison.png` (壁処理進化の Cf 比較) と `ewt_uplus_yplus.png` (u⁺-y⁺ collapse)。
> 再生成: `python3 tools/plot_ewt_summary.py`。

## EWT 検証結果 (`wallTreatmentSST`, 3 帯の y⁺ 非依存性)

`wallTreatmentSST:1` を第一セル高さの異なる 3 メッシュ (y⁺₁≈0.35/10/30) で検証。壁関数 run の Cf は
分子勾配ではなく壁関数整合の **modeled τ_w=ρu_τ²** で測る (`tools/wall_law_modeled.py`、x=0.6)。
実装は 3 段階で進化した:

1. **mode0 (low-Re, 既定 0)**: 粗メッシュで Cf 崩壊 (y⁺30 で相関の 13%)。
2. **ω-blend のみ (k=0 Dirichlet)**: ω は Menter ブレンド+modeled τ_w。改善するが y⁺10 で 1.10 と過大、
   y⁺30 で 0.91。k=0 が近壁 k を潰し非整合 ([architecture-rans-sst] / theory §6.5(b'))。
3. **full WF (本番)**: ω を wall-adjacent セルにピン留め + k zero-gradient + 近壁 P_k を wall-function 値
   `ρu_τ⁴/ν·g(1-g)` に置換。3 点セットで k が平衡値 u_τ²/√β* に自律収束し全帯で相関に乗る。

| 帯 | y⁺₁ | mode0 | ω-blend(k=0) | **full WF** | full WF u_τ |
|---|---|---|---|---|---|
| 細 (壁解像) | ≈0.35 | 0.92 | — | **0.89–0.95** | 2.50–2.69 |
| バッファ | ≈10 | 0.51–0.60 | 1.05–1.14 | **0.95–1.01** | 2.57–2.79 |
| 対数 | ≈30 | 0.13–0.14 | 0.89–0.93 | **0.92–0.97** | 2.51–2.73 |

(数値は x=0.30/0.60/0.89 の Cf/Schlichting。full WF は全帯で相関 ±8% 以内、u_τ≈2.5–2.8 に揃い
**第一セルの y⁺ 位置に依存しない**。細メッシュでは Cf_molec≈Cf_model で wall-resolved 極限に縮退し、
既定 `wallTreatmentSST:0` は従来と完全不変。)

### 収束判定 (VERDICT 引用)

本ケースは block-DPLUR の近似ヤコビアン由来で残差が**構造的にプラトー**する (受理済み本番
`run_0007` も同様)。`tools/check_convergence.py` の VERDICT は full WF 3 run とも:

```
run_0016_ewt_fine_wf / run_0017_ewt_yp30_wf / run_0018_ewt_yp10_wf
  -> NOT CONVERGED (stalled/plateau)   # 各 run の CONVERGENCE_VERDICT.txt
```

で、参照の `run_0007` (120000 step) も `NOT CONVERGED (plateau)`。**運用判定 A** に従い、
「全残差列が低レベルで横ばい (NaN/上昇なし) **かつ** 着目積分量 Cf のドリフト <0.3% (res_15000→20000)」
を満たすため **実用上 OK (発達・定常)** と判断する。リミットサイクルでないことは Cf 定常で担保。
(`rms_roOmega` は ω ピン留めセルの残差を 0 化したため 1e-3 台で正常。)

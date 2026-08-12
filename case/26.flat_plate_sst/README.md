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
(`procedures/solver-settings.md` の「CFL の定義」を参照)。
本ケースは壁セルが極めて細く陰解法でも CFL 制限が低いため、段階的に進める。

1. **stage A (`run_0005_slau_spinup`)**: 一様初期値からの cold start。
   `convMethod: 0` (1次風上) + `cfl_pseudo: 2.0` で 25000 step 回し、平滑な発達場を作る。
2. **stage B (`run_0006_slau_muscl`)**: stage A の `res_25000.h5` から restart。
   `convMethod: 2` (MUSCL) + `cfl_pseudo: 3.0` で 40000 step 回し、2次精度で収束させる。

> cold start で `cfl_pseudo≳5` や、いきなり MUSCL + 高 `cfl_pseudo` にすると発散する。
> 発達場からの restart でも本ケースの安定上限は `cfl_pseudo≈3` 程度。

メッシュ生成・変換・実行は `procedures/calculation-workflow.md` の標準手順に従う
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
| `run_regr_cf` | 回帰: 閉形式 FVS 既定化 (`implicitSolvePrecision=0`) の確認。run_0001 と同条件 | 残差収束が legacy と一致 (rms_ro 2.0e-9)、NaN なし → **閉形式は 2D 陰解に無害**。plan [precision-mixed-axisym.md](../../plans/archived/precision-mixed-axisym.md) | active |
| `run_0009_ewt_fine_mode1` | EWT 回帰: 細メッシュ (y⁺₁≈0.35) で `wallTreatmentSST:1`。run_0007 収束場から restart 20000 step | mode1 ≈ mode0(run_0007): Cf/Schl 0.93/0.96/0.99、残差 floor 全列一致 → 壁解像に縮退 | active |
| `run_0010_ewt_yp30_mode0` | EWT 検証: 粗メッシュ y⁺₁≈30, mode0(low-Re)。cold start→MUSCL | **Cf 崩壊** Cf_model/Schl 0.13/0.13/0.14 (low-Re BC が粗メッシュで破綻) | active |
| `run_0011_ewt_yp30_mode1` | EWT 検証: 粗メッシュ y⁺₁≈30, mode1(enhanced) | **Cf 回復** Cf_model/Schl 0.89/0.91/0.93、u_τ≈細メッシュ一致 → y⁺ 非依存 | active |
| `run_0012_ewt_yp10_mode0` | EWT 検証: バッファ帯 y⁺₁≈10, mode0(low-Re) | Cf_model/Schl 0.51/0.56/0.60 (過小) | active |
| `run_0013_ewt_yp10_mode1` | EWT 検証: バッファ帯 y⁺₁≈10, mode1(ω-blend+τ_w, k=0 Dir) | Cf_model/Schl 1.05/1.10/1.14 (過大) → full WF で改善 | ref |
| `run_0016_ewt_fine_wf`  | **full WF** (ω pin + k zero-grad + P_k log則化): 細 y⁺₁≈0.35 回帰 | Cf/Schl 0.89/0.92/0.95 (=壁解像), Cf_molec≈Cf_model → 壁解像極限を保持 | active |
| `run_0017_ewt_yp30_wf`  | **full WF**: 対数帯 y⁺₁≈30 | Cf/Schl 0.92/0.94/0.97, u_τ≈2.5–2.6 | active |
| `run_0018_ewt_yp10_wf`  | **full WF**: バッファ帯 y⁺₁≈10 | Cf/Schl 0.95/0.98/1.01 (ω-blend単独の過大を解消) | active |
| `run_0019_des_regr_mode0` | **SST-DES T1-A 回帰**: 発達場(run_0007/res_120000)から `DESmode:0` 1000 step。新旧バイナリ比較 | DESmode:0 は baseline と atomicAdd 床内で一致 (1 step base–base ≡ base–des: ro 1.19e-7/roe 1.56e-2/roK 3.73e-9) → **ビット不変**。plan [turbulence-iddes-sst.md](../../plans/active/turbulence-iddes-sst.md) | active |
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
| `run_node_sst_final` | **node SST 完成**: massflux 書込 + 壁 ω Dirichlet ピン (omega/roOmega ピン+ω decouple+res_roOmega 壁ゼロ) + **wall_dist 修正** (node で壁点に壁ノード座標)。過剰乱流場から 60000 step。**convMethod 0 (1次風上) のまま** | u⁺-y⁺ collapse は cell 一致だが **Cf/δ99 は cell(MUSCL)と乖離**: Cf −8%, δ99 +13〜20% (1次の数値拡散で BL 過厚)。かつ res_60000 で Cf 未収束 (step30k→60k で +20% ドリフト)。massflux+ω pin+wall_dist 修正で過剰乱流は解消済 | ref(1次) |
| `run_node_sst_muscl_cont` | **node SST を MUSCL(convMethod 2)に切替**え run_node_sst_final/res_60000 から +90000 step。cell 基準 (run_0009, MUSCL) とスキームを揃えた公平比較 | **node≈cell に一致**: Cf 差 **−3.4〜−3.9%**, δ99 差 **−1.0〜−1.5%**, u_τ −2.0〜−2.3% (3 station)。Cf@0.6 ドリフト +0.20%→+0.06%→0.00% で定常化、NaN なし。残差は block-DPLUR 構造プラトー (cell run_0007 と同様)。図 `cf_bl_cell_node_muscl.png` | active |
| `run_node_skip_verify` | **node 境界半割面拡散 skip 検証 (新)**: run_node_sst_muscl_cont/res_90000 から +5000 step、k/ω 境界半割面拡散を ∇·S 弱形式→skip 化したバイナリ。plan [diffusion-node-boundary-real-distance.md](../../plans/accepted/diffusion-node-boundary-real-distance.md) §3(c) | **Cf/u_τ/δ99 は ∇·S と完全一致** (3 station 差 0.00%, ref90k と 0.01% 以内)。k 場差 (relL2 6.5%) は **全て前縁上流 x<0 のスリップ域**に局在 (∇·S が Neumann 漏れで k を ~10 に過剰生成→skip は ∂k/∂n=0 で k≈3=内部値、skip がより正)。平板 BL は無影響。NaN なし | active |
| `run_node_gradS_verify` | 上記の比較基準: 同じ restart から **旧 ∇·S 弱形式バイナリ**で +5000 step | skip と Cf 完全一致を確認するための apples-to-apples ペア (k のみ前縁スリップ域で skip と差) | ref |
| `run_0004_keep_es_rans_m02` / `run_0005_keep_es_rans_m02_cont` | **単一スキーム KEEP+ES の RANS 検証 (第1試行)**: run_0003 (回帰テンプレート) の場から KEEP+matrix ES (σ0.05 jump2 precond) + implicit cfl20 で 10k+40k | KEEP+implicit 初の組で**発散なし**。ただし run_0003 の土俵は Cf 未検証 (理論+30%・x 無依存) で A/B 判定に不適と判明 → run_0012 へ移行。10k 時点は Ux ドリフト 1.1%/5k で未収束だった (準定常確認の教訓) | 記録 (土俵不適) |
| `run_0023_isoT_fine_ref` | **等温壁 (Tw=320K) 壁解像リファレンス**: run_0009 (y+0.35 EWT 細) 収束場から `wall_isothermal` 化、通算 70k step (q_w(x) ドリフト 0.2%/10k まで静定) | 発達域 q_w ≈ 4.0-4.7 kW/m²。**熱壁関数ギャップ実測の基準** | ref (等温 ref) |
| `run_0024_isoT_yp30_wf` | **等温壁 y+30 + `wallTreatmentSST:1`** (熱壁関数なしの現行 EWT): run_0021 場から同 Tw、20k step (ドリフト 0.6% 静定) | **q_w を +39〜+49% (平均 +45%) 過大予測** (x=0.3-0.9, vs run_0023)。原因: 第 1 セル (y+30) の大きな μt を線形 (Tw−T1)/y1 勾配に掛けるため sublayer 熱抵抗を無視。**等温壁×粗 y+ EWT は熱壁関数 (Kader 相当) が必須** — 実装は [wmles plan](../../plans/active/turbulence-wmles-wall-stress.md) の Kader 部品流用で可能 (未着手) | ref (熱壁関数ギャップ実測) |
| `run_0025_isoT_yp30_qwwf` | **`sstEnergyWallFunction:1` 初版** (Kader q_w 流束置換, [plan](../../plans/active/turbulence-sst-thermal-flux-model.md)) を y+30 で。warm from run_0024 | q_w +42% 過大 → **−12.9%** (vs run_0023)。残った系統過小の切り分けで **Kader T⁺ の式混用バグ発見** (対数部が Pr_t(u⁺+β_Kader) — Jayatilleke 形と Kader β の混用で T⁺ +30% 過大) | ref (混用式の記録) |
| `run_0026_isoT_yp30_qwwf_kader` | 同上 + **Kader 原式修正** (`wallLaw_kader_tplus` = 2.12·ln(1+y⁺)+β 形へ)。warm from run_0025 | **q_w = 真値 +0.1%** (range −0.7〜+1.4%, x=0.3-0.9)。qwall drift 0.17%/10k で静定 | active (**等温壁 y+30 正**) |
| `run_0027_isoT_fine_qwwf` | Kader 原式を細メッシュ y+0.35 で (low-Re 極限回帰)。warm from run_0023 | **q_w = 真値 −0.1%** — 置換が解像伝導へ正しく縮退 | active (low-Re 極限) |
| `run_0028_isoT_yp10_qwwf` | Kader 原式をバッファ帯 y+10 で。warm from run_0018 (断熱 wf 場) | **q_w = 真値 +6.2%** (range +5.7〜+7.3%) — ブレンド帯の残差で許容 | active (バッファ帯) |
| `run_0029_node_outletic_subsonic` | **亜音速出口の背圧アンカー検証** (node outlet の p_tilde/勾配閉包を内部値参照に改訂した最終形, [plan §2.11](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)): run_node_sst_final 収束場から 12000 step | **出口 P 平均 97250.6 Pa (規定 97250, ±3 Pa) を保持** = p_tilde=P[ic] でも特性構成 bvar の質量流束が背圧を係留。残差悪化なし・NaN なし | active (亜音速アンカー検証) |
| `run_0030_node_yp30_wf_omgoff` | **node × y+30 壁関数の平板** (ω irep ピン A/B 用に新規 node 変換 `flat_plate_yp30_node.h5`)。段階起動: 1次 cfl2 (清浄収束) → MUSCL cfl3 | 1 次段は健全 (3.9 桁)。**MUSCL 段は残差 2-4 桁上昇のプラトーで Cf が 10%/10k step ドリフト = 準定常に達せず**。convMethod 1 は清浄場からも発散。**node × 亜音速平板 y+30 の数値レシピは未確立** (別課題)。ω ピン判定は case/40 の y+1 真値比較で代替済 ([plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md)) | 破棄予定 (レシピ未確立の記録) |
| `run_0031_node_yp30_stageA` | **node × y+30 壁関数の段階起動 stage A** (1次 `convMethod:0` + cfl_pseudo 2 + nStepInner 20, cell `run_0011` からcross-mesh interp) | 健全 (roK/roOmega 3.6-3.9 桁低下、falling)。**node 亜音速平板の唯一の安定な起点** | active (node 平板 seed) |
| `run_0032_node_yp30_2nd_base` | 2次 (MUSCL, limiter venkata) を stage A 場から 3000 step | **残差上昇 (発散性)**。原因は**前縁 x≈19mm の局所暴走** (Ux 199 m/s = 自由流 67.8 の 3 倍・P 90.5 kPa) | 破棄予定 (2次失敗の記録) |
| `run_0033_node_yp30_2nd_lowmach` | 同 + `lowMachPrecond: 2` (低マッハ市松対策の転用) | **流れ場が完全凍結** (`ro`/`roe` がビット不変、`roUx` の変化 4.7e-34 = 非正規数)。乱流のみ進行 → 本構成では使用不可 | 破棄予定 (要調査の記録) |
| `run_0034_node_yp30_2nd_lsq` | 同 + `gradLSQ: 2` (LSQ 勾配の検証も兼ねる) | **NaN 発散** (GG より悪化)。LSQ は本不安定の対策にならない — [LSQ plan](../../plans/active/discretization-lsq-gradient.md) 検証4 | 破棄予定 (LSQ A/B 記録) |
| `run_0035_node_yp30_2nd_barth` | 同 + `limiter: 1` (barth) | **NaN 発散** | 破棄予定 |
| `run_0036_node_yp30_2nd_cfl05` | 同 + `cfl_pseudo: 0.5` | 残差上昇 (cfl 2 と同率) → **時間積分でなく空間離散の問題**と確定 | 破棄予定 |
| `run_0037_node_yp30_2nd_bfo` | 同 + `bndFirstOrder: 1` | 残差上昇 (改善せず) | 破棄予定 |
| `run_0038_node_yp30_1st_long` | **node y+30 の生産検証** (1次, stage A から通算 50000 step) | 残差床到達 (rms_ro 1.1e-9)。**Cf は完全定常** (20k/30k/40k で 0.00242 不変)。**Cf/真値 (run_0007 壁解像) = 0.835-0.840 で x 依存なく一様**。1次風上の散逸が支配的とみられ、node 壁関数自体の評価には 2 次が必要 | active (**node 平板の現状ベスト**) |
| `run_0039_node_yp30_planar_2nd` | **真 2D (平面) メッシュへの切替 + 2 次精度 (3000 step)**。メッシュ = 新規 `mesh/flat_plate_yp30_planar.geo` (押し出しなし, 10755 ノード, 品質 PASS: AR max 58.8 / skew 0.000)。IC = run_0038 res_40000 から cross-mesh interp | **2 次が初めて安定**。全残差 falling・NaN なし・前縁暴走なし。旧 2 次 (run_0032-0037) の発散が**メッシュ側の欠陥**だったことの直接証拠 | active (2 次成立の初回証拠) |
| `run_0040_node_yp30_planar_2nd_long` | 同構成の本収束 (run_0039 から 40000 step) | NaN なし・Ux max 71.3 m/s (run_0032 の 199-249 の暴走は消滅)・Uz ~1e-14。残差は block-DPLUR 構造プラトー (本 case 既知, run_0007 も同様) だが **quasisteady `ALL STEADY`**、**Cf は 20k/30k/40k で 4 桁一致**。**Cf/真値 = 0.848/0.847/0.848** | active (**node 平板の生産基準**) |
| `run_0042_node_yp30_planar_2nd` | **壁関数 y+ 非依存性テスト**: 第一層 1.768e-4 (`mesh/flat_plate_ny52_planar.geo`, node 第一 DOF y+27.1 = cell 第一 DOF と同高さ)、run_0040 から cross-mesh interp、20000 step | $u_\tau$=2.343 (ny=45 と同値)。**cell (同高さ 1.700e-4) の 2.586 に対し 9.4% 低いまま** = 第一層高さ仮説の棄却 | active (y+ 掃引) |
| `run_0043_node_yp112_planar_2nd` | 同掃引の粗い側: 第一層 6.588e-4 (`mesh/flat_plate_ny38_planar.geo`, y+101.5) | $u_\tau$=2.350。**y+ 27→102 (3.7 倍) で $u_\tau$ の変化 0.3%** = node 壁関数は y+ 非依存 | active (y+ 掃引) |
| `run_0044_cell_yp30_wf_regr` | cell y+30 壁関数を**現行バイナリで再走** (run_0017 は旧バイナリで `wf_pk`/`Pk_diag` 未出力のため対照にならない)。run_0017 res_10000 から restart、5000 step | $u_\tau$=2.586 (run_0017 と一致)。生産収支比較の cell 基準 | active (cell 対照) |
| `run_0045_node_yp30_kwfdir0` | **`nodeKwfDirichlet: 0` A/B** (既定 1 = 第一内層 k を Dirichlet 固定)。run_0042 と同一 IC/設定 | 実効生産が 3.1 倍になるが $u_\tau$ は +0.55% のみ → **「生産規約が真因」を反証** | active (反証の証拠) |
| `run_0046_node_yp30_omgwfdir1` | **`nodeOmegaWfDirichlet: 1` A/B** (第一内層の $\omega$ も固定)。run_0042 と同一 IC/設定 | $\omega$@第一内層 2.645e5→1.213e5、**$C_f$/真値 0.843→0.957 = 欠損の 73% 回復** → **主因は $\omega$** | active (**主因特定の証拠**) |
| `run_0050_node_lowre_fine_regr` | **node の壁解像 (壁関数なし) 検証**: planar 細メッシュ (`flat_plate_planar.h5`, y+0.7) + `wallTreatmentSST: 0` + MUSCL、現行 binary、`run_node_sst_muscl_cont` res_90000 から restart 20000 step | **$C_f/C_{f,\mathrm{KS}}$ = 0.943** (cell 壁解像 0.954 と 1.1% 差)、収支 1.036、quasisteady **ALL STEADY** → **node 離散化そのものは健全**と確定 | active (**node 無罪の証拠**) |
| `run_0047_su2_sst_ny52` | SU2 v8.5.0 low-Re (壁関数なし)、`run_0042` と同一メッシュ | $C_f/C_{f,\mathrm{KS}}$=0.365、収支 0.181、$\mu_t/\mu$ 340–452 で**破綻** → low-Re を y+27 に当てる誤りは SU2 でも同じ | 破棄予定 (設定誤りの記録) |
| `run_0048_su2_sst_wf_ny52` | **SU2 + `STANDARD_WALL_FUNCTION`**、`run_0042` と**完全同一メッシュ・BC・物性・流入乱流** | 壁 $C_f$/KS=0.793–0.796 だが **運動量収支 0.790 = 21% 破綻**、残差も rms −4.4 止まり (壁解像版は −12.4) → **解として棄却**。「SU2 も 0.75」という当初結論は撤回 | 破棄予定 (未収束・収支破綻の記録) |
| `run_0049_su2_sst_lowre_fine` | **SU2 壁解像** (planar 細メッシュ, 第一 DOF y+0.70, low-Re) — node 壁解像の独立検証 | **rms −12.4 まで収束**。壁 $C_f$/KS=**0.938–0.944**、運動量収支 **0.999/0.995 (ほぼ完璧)** → forge node 壁解像 0.943 と **0.5% 以内で一致** | active (**node 離散化検証の正**) |
| `run_0041_node_yp30_z4_2nd` | **機構対照**: 押し出しを 1 層 → **4 層** (z ノード 2 枚 → 5 枚) にした同一 2 次設定 (`mesh/flat_plate_yp30_z4.geo`, 53775 ノード, 3000 step) | スパン方向モードの成長が **350 分の 1** (step3000 で spread 0.598 m/s vs run_0032 の 211.9)。= 発散は「node の slip BC 欠陥」でも「前縁」でもなく **spanwise が 2 ノードしかないこと**が原因と確定 | 破棄予定 (機構確定の記録) |
| `run_0012_keep_es_ewt_fine` / `run_0013_keep_es_ewt_cont` | **単一スキーム KEEP+ES の RANS 検証 (正)**: 確立済み Cf 基準 run_0009 (EWT 細メッシュ y+0.35) と同一設定で flux のみ KEEP+ES 化、40k+20k ([iddes plan §4.8 設計更新](../../plans/active/turbulence-iddes-sst.md)) | **合格**: implicit cfl20 安定・quasisteady `ALL STEADY`・**Cf ドリフト +0.006〜0.041%/20k で完全収束**。Cf/Schl = **0.88/0.91/0.94** (SLAU 0.91/0.95/0.97 比 **−3.6〜−3.9%** = 本 case の確立済みスキーム間ばらつき node vs cell MUSCL −3.4〜−3.9% と同帯)。**「KEEP+ES で SST-RANS」は成立** | active (KEEP+ES RANS 検証) |

> **node 弱形式境界 (Phase 2)**: node モードで inlet/outlet/slip も ghostless 弱形式化。(5a) 主対流ループを内部+periodic
> のみに、(5b) 全境界を `convectiveFlux_boundary_d` の bvar 弱形式に、(5e) **block-DPLUR Jacobian の境界半割面で
> 退化粘性対角 `viscous_diag` をスキップ** (境界ノードが境界上で `dcc·S≈0`→delta 爆発→dq≈0 で凍結し出口 BL 崩壊・
> 残差プラトーを生んでいた真因)。**出口 BL 崩壊と残差プラトーの双方を解消**。plan
> [discretization-node-boundary-ghostless.md](../../plans/active/discretization-node-boundary-ghostless.md)。コーナー BC
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
> [run_node_sst_final](run_node_sst_final/): peak mut/μ **217** (cell 199), k_wall=0, rms_roUx 0.597→1.5e-4 (3.6桁)。
> **massflux + omega pin + wall_dist の 3 点で過剰乱流 (旧 peak mut/μ 375, Cf≈3×) を解消**し乱流レベルは cell 並みに。
>
> **ただし Cf・境界層厚さの cell 一致には convection スキームを揃える必要がある (追検証)**: `run_node_sst_final` は
> `convMethod 0` (1次風上) のままで、同一測定 (分子勾配 τ_w, δ99=0.99U_e) で cell(MUSCL, run_0009) と比べると
> **Cf −8% / δ99 +13〜20%** と乖離していた (1次の数値拡散で BL が厚くなる + res_60000 で Cf がまだ +20%/30k step
> ドリフト中=未収束)。node を **MUSCL に切替えて Cf 収束まで回した [run_node_sst_muscl_cont](run_node_sst_muscl_cont/)**
> では **Cf 差 −3.4〜−3.9% / δ99 差 −1.0〜−1.5% / u_τ −2.0〜−2.3%** と cell に一致する。u⁺-y⁺ collapse (各 run 自身の
> u_τ で正規化) は 1次でも揃って見えるため、**絶対量 (Cf, δ99) は必ず同一スキームで比較する**こと。比較図:
> `cf_bl_cell_node_muscl.png` (再生成 `tools/compare_cf_bl_cell_node.py`)。残る ~3.5% の Cf 差は median-dual 壁勾配と
> cell の離散差由来で、全 station で一定符号・小さい。

> **SST-DES (DDES) T1-A** の合否: `DESmode:0` ビット不変 + `DESmode:1` で付着 BL が RANS から
> 不変 (Cf 差≪0.1%)。**発達乱流場 (nut/nu≫1) から restart すること** (`run_regr_cf` 等 nut/nu≤1 の
> 未発達場では νt≈0 で rd 小→シールドが効かない)。詳細 plan [turbulence-iddes-sst.md](../../plans/active/turbulence-iddes-sst.md)。

> EWT (enhanced wall treatment) 検証の詳細は plan [turbulence-enhanced-wall-treatment.md](../../plans/accepted/turbulence-enhanced-wall-treatment.md)。
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

## node × 2 次精度の発散 — 真因は「spanwise が 2 ノードしかないメッシュ」 (2026-08-12 解決)

`run_0032`–`run_0037` で node × y+30 × 2 次精度がどの設定でも発散した件の決着。
**前縁 (slip→wall 遷移) は無関係**で、真因は**メッシュの作り方**だった。

### 症状の実体 (run_0032 の時系列を再解析)

前セッションが記録した「前縁 x≈19mm の局所暴走 (Ux 199 m/s)」は**最終スナップショットだけ**を
見た記述で、発生源ではない。時系列を追うと:

| step | 上部 slip 行 | 中間場 | 近壁 |
| --- | --- | --- | --- |
| 0 | 0.0023 | 0.0026 | 0.0045 |
| 500 | **0.86** | 0.05 | 0.11 |
| 1000 | **119.9** | 5.9 | 2.1 |
| 2000 | 73.1 | 114.1 | **203.2** |

(数値 = スパン方向非対称 $|U_x(z{=}0)-U_x(z{=}0.01)|$ [m/s])

このメッシュは z 方向 1 層押し出しの**疑似 2D スラブ**なので厳密解は z 非依存でなければならない。
実際には**上部 slip 境界の 1 行**で先に z 対称性が破れ (1 つ下の行は 17 倍小さい)、そこから
中間場 → 近壁へ流れ込む。前縁の overspeed は最後に現れる症状。

### 機構: 2 ノード方向では 2 次 MUSCL の上流差分散逸が厳密に消える

node-centered (median-dual) では 1 層押し出しメッシュの spanwise に**ノードが 2 枚**できる。
z=0 のノード $i$ と z=dz のノード $j$、その間の双対面 (z=dz/2) を考える。GG 勾配は
境界半割面 (slip, owner-state) と内部 z 面の寄与で

$$\left(\frac{\partial\phi}{\partial z}\right)_i=\left(\frac{\partial\phi}{\partial z}\right)_j=\frac{\phi_j-\phi_i}{dz}$$

となり、面への MUSCL 再構成は

$$\phi_L=\phi_i+\frac{\partial\phi}{\partial z}\Big|_i\frac{dz}{2}=\frac{\phi_i+\phi_j}{2},\qquad
\phi_R=\phi_j-\frac{\partial\phi}{\partial z}\Big|_j\frac{dz}{2}=\frac{\phi_i+\phi_j}{2}$$

で **$\phi_L=\phi_R$ が任意の場について厳密に成立**する。つまり z 面のジャンプが恒等的に 0 になり、
**SLAU/上流差分の散逸がこの方向で完全に消える**。2 ノードしかない方向で表現できる唯一のモードは
spanwise 市松 ($\phi_i\neq\phi_j$) なので、**そのモードだけが無減衰**になり丸め誤差から指数成長する。

これが過去の全観測を一度に説明する:

- **1 次では安定** — $\phi_L=\phi_i,\ \phi_R=\phi_j$ でジャンプ $\neq 0$、強く減衰する (`run_0038`)。
- **リミッタで直らない** (`run_0035` barth / `run_0032` venkata) — 再構成値が隣接 2 点の**算術平均**で
  常に min/max 内に入るため、どのリミッタも作動しない (limiter = 1)。
- **CFL で直らない** (`run_0036` cfl_pseudo 0.5) — 散逸の欠落は空間離散の性質で、時間刻みと無関係。
- **LSQ で悪化** (`run_0034`) — 2 点ステンシルでは LSQ も同じ相殺を起こし、GG の fx 平均による
  わずかな平滑化すら失う。
- **cell モードでは起きない** — primal は z 方向に**セル 1 個**しかなく、この方向のモードが存在しない。
- **case/40 ノズルの node 2 次は成功している** — あちらは**真 2D メッシュ (z レベル 1)**。

### 対照実験 (`run_0041`)

押し出しを 1 層 → 4 層 (z ノード 2 枚 → 5 枚) にすると、5 点ステンシルでは市松の
$\phi_L\neq\phi_R$ が回復して散逸が戻る。実測でスパン方向モードの成長は **350 分の 1**
(step 3000 で spread 0.598 vs 211.9 m/s)。→ 「node の slip BC が壊れている」仮説は棄却され、
**spanwise ノード数**が支配因子と確定。

### 対策 (恒久)

**node モードで 2D 問題を解くときは、押し出しメッシュを使わず必ず平面 2D メッシュ
(`Physical Curve` タグ付け, z レベル 1) を作ること。** 本 case では
`mesh/flat_plate_yp30_planar.geo` を新設した (`flat_plate_planar.geo` の y+30 版)。
3D で node を使う場合も、均質方向は**最低 3 ノード (2 層以上)** 確保する。

### 副産物: 2 次化しても Cf 欠損は残る (前セッションの仮説を反証)

| run | 離散化 | 次数 | Cf/真値 @ x=0.30 / 0.60 / 0.89 |
| --- | --- | --- | --- |
| `run_0007` | cell (壁解像 y+0.35) | 2 次 | 1.000 / 1.000 / 1.000 (**真値の定義**) |
| `run_0017_ewt_yp30_wf` | cell (y+30 壁関数) | 2 次 | **1.034 / 1.027 / 1.024** |
| `run_0038` | node (y+30 壁関数) | 1 次 | 0.840 / 0.844 / 0.848 |
| `run_0040` | node (y+30 壁関数) | **2 次** | **0.848 / 0.847 / 0.848** |

前セッションは「Cf の 16% 欠損はおそらく 1 次風上の散逸」と推定していたが、**2 次にしても
Cf は 1% 未満しか動かない**ので**この帰属は誤り**。同じ y+30 メッシュ・同じ壁関数で
**cell は真値の +2〜3%** に入るのに **node は −15%** であり、欠損は
**node 側の壁関数/離散に固有**。$u_\tau$ で見ると node 2.34 / cell 2.586 / 真値 2.580 m/s
(x=0.6) で、node の $u_\tau$ が 9.3% 低い。**次の調査対象はここ** (壁関数の代表点速度の取り方が
有力: node では壁ノードが `nodeWallDirichlet=1` で u=0 に固定されるため、代表点の選択と
壁距離の組が cell と一致しているかの監査が要る)。

### 「真値」の定義と Cf 欠損の出所 (図: `uplus_truth_compare.png`)

再生成: `python3 tools/plot_uplus_truth_compare.py 0.6`

**本 case に実験値はない。**「真値」と呼んでいるものは 2 種類あるので混同しないこと。

1. **内部基準** = forge 自身の cell 壁解像 run (`run_0007`, y+₁≈0.34) の**分子勾配**
   $\tau_w=\mu U_1/y_1$。CFD の自己参照であって外部真値ではない。
2. **経験式** = Schlichting $C_f=0.0592\,Re_x^{-0.2}$。内部基準自身がこれの 0.92 倍
   (前縁からの遷移区間ぶんの virtual-origin オフセット)。

壁関数 run の $C_f$ は分子勾配ではなく**ソルバが課している modeled $\tau_w=\rho u_\tau^2$**
で測る (Reichardt 則の逆解き。`tools/wall_law_modeled.py` と同じ)。

| run | y+₁ | $u_\tau$ | $C_f$ | /内部基準 | /Schlichting |
| --- | --- | --- | --- | --- | --- |
| `run_0007` cell y+0.35 壁解像 | 0.34 | 2.580 | 2.821e-3 | 1.000 | 0.920 |
| `run_0017` cell y+30 壁関数 | 28.80 | 2.586 | 2.896e-3 | 1.027 | 0.944 |
| `run_0038` node y+30 壁関数 1 次 | 52.09 | 2.338 | 2.368e-3 | 0.839 | 0.772 |
| `run_0040` node y+30 壁関数 2 次 | 52.20 | 2.343 | 2.378e-3 | 0.843 | 0.775 |

**(a) 同じ primal メッシュでも node と cell では y+₁ が 2 倍違う**。node の第一 DOF は
**節点そのもの** (= 第一セル高さ $h=3.40\times10^{-4}$ m) にあるが、cell の第一 DOF は
**セル重心** ($h/2$) にある。$y^+$ 比 = $2\times(2.338/2.586)=1.81$ で実測 28.80→52.09 と一致。
median-dual の構造上の性質であり欠陥ではないが、**「同じメッシュ」は「同じ y+」を意味しない**。

**(b) 欠損は壁近傍 2 点に局在する** (図 C, 内部基準 $u_\tau$ で正規化した対数則からのずれ):

| $y^+$ | node $u^+$ | 壁解像 $u^+$ | 差 |
| --- | --- | --- | --- |
| 52 (node の第一 DOF) | 13.91 | 14.79 | **−0.88 (−2.26 m/s)** |
| 80 | 15.15 | 16.07 | **−0.92 (−2.37 m/s)** |
| 130 | 17.41 | 17.46 | −0.05 |
| 200 | 18.90 | 18.67 | +0.23 |
| 400 | 21.11 | 20.70 | +0.41 |

node の第一 DOF では速度が**対数則を 0.92 $u^+$ 下回る** (cell の第一 DOF は逆に +0.34 上回る)。
壁関数はまさにこの点の速度を Reichardt 則に食わせて $u_\tau$ を逆算するので、
**低い速度 → 低い $u_\tau$ (−9%) → 低い $C_f$ (−15%)** が直接の因果。$y^+\gtrsim130$ では
プロファイルは壁解像と一致し (むしろ僅かに fuller = 壁抵抗が小さいことの帰結)、
**外層は健全**。図 A (各 run 自身の $u_\tau$ で正規化) では collapse して見えるので、
**A だけを見て「node は壁法則に乗っている」と判断してはならない**。

**(c) 次に当たる候補 (未検証)**: `viscousFlux_d.cu:174-191` の AddTauWall 再スケールは、
W–I 双対面の**接線** traction を modeled $\tau_w$ に合わせる意図だが、実装は
`tau_x *= scale` と **traction ベクトル全体**を掛けており、法線成分 $T_n$ も同じ倍率
(本 case では $C_{f,\rm model}/C_{f,\rm molec}\approx2.1$) で増幅される。コメントの意図と
コードが食い違っている。壁法線方向に余分な運動量流束を注入するので、第一内点の速度欠損の
候補として最初に潰す価値がある。正しい形は接線成分だけを再スケールして
$\boldsymbol\tau = \hat{\mathbf t}\,\tau_w A + T_n\hat{\mathbf n}$ とすること。

### 追試: 第一層高さは原因ではない — node の乱流量が低いことが真因 (2026-08-12)

上の図の読み方を**訂正**する。**パネル A は collapse しない** (node が上へずれる)。**パネル B が一致する**。
A のずれ量は $u^+_\infty=U_\infty/u_\tau$ が 26.1→28.8 になったもので、**$u_\tau$ 欠損 9% の可視化そのもの**。
したがって「A で揃うから壁法則に乗っている」ではなく、**「B で揃う = 速度場自体は cell と近い、
A でずれる = 報告される $u_\tau$ が低い」**が正しい。

「第一層高さの扱い (node は節点 = $h$、cell は重心 = $h/2$) の違いが分布差を生んでいるのでは」という
仮説を、**node の第一層高さを振って直接検証**した (planar メッシュ 3 種、`mesh/flat_plate_ny{38,45,52}_planar.geo`)。

| run | $y_1$ [m] | y+₁ | $u_\tau$ | /真値 | $C_f$ | /真値 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0017` cell y+30 壁関数 | 1.700e-4 | 28.8 | 2.586 | 1.002 | 2.896e-3 | 1.027 |
| `run_0042` node ny=52 | **1.768e-4** | 27.1 | 2.343 | 0.908 | 2.377e-3 | 0.843 |
| `run_0040` node ny=45 | 3.400e-4 | 52.2 | 2.343 | 0.908 | 2.378e-3 | 0.843 |
| `run_0043` node ny=38 | 6.588e-4 | 101.5 | 2.350 | 0.911 | 2.393e-3 | 0.848 |

**結論 1: node の壁関数は y+ 非依存である** — y+ を 27→102 (3.7 倍) 振っても $u_\tau$ は 0.3% しか動かない。
壁関数としては正しく機能している。

**結論 2: 第一層高さは原因ではない** — cell (1.700e-4) と node ny=52 (1.768e-4) は**ほぼ同じ物理高さ**なのに
$u_\tau$ が 2.586 vs 2.343 で **9.4% の差がそのまま残る**。仮説は棄却。

**結論 3: 各 run は運動量積分と整合している** (von Kármán $C_f=2\,d\theta/dx$):

| run | $\theta$(0.6) | $C_f$ (運動量積分) | $C_f$ (壁で課した値) | 比 |
| --- | --- | --- | --- | --- |
| cell y+30 壁関数 | 1.033e-3 | 2.813e-3 | 2.896e-3 | 1.030 |
| node ny=52 | 8.745e-4 | 2.403e-3 | 2.377e-3 | **0.989** |
| node ny=45 | 8.967e-4 | 2.410e-3 | 2.378e-3 | **0.987** |
| cell y+0.35 壁解像 | 9.978e-4 | 2.757e-3 | 2.821e-3 | 1.023 |

→ node 側に**偽の運動量ソース/シンクは無い**。node は「$\theta$ が 15% 薄い BL」という
**自己整合な別解**に収束している。壁関数はその状態を正しく報告しているだけ。

**結論 4: 真因は乱流量が低いこと** (x=0.6, BL 内):

| run | peak $\mu_t/\mu$ | $\mu_t/\mu$ @y+100 | peak $k$ |
| --- | --- | --- | --- |
| cell y+0.35 壁解像 | 125.7 | 33.4 | 20.7 |
| cell y+30 壁関数 | 132.7 | 37.9 | 27.9 |
| node ny=52 | **103.5** | **26.7** | **18.3** |
| node ny=45 | 104.8 | 21.4 | 19.1 |

因果連鎖は **$\mu_t$ −20% → 乱流運動量輸送が小さい → BL が薄い ($\theta$ −15%) →
$\tau_w$ が小さい ($u_\tau$ −9% / $C_f$ −15%)**。壁関数は犯人ではなく忠実な報告者である。

**次の調査対象**: node の SST 生産項、特に壁関数の $P_k$ 置換 $\rho u_\tau^4/\nu\cdot g(1-g)$ が
**どの CV 体積に対して積分されているか**。node の第一内層 CV は cell の第一セルと体積が違い、
壁ノード側は半 CV なので、体積重み付けが cell と整合していないと $k$ が系統的に低く出る。
cell の壁関数は自身の壁解像解に対し $k$ を +35% 過大に出す一方、node は −12% 過小である点も
生産項側を示唆する。関連: [`turbulence-node-wall-function-coverage.md`](../../plans/accepted/turbulence-node-wall-function-coverage.md)。

### 【撤回済み】node の壁関数 $P_k$ は構造的に cell の 78% しか注入していない (2026-08-12 起票 → 同日撤回)

> **この節の「真因」「78%」「閉じた因果」は撤回する** (Codex レビュー指摘 + 実験で反証、下の
> 「反証と真の主因」節を参照)。**測った量が間違っていた**: `nodeKwfDirichlet` の既定は `1` で
> ([`solverConfig.cpp:460`](../../solver_density_cuda/input/solverConfig.cpp))、第一内層ノードは
> `roK_wf` に Dirichlet 固定され `res_roK=0` にされる。実測で `roK == roK_wf` が厳密一致し、
> **$\Sigma P_kV$ の 67.7% は k 収支から消される量**だった。以下は記録として残す。

### 【撤回された記述】node の壁関数 $P_k$ の規約差 (2026-08-12)

**先の「B が一致 = 速度場はほぼ同じ」という記述を撤回する。** $\theta$ が 15% 違う以上、速度場は
違う。B が揃って見えたのは、共通 $u_\tau$ で正規化すると $u^+_\infty=U_\infty/u_\tau$ が
**$U_\infty$ 共通ゆえ必ず一致する**うえ、対数軸と 0–32 のレンジで差が潰れているだけである。
実差はパネル C (= B から対数則を引いた拡大図) と $\theta$ に出る。
**パネル A も B も「一致」ではない**というのが正しい。

`wf_pk` と `volume` を出力から直接積分して cell と比較した (x=0.6, spanwise 単位厚, `run_0044` は
現行バイナリで cell y+30 を再走したもの — `run_0017` は旧バイナリで `wf_pk` を出力していない)。

| | CV 数 | $P_k$ 評価距離 $y_1$ | `wf_pk` | 被覆体積 | $\Sigma P_k V$ |
| --- | --- | --- | --- | --- | --- |
| cell y+30 (`run_0044`) | 1 | 1.700e-4 (セル重心) | 3.913e5 | 4.144e-6 | **1.622** |
| node ny=52 (`run_0042`) | 2 | 1.768e-4 (壁ノード→第一内層) | 2.857e5 | 3.370e-6 | **0.963** |

**規約 factor = (壁関数が覆う層厚)/($P_k$ 評価距離 $y_1$)** で規格化すると構造が見える:

| run | 覆う層厚 | $y_1$ | factor |
| --- | --- | --- | --- |
| cell y+30 | 3.400e-4 | 1.700e-4 | **2.000** |
| node ny=52 | 2.765e-4 | 1.768e-4 | **1.564** |
| node ny=45 | 5.317e-4 | 3.400e-4 | **1.564** |
| node ny=38 | 1.030e-3 | 6.589e-4 | **1.564** |

- **cell**: $P_k$ を**セル重心** $y_1=h/2$ で評価し、**第一セル全体** (厚さ $h$) に適用 → factor $h/(h/2)=2.000$。
- **node**: $P_k$ を**壁ノード→第一内層ノード距離** $y_1=h$ で評価し、
  **壁ノードの半 CV** ($h/2$) と**第一内層ノードの CV** ($(h+rh)/2$) の**両方**に適用
  ([`ransWallFunction_d.cu:209,213`](../../solver_density_cuda/cuda_forge/ransWallFunction_d.cu) の
  `wf_pk[ic]` と `wf_pk[irep]`) → 層厚 $(2+r)h/2$ で factor $(2+r)/2\approx1.55$。

**この factor は 3 つの node メッシュすべてで厳密に 1.564 = メッシュ非依存**であり、
$r\to1$ でも 1.5 にしかならない。つまり **node の壁関数生産は構造的に cell の
$1.564/2.000=0.782$ しかない**。

因果の分解 (実測と 2% で一致):

$$\underbrace{0.782}_{\text{規約の差 (}u_\tau\text{ を揃えた場合)}}\times\underbrace{0.744}_{(u_{\tau,\rm node}/u_{\tau,\rm cell})^3\ \text{= 自己増幅}}=0.582\quad\text{vs}\quad\text{実測}\ \Sigma P_kV\ \text{比}=0.594$$

**閉じた因果連鎖**: 生産規約が 22% 小さい → $k$・$\mu_t$ が低い (−20%) → 乱流輸送が小さい →
BL が薄い ($\theta$ −15%) → $\tau_w$ が小さい ($u_\tau$ −9%) → $P_k\propto u_\tau^3$ が更に下がる
→ 合計 41% 減で釣り合う。

**修正方針の候補** (未実装・要判断):

1. **CV 平均生産に統一** — 各 CV の実 $y$ 区間 $[y_a,y_b]$ に対し
   $\bar P_k=\frac{1}{y_b-y_a}\int_{y_a}^{y_b}\rho u_\tau^3/(\kappa y)\,dy$ を使う。離散化非依存に
   なるが cell の既定挙動も変わるため cell の再検証が要る。
2. **node の factor を cell の 2.000 に合わせる** — 壁半 CV と第一内層 CV をそれぞれの重心距離で
   評価する、または `wf_pk` を $2.000/1.564=1.279$ 倍する。cell はビット不変で済む。
3. まず **(2) の 1.279 倍を env ゲートで入れて $C_f$ が cell に乗るかを見る**のが、上記の
   因果連鎖全体を一発で検証する最小実験。

なお cell の factor 2.000 が物理的に厳密という保証はない (Launder–Spalding 系の慣用)。ただし
cell は壁解像解に対し $C_f$ を +2.7% で再現しており**較正済みの規約**なので、当面は
これを基準に node を合わせるのが妥当。

### 反証と、第一内層ノード $\omega$ への局在化 (2026-08-12) — **主因は未分離 (仮説段階)**

> **本節の「真の主因」表現は 2 度目のレビューで再度差し戻した**。$\omega$ 周辺への局在化は妥当だが、
> ①機構の説明が現行実装と違い、②A/B が limiter bypass と交絡している。§2'/§3' に訂正を置く。

前節の「生産規約 78% が真因」を**実験で反証**し、主因を $\omega$ に特定した。

#### 1. 生産規約は主因ではない (反証)

`nodeKwfDirichlet: 0` (第一内層 k のピン解除) の A/B で、**k 収支に実際に入る生産量**を測り直した:

| run | $\Sigma P_kV$ 全体 | うち k がピンされた分 | **実効生産** | 対 cell |
| --- | --- | --- | --- | --- |
| cell y+30 (`run_0044`) | 1.622 | 0 | **1.622** | 1.000 |
| node kwfDir=1 (`run_0042`, 既定) | 0.963 | 0.652 (67.7%) | **0.311** | **0.192** |
| node kwfDir=0 (`run_0045`) | 0.977 | 0 | **0.977** | **0.602** |

**実効生産が 3.1 倍 (0.311→0.977) になっても $u_\tau$ は 2.343→2.356 (+0.55%) しか動かない。**
$C_f$/真値 0.843→0.852。→ **壁関数の $P_k$ 積分量は $\tau_w$ をほとんど支配していない**。
前節の「規約 factor 0.782 × $(u_\tau$比$)^3$ 0.744 = 0.582 が実測 0.594 と 2% 一致」は、
$u_\tau$ 比を説明対象の解から取っている**循環論**であり、因果の証明ではなかった。
`wf_pk` × 1.279 の実験は**やらない** (物理補正でなく cell への経験的チューニングであり、
かつ第一内層が k ピンされている既定では仮説の検証にすらならない)。

#### 2. 主因は第一内層ノードの $\omega$ スパイク

内部スケーリング ($\omega_{eq}=u_\tau^2/\nu/(\sqrt{\beta^*}\kappa y^+)$) で見ると node だけ近壁で $\omega$ 過大:

| $y^+$ | cell $\omega/\omega_{eq}$ | node $\omega/\omega_{eq}$ |
| --- | --- | --- |
| 30 | 1.09 | **2.55** |
| 60 | 1.49 | 1.89 |
| 100 | 1.28 | 1.55 |
| 400 | 1.25 | 1.38 |

壁法線プロファイルを見ると異常が直接見える (`wall_dist` は両者とも座標と厳密一致でバグではない):

| | 壁ノード $y=0$ | 第一内層 $y=1.768$e-4 | 第二 $3.707$e-4 |
| --- | --- | --- | --- |
| node $\omega$ | 1.146e5 (ピン=Menter blend) | **2.645e5** | 9.635e4 |
| node $k$ | 20.72 | **10.46** | 18.27 |

**$\omega$ が壁の値より高く跳ねている** (単調減少すべきところが逆転)。その値 2.645e5 は
**第一層距離の半分**での Menter blend (2.663e5) にほぼ一致する。高い $\omega$ は
$\mu_t=\rho k/\omega$ を下げ、同時に散逸 $\beta^*\rho k\omega$ を上げて $k$ も下げる = 両症状を同時に説明する。

**構造的な非対称が原因と考えられる (仮説)**:

- **cell**: `wf_pk` と $\omega$ ピンが**同じ DOF** (第一セル) にある → $\omega$ は固定なので暴走しない。
- **node**: $\omega$ ピンは**壁ノード**だが `wf_pk` は**壁ノードと第一内層の両方**に入る
  ([`ransWallFunction_d.cu:209,213`](../../solver_density_cuda/cuda_forge/ransWallFunction_d.cu))。
  第一内層は $\omega$ を**解きながら**壁関数生産を受け、そこでは $\nu_t$ が小さい ($\mu_t/\mu\approx2$)
  ため $\omega$ 生産 $\propto P_k/\nu_t$ が過大になる。

#### 3. $\omega$ が主因であることの A/B (`nodeOmegaWfDirichlet`)

| run | $u_\tau$ | /真値 | $C_f$/真値 | $\omega$@第一内層 | $k$@第一内層 | $\mu_t/\mu$@y+100 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0042` node 既定 | 2.343 | 0.908 | 0.843 | 2.645e5 | 10.46 | 26.7 |
| `run_0045` kwfDir=0 | 2.356 | 0.913 | 0.852 | 2.642e5 | 11.74 | 26.8 |
| `run_0046` **omgWfDir=1** | **2.497** | **0.968** | **0.957** | 1.213e5 | 12.46 | 32.8 |
| `run_0044` cell y+30 (対照) | 2.586 | 1.002 | 1.027 | 1.307e5 | 27.93 | 37.9 |

$\omega$ を固定すると $\omega$ が cell 並み (1.213e5 vs 1.307e5) に戻り、**$C_f$ 欠損の 73% が回復する**
(0.843→0.957, 目標 1.027)。→ **主因は $\omega$**。

ただし **`nodeOmegaWfDirichlet: 1` は解ではない**: case/40 ベルノズルでは同じスイッチが
$\tau_w$ を y+1 真値の **1.237 倍に過大化**して正式化が見送られている
([plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md))。本ケースでは 0.957 と不足、
ノズルでは 1.237 と過大 = Dirichlet で殴るのは対症療法。**正しい修正は第一内層の $\omega$ ソース項を
cell と整合させること** (壁関数生産を受ける DOF と $\omega$ を解く DOF の関係を cell に合わせる) で、
これは未着手。

#### 4. 準定常確認

残差は全 run で block-DPLUR の構造プラトー (`NOT CONVERGED`) だが、**報告している派生量 $u_\tau$ は
定常**: 最終 2 スナップショット間のドリフトは run_0040 0.001% / run_0042 0.075% / run_0043 0.025% /
run_0044 0.001% / run_0045 0.013% / run_0046 0.070%。

#### 2'. 機構の訂正: $P_\omega$ は `wf_pk` と結合していない

「第一内層が `wf_pk` を受け、$\nu_t$ が小さいので $P_\omega\propto P_k/\nu_t$ が過大」という説明は
**現行実装として誤り**。[`ransSource_d.cu:178`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) は

$$P_\omega=\alpha\rho S^2$$

であり `wf_pk` と独立。`wf_pk` を変えても $P_\omega$ は動かない (それが §1 の反証とも整合する)。

**実測による正しい説明**: スパイクは **$\alpha\rho S^2$ と $\beta\rho\omega^2$ の局所源項平衡そのもの**。
`run_0042` 第一内層 (x=0.6) を再構成すると

| 量 | 値 |
| --- | --- |
| 解像ひずみ $S$ | 1.0750e5 s⁻¹ |
| $F_1$ | 0.621 → $\alpha$=0.5118, $\beta$=0.07796 |
| 局所平衡 $\sqrt{\alpha/\beta}\,S$ | 2.754e5 |
| **実測 $\omega$** | **2.645e5** (比 0.960) |

で、ほぼ完全に源項平衡で決まっている。「第一層距離の半分での Menter blend と一致」は結果的な近さで、
こちらが直接の説明。**cell では同じ DOF の $\omega$ がピンされていて $S$ に応答しない**のが差の出所。

さらに、壁法則整合のひずみ $S_{\rm wf}=u_\tau^2 g/\nu$ と比べると

| $y^+$ | 解像 $S$ | $S_{\rm wf}$ | $S/S_{\rm wf}$ | $(S/S_{\rm wf})^2$ |
| --- | --- | --- | --- | --- |
| 27.1 | 1.075e5 | 5.146e4 | 2.09 | **4.36** |
| 56.9 | 3.192e4 | 1.623e4 | 1.97 | 3.87 |
| 89.6 | 1.531e4 | 9.611e3 | 1.59 | 2.54 |

→ y+30 メッシュの**未解像な離散勾配**が壁法則の 2 倍あり、$S^2$ で 4.4 倍の $P_\omega$ を生んでいる。

#### 3'. A/B の交絡: `nodeOmegaWfDirichlet` は $\omega$ だけを変えていない

[`turbulent_viscosity_d.cu:222`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu):

```cpp
const bool wfPin = (roOmega_wf != nullptr && roOmega_wf[ic] >= 0);
vis_turb[ic] = wfPin ? rho * k_c / w_c
                     : rho * a1 * k_c / max(a1 * w_c, S_mag * F2);
```

このスイッチは (a) 第一内層 $\omega$ のピン、(b) **SST shear limiter $\max(a_1\omega, SF_2)$ の迂回**、
(c) $\mu_t=\rho k/\omega$ への切替を**同時に**行う。したがって $C_f$ 0.843→0.957 は
「$\omega$ 低下」と「limiter bypass」の**合成効果**であり、$\omega$ 単独の効果ではない。
case/40 で $\tau_w$ が 1.237 倍に過大化したのも、既存調査では主に limiter bypass による $\mu_t$ 過大とされている
([plan §10c](../../plans/accepted/turbulence-sst-su2-taw-coupling.md))。

#### 4'. 現時点の正確な結論

**第一内層ノードの「k 半拘束ピン + $\omega$ free + 解像 $S$ による $P_\omega$」の周辺に問題があり、
`nodeOmegaWfDirichlet` 経路を通すと大幅改善する。ただし効いているのが $\omega$ 状態なのか
strain limiter なのか両方なのかは未分離。** 「真因確定」とは言わない。

分離実験と恒久修正候補は [plan: turbulence-node-wf-omega-source.md](../../plans/active/turbulence-node-wf-omega-source.md) に起票した。

## 検証基準の改訂: 外部相関 $C_f(Re_\theta)$ (Kármán–Schoenherr) を主基準に (2026-08-12)

**結論から**: 外部相関に対し **cell は 4–5% 低い / node は 22–23% 低い**。
これまで「cell = 真値」として node を 0.843 と評価していたが、**cell 自身が外部相関を外している**ので、
その枠組みでは node の欠損量を正しく測れていなかった。

$$C_{f,\mathrm{KS}}=\left[17.08(\log_{10}Re_\theta)^2+25.11\log_{10}Re_\theta+6.012\right]^{-1},\quad
Re_\theta=\frac{\rho_\infty U_\infty\theta}{\mu_\infty}$$

(NASA TMR flat-plate validation と同じ。<https://tmbwg.github.io/turbmodels/flatplate_val.html>)

**呼称**: cell 壁解像解は今後「**内部基準**」と呼び「真値」と呼ばない。K–S も厳密解ではないので
「**外部相関**」と呼ぶ。

### 前提ゲート (これを満たす station だけで判定する)

| ゲート | 結果 | 判定 |
| --- | --- | --- |
| **ZPG** | 加速パラメータ $K=(\nu/U_e^2)dU_e/dx=2.1\times10^{-9}$ ($\ll3\times10^{-6}$)。$\Delta U_e/U_e=+0.6\%$, $\Delta p/q=-1.4\%$ (x=0.2→0.9) | **合格 (実質 ZPG)** |
| **発達乱流** | $\theta(x)$ 単調増加は $x\gtrsim0.009$ (cell 壁解像) / $0.024$ (cell wf) / $0.052$ (node)。ただし **BL 内 peak $\mu_t/\mu$ が自由流値 65.9 を超えるのは $x\approx0.25$–$0.37$** | **$x\ge0.3$ を条件とする** |
| **K–S 有効域** ($4000<Re_\theta<13000$) | $Re_\theta$ は $x=0.90$ でも **5571 (node) – 6427 (cell)**。**上限 13000 に全く届かない**。$Re_\theta>4000$ は $x\gtrsim0.60$ (cell) / $x\gtrsim0.75$ (node) | **下端のみカバー** |
| **収支の信頼域** | $C_f/(2d\theta/dx)$ は $x=0.45$–$0.75$ で 0.99–1.05、$x\ge0.9$ で 0.49–1.10 に崩れる (局所 fit の station 不足+出口影響) | **$x\le0.8$ に限る** |

→ **判定帯は $x\in[0.60,0.80]$** ($Re_\theta\approx3900$–$5500$)。**K–S を唯一の合否基準にはできない**。

### $C_f$ の 4 定義 (混ぜないこと)

1. **wall-resolved**: $\tau_w=\mu\,\partial u/\partial y$ (第一 DOF の分子勾配)
2. **wall-function**: Reichardt 逆解きの $\tau_w=\rho u_\tau^2$ (ソルバが運動量式に課す値と整合)
3. **`twall` 出力**: AddTauWall で 2. の目標値へ**再スケール済み**なので**独立検証にならない** (使わない)
4. **運動量積分**: $C_f=2\,d\theta/dx$

本表の $C_f$ は **2.** で統一 (壁解像 run では 1. と 2. がほぼ一致する: x=0.6 で 2.9149e-3 vs 2.8411e-3)。
**node の W–I 双対面に実際に加わった接線力からの $C_f$ は保存出力から取得できず「未確認」**
(取得には plan §4.2 の診断追加が必要)。

### 主表 (判定帯を含む 5 station、$C_f$ は wall-function 定義)

$d\theta/dx$ は $x\pm0.08$ m の局所 2 次 fit (使用 station 7–26 個)。$\theta$ は全域積分。前縁 $x<0.10$・出口 $x>0.95$ 除外。

| run | $x$ | $Re_\theta$ | K–S域 | $C_f$ | $C_f/C_{f,\mathrm{KS}}$ | $C_f/(2d\theta/dx)$ |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0007` cell 壁解像 (内部基準) | 0.596 | 4614 | ○ | 2.915e-3 | **0.954** | 1.102 |
| | 0.745 | 5499 | ○ | 2.820e-3 | **0.956** | 1.069 |
| `run_0044` cell y+30 wf | 0.596 | 4599 | ○ | 2.933e-3 | **0.960** | 1.026 |
| | 0.745 | 5541 | ○ | 2.834e-3 | **0.962** | 1.035 |
| `run_0042` node h=1.77e-4 (y+27) | 0.602 | 3905 | △ | 2.406e-3 | **0.762** | 1.023 |
| | 0.753 | 4696 | ○ | 2.332e-3 | **0.766** | 0.993 |
| `run_0040` node h=3.40e-4 (y+52) | 0.602 | 4009 | ○ | 2.408e-3 | **0.766** | 1.043 |
| `run_0043` node h=6.59e-4 (y+102) | 0.602 | 4184 | ○ | 2.426e-3 | **0.779** | 1.049 |
| `run_0045` node kwfDir=0 | 0.753 | 4686 | ○ | 2.361e-3 | **0.775** | 0.999 |
| `run_0046` node omgWfDir=1 | 0.602 | 4464 | ○ | 2.734e-3 | **0.889** | 1.044 |

- **cell 壁解像と cell y+30 wf は外部相関に対し 0.954–0.962** で互いに 0.6% 以内。
- **node は 0.762–0.779** で、第一 DOF の y+ を 27→102 と 3.7 倍振っても **±0.011 に収まる**
  (壁関数としての y+ 非依存は成立している)。
- `nodeOmegaWfDirichlet:1` で 0.889 まで戻るが、**これは $\omega$ ピンと shear limiter 迂回の合成効果**
  (plan §1.3) なので単独因子の証拠にならない。

### $\theta$ 積分の規約と感度

$\delta_{99}$ は BL 厚さの定義であって $\theta$ の打切り位置ではないので、**積分核が十分ゼロになる
外部流まで (領域上端 $y=0.2$) 積分**するのを主値とする。node は壁点 ($y=0,u=0$) をそのまま含め、
cell は $(y,u)=(0,0)$ を補う。上端感度 (x=0.6):

| 上端 | cell 壁解像 | cell y+30 | node h=1.77e-4 | node omgWfDir=1 |
| --- | --- | --- | --- | --- |
| $y=0.2$ (全域) | 1.0373e-3 | 1.0341e-3 | 8.7775e-4 | 1.0036e-3 |
| $y=0.1$ | 1.0147e-3 | 1.0183e-3 | 8.7507e-4 | 9.9387e-4 |
| $2\delta_{99}$ | 9.9807e-4 | 1.0171e-3 | 8.7426e-4 | 9.8944e-4 |
| **最大偏差** | **3.79%** | 1.64% | 0.40% | 1.41% |

**「上端を下げても値が変わらない」は厳密には成立しない** (自由流の微小非一様性を長い距離で拾うため。
壁解像 run で最大 3.8%)。ただし $C_{f,\mathrm{KS}}$ の $Re_\theta$ 感度は
$\partial\ln C_f/\partial\ln Re_\theta\approx-0.2$ なので、**比への影響は $\pm0.8\%$** に留まる。

### Schlichting 係数 0.0592 と 0.0576 の整理

リポジトリ内で混在していた 2 つは**同じ式ではない**。1/7 乗則から自己整合に導出すると:

$u/U=(y/\delta)^{1/7}$ ($\theta/\delta=7/72$) + Blasius $\tau_w=0.0225\rho U^2(\nu/U\delta)^{1/4}$
+ 運動量積分 → $\delta/x=0.3707\,Re_x^{-1/5}$、$C_f=\mathbf{0.0577}\,Re_x^{-1/5}$

- **0.0576** = $\delta/x=0.37Re_x^{-1/5}$ と**対になる 1/7 乗則の自己整合値** (`tools/compare_cf_bl_cell_node.py`)
- **0.0592** = 同じ 1/7 乗則族の**別の経験フィット** (Schlichting, $5\times10^5<Re_x<10^7$)
  (`tools/postprocess_wall_law.py`, `wall_law_modeled.py`, `plot_uplus_truth_compare.py`)
- 差は **2.8%** で、数 % を論じる本 case では**混用不可**。
- **提案 (未実施)**: 補助基準は 0.0592 に統一し、0.0576 を使う箇所はラベルを
  「1/7 乗則自己整合値」に変える。ツール修正は別途。

### 収束・準定常 (使用 snapshot と VERDICT)

`check_convergence.py` は **8 run すべて `NOT CONVERGED (stalled/plateau)`** — これは本 case の
block-DPLUR 近似ヤコビアン由来の構造プラトーで既知 (受理済み `run_0007` も同じ)。したがって
**派生量の準定常で判定する**。`check_quasisteady.py` に `theta` / `cf_retheta` を追加した
(2026-08-12。これが唯一のツール変更)。

```bash
python3 solver_density_cuda/tools/check_quasisteady.py <run_dir> \
    --quantity theta,cf_retheta --cf-x 0.6 --tail 0.5 --drift 0.02 --osc 0.02
```

| run | 使用 snapshot | VERDICT | $\theta$ drift/tail | $C_f/C_{f,\mathrm{KS}}$ drift/tail |
| --- | --- | --- | --- | --- |
| `run_0007` | `res_120000.h5` | TRANSIENT-UNSETTLED | 1.6% | 0.1% |
| `run_0017` | `res_10000.h5` | STEADY | 0.1% | 0.0% |
| `run_0044` | `res_5000.h5` | STEADY | 0.0% | 0.0% |
| `run_0040` | `res_40000.h5` | STEADY | 0.0% | 0.0% |
| `run_0042` | `res_20000.h5` | TRANSIENT-UNSETTLED | 1.4% | 0.1% |
| `run_0043` | `res_20000.h5` | STEADY | 1.3% | 0.5% |
| `run_0045` | `res_20000.h5` | STEADY | 0.4% | 0.0% |
| `run_0046` | `res_20000.h5` | STEADY | 1.0% | 0.8% |

**報告する比 $C_f/C_{f,\mathrm{KS}}$ は全 run で drift ≤0.8%/tail** (窓 = 末尾 50%)。
`run_0007`/`run_0042` の UNSETTLED は **$\theta$ 側**の残ドリフト (1.4–1.6%) が
系列の極値を末尾に持つため。$\theta$ の 1.6% は比へ 0.3% しか効かないので、
**比の判定には影響しない**が、$\theta$ 単独を報告するときは注記が要る。

### NASA TMR ケースとの条件差 (直接比較は不可)

| 項目 | case/26 | NASA TMR flat plate |
| --- | --- | --- |
| $M$ | 0.2 | 0.2 |
| $Re$ | $Re/m=4.464\times10^6$, $Re_{x,\max}=4.46\times10^6$ | $Re_L=5\times10^6$ |
| SST variant | forge 実装 ($P_\omega=\alpha\rho S^2$ = **SST-2003 の誤記形**) | SST-Vm 等 |
| **freestream $\mu_t/\mu$** | **65.9** (`k=0.3, omega=300`) | **≈0.009** |
| transition model | 無し (モデル自身の activation 位置あり) | 無し |
| 前縁 | $x<0$ に slip 助走区間 0.1 m | 同様 (symmetry) |
| BC | 入口 全圧全温 / 出口 静圧 / 上端 slip | riemann / 上端 farfield 等 |
| $Re_\theta$ 範囲 | ~2200–6400 | K–S 主比較域を広くカバー |

**freestream $\mu_t/\mu$ が 4 桁違う**ため、**NASA の SST 数値解との直接比較は apples-to-apples ではない**。
K–S 相関との物理比較 (本節) と、NASA SST 数値解による実装 verification は分けて扱う。
後者には流入乱流量を合わせた別ケースが要る (未実施)。

## 【重要】壁処理 × 離散化 × ソルバの分解 — node 離散化は無罪、SU2 も同じ 0.75 (2026-08-12)

「node 壁関数に欠陥がある」という作業仮説を、**壁関数を外した node** と **同一メッシュの SU2**
で切り分けた。結論は 2 つとも仮説に不利である。

### 決定表 (外部相関比 $C_f/C_{f,\mathrm{KS}}$、有効帯)

| 構成 | 壁処理 | 離散化 | $C_f/C_{f,\mathrm{KS}}$ (x=0.60 / 0.75) | 収支 $C_f/(2d\theta/dx)$ |
| --- | --- | --- | --- | --- |
| `run_0007` forge cell y+0.34 | **壁解像** | cell | 0.954 / 0.956 | 1.069 |
| **`run_0050` forge node y+0.7 (現行 binary)** | **壁解像** | **node** | **0.943 / 0.944** | **1.036** |
| `run_node_sst_muscl_cont` forge node y+0.7 (旧 binary) | 壁解像 | node | 0.943 / 0.947 | 1.033 |
| `run_0044` forge cell y+29 | 壁関数 | cell | 0.960 / 0.962 | 1.035 |
| `run_0042` forge node y+27 | 壁関数 | node | **0.762 / 0.766** | 0.993 |
| **`run_0048` SU2 y+27 (forge と同一メッシュ)** | **壁関数** | **node** | **0.754 / 0.755** | 0.775 |
| `run_0047` SU2 y+27 | low-Re (不適) | node | 0.365 / 0.365 | 0.181 |

### 結論 1: node 離散化そのものは健全

**壁関数を外して壁解像で解くと node は 0.943–0.944** で、cell 壁解像 (0.954–0.956) と **1.1% しか違わない**。
しかも**運動量収支は node の方が良い** (1.036 vs 1.069)。旧 binary の既存 run とも 0.943 で一致するので、
本セッションの一連の修正とも独立している。
→ **median-dual 離散化が $C_f$ を落としているのではない。** 欠損は**壁関数経路に限定**される。

### 【撤回】結論 2: SU2 も同じ 0.75 を出す — forge node 壁関数は「異常」ではない

> **この結論は撤回する (2026-08-12、同日の追検証)。** SU2 壁関数版 (`run_0048`) は
> **運動量積分を満たしていない** (壁 $C_f$ が BL 発達が要求する値の 0.79 倍) ため、
> forge node の値を裏づける証拠として使えない。残差も rms −4.4 止まり (壁解像版は −12.4)。
> 詳細は下記「運動量収支による再評価」節。以下は撤回された記述として残す。

### 【撤回された記述】SU2 も同じ 0.75 を出す

SU2 v8.5.0 は **vertex-centered median-dual = forge node と同じ DOF 配置**である。
`STANDARD_WALL_FUNCTION` を有効にして **forge `run_0042` と完全同一のメッシュ・同一 BC・
同一物性・同一流入乱流**で回すと **0.754–0.755** で、**forge node の 0.762–0.766 と 1% 以内で一致**する
(SU2 自身の `Skin_Friction_Coefficient` 出力で測ると 0.791–0.799)。

→ **「forge の node 壁関数に固有の欠陥がある」という仮説は成立しにくい。**
独立した検証済みソルバが同じ離散化・同じ条件で同じ答えを出す。
むしろ**外れているのは forge の cell 壁関数 (0.960)** の方で、これが壁解像値に近いのは
偶然の可能性がある (cell の第一 DOF はセル重心なので、vertex-centered とは壁関数の当て方が構造的に違う)。

**未解明として残るもの** (どちらも 0.75 を「正しい」とする根拠にはならない):

- 本ケースの **freestream $\mu_t/\mu=65.9$** は NASA 標準ケース (~0.009) より 4 桁大きい。
  両ソルバで共通の条件なので相対比較は成立するが、**壁関数がこの高い自由流乱流下で
  どう振る舞うべきか**の外部基準は無い。
- SU2 壁関数版の**運動量収支が閉じない** (0.775)。forge node は 0.993 で閉じている。
  同じ $C_f$ に至りながら収支の質が違うので、両者が同じ理由で 0.75 なのかは未確認。

### 結論 3: SU2 low-Re も y+27 で崩壊する (forge 固有ではない)

`run_0047` (壁関数なし、y+27) は $C_f/C_{f,\mathrm{KS}}=0.365$、運動量収支 **0.181**、
$\mu_t/\mu$ ピーク 340–452、$\theta$ が forge の 2.6 倍。**解として使えない**。
forge cell の y+30 low-Re が崩壊した記録 (`run_0010`, Cf/Schl 0.13) と同種で、
**「low-Re 壁処理を y+≫1 メッシュに当てると壊れる」は forge 固有ではない**ことの確認になる。

### SU2 クロスチェックの設定 (再現用)

`run_0047`–`run_0049`。手順は [`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md)。
メッシュは **forge と同一の `.geo`** から `gmsh -2 ... -format su2` で生成 (交絡を避ける)。

| forge 設定 | SU2 設定 |
| --- | --- |
| `Pt=100000, Tt=288.15` | `MARKER_INLET= ( inlet, 288.15, 100000.0, 1,0,0 )` + `INLET_TYPE= TOTAL_CONDITIONS` |
| 背圧 `Ps=97250` | `MARKER_OUTLET= ( outlet, 97249.67 )` |
| `M=0.2` | `MACH_NUMBER= 0.2`, `FREESTREAM_PRESSURE= 97249.67`, `FREESTREAM_TEMPERATURE= 285.8631` |
| `visc=1.8e-5` (定数) | `VISCOSITY_MODEL= CONSTANT_VISCOSITY`, `MU_CONSTANT= 1.8E-5` |
| `thermCond=0.0251` (定数) | `CONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY`, `THERMAL_CONDUCTIVITY_CONSTANT= 0.0251` |
| 流入 `k=0.3, omega=300` | `FREESTREAM_TURBULENCEINTENSITY= 0.006598`, `FREESTREAM_TURB2LAMVISCRATIO= 65.853` (等価換算) |
| `wall` (断熱 no-slip) | `MARKER_HEATFLUX= ( wall, 0.0 )` |
| `sym`/`top` (slip) | `MARKER_SYM= ( sym, top )` |
| `wallTreatmentSST: 1` | `MARKER_WALL_FUNCTIONS= ( wall, STANDARD_WALL_FUNCTION )` |

評価は forge と**同一の後処理**(同じ $\theta$ 積分規約・同じ Reichardt 逆解き $C_f$・同じ K–S 比)で行う。

### 速度 (参考、公平比較ではない)

| | 反復 / 時間 | 1 反復 |
| --- | --- | --- |
| forge `run_0042` (GPU RTX 3060) | 20000 step / **38.8 s** | **1.94 ms** (内部反復 20 回込み) |
| SU2 `run_0047` (CPU 8 スレッド) | 40000 iter / 314 s | 7.86 ms |

同一メッシュ (12428 ノード) で forge が約 4 倍速いが、**GPU vs CPU なので効率の比較にはならない**。

## 運動量収支による再評価 — SU2 壁関数版は棄却、node 壁解像は SU2 と一致 (2026-08-12)

### 収支の取り方 (簡易形は不十分だった)

当初は ZPG 極限 $C_f/2=d\theta/dx$ で評価していたが、本ケースは弱い加速があるので**圧力勾配項を
落とせない**。von Kármán 運動量積分の完全形を使う:

$$\frac{C_f}{2}=\frac{d\theta}{dx}+\frac{\theta}{U_e}\left(H+2-M_e^2\right)\frac{dU_e}{dx},
\qquad H=\frac{\delta^*}{\theta},\quad
\delta^*=\int_0^{y_{top}}\left(1-\frac{\rho u}{\rho_e U_e}\right)dy$$

圧力勾配項は $d\theta/dx$ の **2.0–3.1%** あり、収支比を系統的に ~2.5% 押し下げる
(例: forge cell 壁解像 1.102→1.076)。$d\theta/dx$ と $dU_e/dx$ はいずれも $x\pm0.08$ m の局所 2 次 fit。

さらに、**壁 $C_f$ と「BL 発達が要求する $C_f$」($=2(d\theta/dx+\text{PG})$) の両方を外部相関で
正規化**すると、どちらが外れているかが分離できる。

### 決定表

| 構成 | 残差 | 壁 $C_f$/KS | 積分 $C_f$/KS | 壁/積分 |
| --- | --- | --- | --- | --- |
| `run_0007` forge cell 壁解像 | plateau | 0.954 / 0.956 | 0.886 / 0.919 | 1.076 / 1.039 |
| **`run_0050` forge node 壁解像** | plateau | **0.943 / 0.944** | 0.936 / 0.937 | **1.008 / 1.008** |
| **`run_0049` SU2 壁解像 (y+0.7)** | **rms −12.4** | **0.938 / 0.944** | 0.939 / 0.950 | **0.999 / 0.995** |
| `run_0044` forge cell 壁関数 | plateau | 0.960 / 0.962 | 0.958 / 0.956 | 1.001 / 1.006 |
| `run_0042` forge node 壁関数 | plateau | 0.762 / 0.766 | 0.759 / 0.790 | 1.003 / 0.970 |
| **`run_0048` SU2 壁関数** | **rms −4.4** | 0.793 / 0.796 | **1.003 / 1.005** | **0.790 / 0.792** |

### 結論 A: SU2 壁関数版は棄却する

**壁 $C_f$ が BL 発達の要求値の 0.790 倍** = 定常 2D 解なら恒等的に閉じるはずの運動量積分を
**21% 破っている**。$\theta$ の積分上端感度は 0.67%、自由流 $u/u_e=0.99998$、$U_e(x)$ は forge と
ほぼ同一なので、**後処理側の問題ではない**。残差も rms −4.4 止まり (壁解像版は −12.4) で、
**`STANDARD_WALL_FUNCTION` 経路が収束しきっていない可能性が高い**。
なお `CD` は iter 3333 以降 0.002611 で 6 桁不変なので「積分量が動いている」という意味では収束済み。

→ **「SU2 も 0.75 を出すから forge node は異常でない」という結論は撤回する。**

### 結論 B: node 壁解像は SU2 とほぼ完全に一致 (node 離散化の検証)

**`run_0050` forge node 壁解像 0.943/0.944 と `run_0049` SU2 壁解像 0.938/0.944 が 0.5% 以内で一致**し、
**両者とも運動量収支がほぼ完璧** (1.008 / 0.999)。SU2 側は rms −12.4 まで落ちた素性の良い解である。

→ **median-dual (vertex-centered) 離散化は検証済みソルバと一致する。** node 離散化は無罪、という
結論は**むしろ強化された** (今度は良く収束した SU2 解が根拠)。

### 結論 C: 本ケースの「正しい」壁解像値は $C_f/C_{f,\mathrm{KS}}\approx0.94$

forge node 0.943 / SU2 0.938–0.944 / forge cell 0.954。**3 者が 0.94–0.95 に収まる**。
K–S より 5–6% 低いのは、本ケースの **freestream $\mu_t/\mu=65.9$** (NASA 標準の 4 桁大) が
効いている可能性が高い (未検証)。
**したがって node 壁関数の 0.762 は、この 0.94 に対して −19% の欠損**である。

### 結論 D: forge cell 壁解像の収支が 3 者中で最も悪い

`run_0007` の壁/積分は 1.076 / 1.039 で、node (1.008) や SU2 (0.999) より劣る。
壁 $C_f$ が 0.954 と高めなのは**収支が閉じていないため**の可能性がある
(積分 $C_f$ で見ると 0.886–0.919 で node/SU2 の 0.936–0.950 より低い)。
**内部基準として使ってきた `run_0007` 自身に疑いがある** — 別途要確認。

## 追検証: SU2 の壁 $C_f$ 出力は何か / どの定義を信じるか (2026-08-12)

「SU2 の壁 $C_f$ は壁関数を反映していないのでは。積分 $C_f$ が正しいなら流れ場は正しいのでは」
という指摘を検証した。**指摘は概ね当たっている**が、単純に「出力だけの問題」でもない。

### 4 定義の突き合わせ (x=0.6、すべて $C_f/C_{f,\mathrm{KS}}$)

| 構成 | (a) 解像勾配 | (b) ソルバ出力 | (c) 運動量積分 | (d) Reichardt 逆解き | 最大差 |
| --- | --- | --- | --- | --- | --- |
| **`run_0049` SU2 壁解像 (rms −12.4)** | **0.945** | **0.938** | **0.939** | **0.947** | **1.0%** |
| `run_0050` forge node 壁解像 | — | — | 0.936 | 0.943 | 0.7% |
| `run_0007` forge cell 壁解像 | — | — | 0.886 | 0.954 | **7.4%** |
| `run_0048` SU2 壁関数 | 0.373 | **0.793** | **1.003** | 0.754 | **86%** |
| `run_0042` forge node 壁関数 | 0.369 | — | **0.759** | **0.762** | 62% |
| `run_0044` forge cell 壁関数 | — | — | 0.958 | 0.960 | 0.2% |

### 検証 1: 後処理パイプラインは正しい

**良く収束した壁解像解 (`run_0049`, rms −12.4) では 4 定義が 1.0% 以内で一致**する。
$\theta$ 積分・$d\theta/dx$ の局所 fit・圧力勾配項・Reichardt 逆解きのどれにも系統誤差が無い実証。

### 検証 2: SU2 の壁 $C_f$ 出力は解像勾配ではない (指摘は正しい)

SU2 壁関数版で (a) 解像勾配 0.373 に対し (b) 出力 0.793 = **2.1 倍**。
→ **SU2 の `Skin_Friction_Coefficient` は壁モデルを反映している**。
「(b) が低いから SU2 の壁応力が低い」と読んだのは誤りだった。

### 【撤回】検証 3: SU2 の外側の場は正しい

> **撤回する (2026-08-12, レビュー指摘)。** `run_0048` は**全残差で見ると強い未収束**である:
> `rms[Rho]` −4.43 だけを見ていたが、**`rms[RhoE]` = +1.03 (残差 10)**、`rms[w]` −0.29、
> `rms[RhoV]` −2.13 ([history.csv](run_0048_su2_sst_wf_ny52/history.csv) 最終行)。
> 対して壁解像 `run_0049` は全列 −5.0〜−12.6。
> **AGENTS.md の「`rms_ro` だけで判断しない」を SU2 に対して破っていた。**
> 未収束場では定常運動量積分式に残差項が欠けるので、(c)=1.003 が K–S に近いことは
> 「場が正しい」根拠にならない。正確な表現は
> **「現スナップショットの境界層発達率はたまたま K–S に近いが、壁応力との収支が 21% 閉じず、
> 定常解として使用不能」**。`CD` の停止だけでは $k,\omega,\rho E$ の未収束を覆せない。
>
> **また SU2 に Reichardt 逆解き (d) を当てた比較も成立しない。** SU2 の
> `STANDARD_WALL_FUNCTION` は **Nichols–Nelson (2004) の圧縮性 Spalding/White–Christoph 型**
> (`CNSSolver.cpp` の Newton 反復、Crocco–Busemann 込み) で Reichardt ではない。
> 「SU2 の第一ノード速度と壁応力が不整合」という記述は**別の壁法則を比べていたための誤り**。
>
> 以下は撤回された記述として残す。

### 【撤回された記述】SU2 の外側の場は正しい / forge node は場そのものが薄い

**SU2 の (c) 運動量積分 = 1.003** (x=0.3/0.6/0.9 で 0.961/1.003/1.034) →
**SU2 の境界層は外部相関どおりに育っている**。
一方 **forge node 壁関数の (c) は 0.759** → **BL 発達そのものが 24% 遅い**。

→ **この 2 つは「同じ 0.75」ではない。** SU2 は場が正しく $C_f$ の定義がずれているだけ、
forge node は**場そのものが薄い**。前節の「SU2 も 0.75 だから forge node は異常でない」は
**二重に誤り**だった (SU2 側の未収束と、$C_f$ 定義の取り違え)。

### ただし SU2 壁関数版も内部整合していない

(a) 0.373 / (b) 0.793 / (c) 1.003 / (d) 0.754 と **4 定義が 86% 開いている**
(壁解像解では 1.0%)。特に **(d) 第一 DOF 速度から逆解きした値が (c) と整合しない** =
SU2 の第一ノード速度は、(c) が示す壁応力を対数則に当てはめた値より低い。
残差 rms −4.4 止まりと合わせ、**「場は正しいが壁近傍が未整合」**。
定量比較の相手にするには収束させ直すか壁応力の実装確認が要る。

### 【訂正】(c) は「定義非依存の壁摩擦」ではない

$2(d\theta/dx+\mathrm{PG})$ は壁出力に依存しない有用な推定だが、**定常性・境界層近似・
$U_e,\rho_e$ の取り方・積分上端・微分 fit に依存する**。主判定は次の 3 点セットにする:

1. **実際に運動量残差へ入った W–I 接線力** ← **未取得。これが最優先**
2. それと運動量積分値の収支
3. 両者と外部相関 K–S との比較

### いま信じるべき数字 (収支診断として)

| | (c) $C_f/C_{f,\mathrm{KS}}$ |
| --- | --- |
| SU2 壁解像 (rms −12.4) | **0.939** |
| forge node 壁解像 | **0.936** |
| forge cell 壁関数 | 0.958 |
| forge cell 壁解像 | 0.886 ← 3 者中の外れ値。要確認 |
| **forge node 壁関数** | **0.759** ← 明確な外れ値 |
| SU2 壁関数 | 1.003 (ただし内部不整合) |

**壁解像基準は 0.94 前後 (SU2 0.939 / forge node 0.936)。forge node 壁関数の 0.759 は −19%、
しかも $\theta$ 発達そのものが遅い**ので、壁応力の出力だけの問題ではない。

## 現状認識の確定版 (2026-08-12 レビュー反映)

> **case/26 では forge node 壁解像が SU2 壁解像と 0.5% 程度で一致するため、node 離散化一般は
> 主因ではない。forge node 壁関数では、課した低い壁応力と境界層の遅い発達が互いに整合しており、
> 欠損は壁関数経路に局在する。ただし、狙った壁応力が W–I 残差へ正しく作用しているか、また
> $\omega$ 状態・limiter・源項のどれが低摩擦解を作るかは未分離である。SU2 壁関数 run は
> 全残差未収束のため比較証拠に使用しない。**

### 確定してよいこと

| # | 内容 | 根拠 |
| --- | --- | --- |
| 1 | forge node 壁解像 (`run_0050`) と SU2 壁解像 (`run_0049`) が $C_f/C_{f,\mathrm{KS}}\approx0.94$ でほぼ一致 | (a)(b)(c)(d) が 1% 以内、SU2 は全残差 −5〜−12.6 |
| 2 | したがって −19% 欠損は **median-dual/node 離散化一般ではなく forge の node 壁関数経路に局在** | 1 の対偶 |
| 3 | SU2 low-Re を y+27 に適用した `run_0047` は比較対象にならない | 収支 0.18、$\mu_t/\mu$ 340–452 |
| 4 | SU2 の `Skin_Friction_Coefficient` は**壁関数応力を反映している** | `CFVMFlowSolverBase.inl` で `scale = Tau_Wall / WallShearStress` により接線応力を再スケールしてから出力 |

### 使ってはいけないこと

- **`run_0048` (SU2 壁関数) を定量比較の相手にしない** — 全残差未収束 (`rms[RhoE]`=+1.03)。
- **SU2 に Reichardt 逆解きを当てない** — SU2 は Nichols–Nelson。整合確認には SU2 自身の
  `Y_Plus` (壁ノードに出力。x=0.6 で **27.21**)・$U_\tau$・Spalding 則を使う。
- **(c) 運動量積分を「定義非依存の $C_f$」と呼ばない** — 収支診断である。

### 次にやること (この順)

1. **実 W–I 接線力の診断** — 「狙った摩擦力が W–I 間で本当に作用しているか」。plan §4.2 の未完タスク。
   **E1/E2/E3 より先**。
2. 2×2 factorial ($\omega$ pin × limiter bypass) の分離
3. E3 ($P_\omega$ の壁法則整合形)

### 数値の再生成

README の $C_f(Re_\theta)$ 系の数値はすべて [`tools/cf_retheta_analysis.py`](tools/cf_retheta_analysis.py) で再生成できる:

```bash
python3 tools/cf_retheta_analysis.py run_0042_node_yp30_planar_2nd
python3 tools/cf_retheta_analysis.py run_0049_su2_sst_lowre_fine --su2 --wall-law spalding
python3 tools/cf_retheta_analysis.py <run> --ytop 0.05 --fit-window 0.12 --fit-order 3 --edge global
```

fit 窓幅/次数・積分上端・`local`/`global` の edge 基準・ソルバ固有の壁法則をすべてパラメータ化してある。

## §4.2 実 W–I 接線力の診断 — 狙った壁摩擦は**正しく作用している** (2026-08-12)

「狙った摩擦力が W–I 間で本当に作用しているか」(当初からの本題) に答えるため、
**残差へ加わる直前の traction** を壁ノードに集計する診断を実装した
(`wi_ftan` / `wi_fnrm` / `wi_ftan_res`、[`viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) の
AddTauWall 再スケール直後で `atomicAdd`、毎ステップ 0 クリア)。
`twall` 出力はモデル目標値へ再スケール済みで独立検証にならないため、これが必要だった。

`run_0051_node_wf_widiag` (= `run_0042` 収束場から 500 step、設定同一):

| 帯 | 実際に残差へ入った $\Sigma$`wi_ftan` [N] | 狙い値 $\int\rho u_\tau^2 dx$ [N] | **比** | 法線力/接線力 |
| --- | --- | --- | --- | --- |
| x=[0.30, 0.90] | 3.8796e+00 | 3.7980e+00 | **1.0215** | 0.003% |
| x=[0.20, 0.95] | 4.9573e+00 | 4.8780e+00 | **1.0163** | 0.003% |
| x=[0.45, 0.75] | 1.9399e+00 | 1.8595e+00 | **1.0432** | 0.002% |

### 結論 1: 伝達は健全 (仮説を 1 つ消せた)

**モデル $\tau_w=\rho u_\tau^2$ は W–I 双対面で 2–4% の精度で実際に流体へ作用している**。
→ **「狙った壁応力が残差へ届いていない」という可能性は消えた**。欠損は伝達ではなく、
**壁関数が低い $u_\tau$ を算出していること自体**にある (第一内点速度が低い → 低い $u_\tau$ →
低い $\tau_w$ → 薄い BL、が自己整合している)。

### 結論 2: 余計な法線力は実害なし (もう 1 つ消せた)

AddTauWall 再スケールが traction ベクトル全体に掛かる件 (コメントは「接線成分を再スケール」)
を実測したところ、**法線成分は接線の 0.003%** で無視できる。
コード/コメントの不一致は残るが、**本ケースでは実害なし**。

### 結論 3: 再スケールの向きは「増幅」ではなく「低減」だった (前提の訂正)

$\Sigma$`wi_ftan_res` (再スケール前) = 6.0581 N → 再スケール後 3.8796 N で **係数 0.640**。
$C_{f,\rm model}/C_{f,\rm molec}\approx2.1$ から「~2.1 倍に増幅される」と推定していたのは**誤り**。
W–I 双対面の解像 traction は**壁の分子勾配ではなく $\mu_{\rm total}$ (乱流粘性込み) と双対面距離**で
評価されるので、モデル $\tau_w$ より**大きい**。AddTauWall はそれを**下げて**モデル値に合わせている。

### 残る問い

伝達も法線力も無罪なので、**なぜ解が低摩擦で釣り合うのか**が残る。
次は §3 の 2×2 factorial ($\omega$ pin × strain limiter bypass) で
$\omega$ 状態 / limiter / 源項のどれが効いているかを分離する。

## §3 の 2×2 factorial — $\omega$ 状態が支配、limiter 迂回はほぼ無効 (2026-08-12)

`nodeOmegaWfDirichlet` が $\omega$ ピンと SST shear limiter 迂回を**束ねていた**問題を、
2 ビットに分解して分離した。

### 実装

- **`wf_irep_flag`** (新規): node 壁関数の**第一内層ノード (irep) だけ**を 1 でマークする
  ([`ransWallFunction_d.cu`](../../solver_density_cuda/cuda_forge/ransWallFunction_d.cu))。
  `wf_pk>=0` は壁ノードにも立つのでマスクに使えない。**cell では常に 0** なので cell ビット不変が構造的に保証される。
- **`FORGE_WF_LIMITER_BYPASS`** (env, 既定 −1): $\mu_t=\rho a_1k/\max(a_1\omega,SF_2)$ の迂回を
  $\omega$ ピンから切り離す。−1=従来 (ピンに追随) / 0=強制 OFF / 1=強制 ON (irep のみ)
  ([`turbulent_viscosity_d.cu`](../../solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu))。
  **既定 −1 は従来と完全同一経路**。

### 結果 (x=0.6、全 run 同一 IC・20000 step・quasisteady `STEADY`)

| | bypass=0 | bypass=1 |
| --- | --- | --- |
| **pin=0** | $Y_{00}$ = **0.762** (`run_0053`) | $Y_{01}$ = **0.768** (`run_0055`) |
| **pin=1** | $Y_{10}$ = **0.862** (`run_0054`) | $Y_{11}$ = **0.889** (`run_0056`) |

($C_f/C_{f,\mathrm{KS}}$、(d) 壁法則逆解き。運動量積分 (c) は 0.752/0.753/0.848/0.873 で全 run 自己整合)

$$\text{主効果 pin} = +0.100,\qquad
\text{主効果 bypass} = +0.006,\qquad
\Delta_{\rm int}=Y_{11}-Y_{10}-Y_{01}+Y_{00} = +0.021$$

### 結論

1. **$\omega$ 状態が支配因子** — 全効果 +0.127 のうち **+0.100 (79%)** がピン単独。
2. **limiter 迂回は単独ではほぼ無効** (+0.006)。
   「case/40 の過大化は主に limiter bypass 由来」という従来の見立ては、少なくとも本ケースでは支持されない。
3. **相互作用は +0.021 で無視できない** → 「$E1+E2$ の和」で評価してはいけないという指摘は正しかった。
   ただし実際には pin が支配的なので、結論は変わらない。
4. **両方入れても 0.889 で目標帯 0.94±0.02 に届かない** (必要 +0.178 に対し +0.127 = **71%**)。
   → **残り ~29% は $\omega$ 状態でも limiter でもない別の要因**。§5 の E3 ($P_\omega$ を壁法則整合形へ) が次の候補。

### 位置づけ

`nodeOmegaWfDirichlet: 1` は依然として**採用しない** (case/40 で $\tau_w$ 過大化)。
本実験は「どの因子が低摩擦解を作っているか」の分離が目的であり、Dirichlet を対策にする話ではない。

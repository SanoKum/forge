# KEEP 用 entropy-stable 散逸レイヤ (段階導入)

## メタ

- **area**: `convection`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/convection/implementation.md`](../../methods/convection/implementation.md) (`KEEP_d` 節)
  - [`methods/convection/theory.md`](../../methods/convection/theory.md)
- **related_plans**:
  - [`convection-central-scheme-oscillation-control.md`](../../notes/investigations/convection-central-scheme-oscillation-control.md) — 本計画の背景サーベイ (振動源4分類・§10 検証ラダー)
  - [`turbulence-iddes-sst.md`](turbulence-iddes-sst.md) — Phase 1.5 低散逸 flux の受け皿 (KEEP を LES 域で使う動機)
  - [`time_integration-lowmach-preconditioning.md`](../accepted/time_integration-lowmach-preconditioning.md) — `lowMachCprime` (c') の出典
  - [`time_integration-general-eos-jacobian.md`](../accepted/time_integration-general-eos-jacobian.md) — Step 2 matrix 版の固有系素材
- **created**: `2026-07-07`
- **owner**: `CFD Dev`

## 1. 目的

純粋 KEEP (`solver: KEEP`) は散逸ゼロの中心流束であり、(a) 低マッハで圧力 odd-even 市松 (checkerboard) を減衰できない (離散 null-mode)、(b) EC/KEP 系固有の線形不安定 (anti-diffusion, Ranocha-Gassner) を持つ。**`lowMachPrecond` は SLAU 散逸の中の c' 置換であり、散逸を持たない KEEP には作用しない** (前処理する対象が無い)。

本計画は KEEP の中心流束を汚さずに、**独立な opt-in 散逸レイヤ**を段階導入する:

$$F_f = F_f^\mathrm{KEEP} - \tfrac12\,\sigma\,D_f$$

完了時: LES/ILES 用の KEEP が低マッハでも市松らず、KE 保存性の劣化が定量管理された状態で回る。

## 2. スコープ

**やる (段階導入)**:

| Step | 内容 | 状態 |
| --- | --- | --- |
| 1 | **単成分 CPG・スカラー ES 散逸**: $D_f=\lambda'\,\Delta U$, $\lambda'=\|U_n\|+c'$ (`lowMachPrecond>=1` で `lowMachCprime`、else $c$)。config `keepDissType`(0=off 既定・ビット不変/1=scalar)・`keepDissCoeff`(σ) | 本実装 |
| 2 | 単成分 CPG・**matrix 版**: $D_f=R\|\Lambda'\|SR^{\mathsf T}\Delta w$ (Chandrashekar 2013, entropy-scaled)。音響のみ $\|U_n\|+c'$、せん断/エントロピーは $\|U_n\|$ | **実装+検証済** |
| 3 | TP 単成分 ($s^0(T)$ 多項式=`thermo_s0_mass` 流用) | 将来 |
| 4 | **多成分** (Chalot-Hughes-Shakib のエントロピー変数 + Gouasmi スケーリング + `thermoHrefTemp` datum/ゼロ濃度種対策) | 将来 |

**やらない**:
- 衝撃捕獲 (衝撃レイヤは当面 SLAU ブレンド、中期に directional LAD を別 plan で)
- KEEP 中心流束自体の変更 (KEEP-PEP/AEC 化は別 plan)
- RANS↔LES f_d ブレンド ([turbulence-iddes-sst](turbulence-iddes-sst.md) §4.8 側)

## 3. 理論 (Step 1)

スカラー LLF/Rusanov 型散逸 $D=\lambda'\Delta U$ ($\lambda'>0$) は、エントロピー $\eta$ の凸性より
$\Delta w^{\mathsf T}\Delta U = \Delta w^{\mathsf T}\bar H\Delta w \ge 0$ ($\bar H$=経路平均の $\partial U/\partial w$、SPD) で**エントロピー散逸的** (Tadmor)。低マッハスケール $c\to c'$ は $\lambda'>0$ を保つので ES 性は不変。

- **KE 非増加は保証しない** (Chandrashekar の釘: KE 安定には別条件)。σ を小さく取り、L2 (TGV) で KE 劣化を実測管理する。
- **市松への効果**: 面ジャンプ $\Delta U$ は 2Δ モードを最大感度で見る (中心流束の null-mode と対照的)。L1 の市松圧力摂動は $\Delta(\rho e)=\Delta p/(\gamma-1)$ 経由で減衰する。
- σ の既定 0.15 は仮置き。L1 (減衰) と L2 (KE 保存) のトレードオフで較正する。

## 4. 実装設計 (Step 1)

- **config** (`solverConfig.hpp/.cpp`): `keepDissType` (int, 0 既定, 0..1)・`keepDissCoeff` (flow_float, 0.15 既定)。`keepDissType==0` で既存 `KEEP_d` とビット不変。
- **カーネル** (`convectiveFlux_keep_d.inc.cuh`): 引数追加 `(keepDissType, keepDissCoeff, lowMachPrecond, precondEps)`。face で保存量ジャンプ $\Delta U$ (primitives から構成)・$\lambda'=|U_n|+c'$ を計算し、`temp -= 0.5*σ*λ'*ΔU*sss` を全 5 式に加算。`massflux[ip]` は散逸込みの総質量流束 (スカラー輸送と整合)。
- **cell/node 両対応**: CV 抽象のまま (面ジャンプ+面法線のみ使用)。node は内部双対面のみ (境界は既存どおり別カーネル)。
- **陰解法**: RHS のみの変更。LHS は defect-correction で凍結のまま (収束解は σ に依存して変わる=スキームが変わるため当然)。

## 5. 検証 (サーベイ §10 ラダー)

- **L1 (case/35, 激軽)**: 一様低マッハ流 ($M_0=0.1$, $\rho_0=1$, $P_0=1/\gamma$) + 市松圧力摂動 $P=P_0[1+\epsilon(-1)^{i+j+k}]$, $\epsilon=10^{-3}$ (python で input h5 の VALUE を直接生成、パリティは重心から)。
  - **L1a**: pure KEEP (σ=0) で初期 RHS が機械ゼロ (null-mode の実証)。
  - **L1b**: σ>0 で市松振幅 $E^{HF}_p$ が減衰。
- **L2 (case/09, 軽)**: TGV M0.4 32³ explicit RK4。σ=0 が既存とビット不変・σ>0 の KE 減衰率を定量化 (SLAU 比で桁小)。
- 全 run: 新規 `run_*`・README run 表・NaN 監視・`residual_history.png`。

## 6. リスク・注意

- **stale build**: config struct 変更 → **full rebuild 必須** (既知の罠)。
- σ 過大は L2 で KE を殺し LES の意味を失う。σ 過小は L1b で市松が残る。両ゲートで挟む。
- Step 4 (多成分) の float32 桁落ちは `thermoHrefTemp` datum + thermo 層が double であることで対処見込み (サーベイ open question 7 の評価済み)。

## 変更ログ

- `2026-07-07` — 計画作成 (Step 1 実装開始)。背景は [convection-central-scheme-oscillation-control.md](../../notes/investigations/convection-central-scheme-oscillation-control.md) §5.5/§6/§10。
- `2026-07-19` — **Step 1 実装+検証完了**。config `keepDissType`/`keepDissCoeff` + `KEEP_d` に scalar ES 散逸 (λ'=|Un|+c', lowMachPrecond>=1 で `lowMachCprime`)。full rebuild 後 L1/L2 検証:
  - **L1a (null-mode 実証, `case/35...run_0014_cell_keep_cbd_pure`)**: 市松圧摂動に対し pure KEEP の全 rms 残差が**厳密ゼロ** = 市松は純 KEEP の厳密な離散 null-mode (机上予言と一致)。SLAU 対照 (`run_0016`) は即座に非ゼロ。
  - **L1b (減衰, `run_0015` σ=0.15 / `run_0017` σ=0.05)**: A_cb 1e-3 → 8.0e-10 (6桁) / 1.3e-7 (~4桁)。散逸は rms_roe のみ非ゼロ (Δp 経由で設計通り選択的)。
  - **L2 (KE cost, `case/09...run_0022/0023/0024`)**: σ=0 は既存参照と同挙動 (−0.33%, 不変確認)。σ=0.15: KE −8.37%、σ=0.05: −2.71%。
  - **σ 既定を 0.05 に決定** (市松~4桁減衰/400step と KE cost 2.7% のバランス)。ケース次第で 0.15 まで上げてよい。
  - 全 run NaN なし。これらは非定常保存テストであり定常収束は主張しない (VERDICT は record として各 run に生成済)。
  - 残: Step 2 (matrix 版)・Step 3 (TP)・Step 4 (多成分) は未着手。node モードでの L1 相当も未実施 (cell のみ検証)。
- `2026-07-19` — **Step 2 (matrix ES) 実装+検証完了** (`keepDissType: 2`)。
  - **式検証を実装前に実施** (`solver_density_cuda/tools/verify_entropy_scaling.py`): $w=\partial\eta/\partial U$・$S=R^{-1}HR^{-\mathsf T}$ が対角で文献値 (音響 $\rho/2\gamma$・エントロピー $\rho(\gamma-1)/\gamma$・せん断 $p$) に一致・小ジャンプで Roe 型と漸近一致・$Q$ SPD、を数値確認 (PASS)。
  - 実装: $D=R|\Lambda'|SR^{\mathsf T}\Delta w$、$|\Lambda'|$ は**音響のみ** $|U_n|+c'$ (lowMachPrecond>=1 で `lowMachCprime`)、せん断/エントロピーは $|U_n|$。$\Delta w$ は logf×2/face のみの追加コスト。
  - **L1b (`case/35...run_0018_cell_keep_cbd_dissmat005`)**: A_cb 1e-3 → **7.1e-8** (scalar σ=0.05 の 1.3e-7 より良)。
  - **L2 (`case/09...run_0025_cell_keep_dissmat005_l2`)**: KE cost **1.36%** = scalar (2.71%) の**半分**。市松減衰同等以上で渦を守る選択性を実証。
  - **推奨**: LES 用途は `keepDissType: 2` を第一候補 (σ=0.05)。scalar (1) は頑健フォールバック。
  - 残: Step 3 (TP 単成分)・Step 4 (多成分)・node モード検証・L3 (粘性 TGV Re=1600 で WALE との共存確認)。
- `2026-07-19` — **c' の要否切り分け + `keepDissCprime` ノブ分離** (ユーザー指摘「lowMachPrecond=1 への相乗りが不安」を受け):
  - **フル c 実測** (`case/35...run_0019` / `case/09...run_0026`, matrix σ=0.05): 市松掃除は 9.2e-9 と最速だが **TGV KE cost 4.35% = c' 比 3.2 倍悪化** (scalar+c' の 2.71% より悪い)。渦の ΔUn が音響固有ベクトルに射影されフル c で食われる = Guillard-Viozat/Thornber の c·Δu 汚染の実測。**c' は市松退治には不要だが KE 保護に実質的に効く**。
  - **設計修正**: KEEP 散逸レイヤの音響波速を `keepDissCprime` (既定 1=c') で制御し **`lowMachPrecond` から完全独立化** (グローバル旗への相乗り廃止)。c' は散逸係数のみに入り伝播・中心flux・LHS は不変。SLAU の c' への不信はこの用途には転移しない (opt-in σ=0.05 の追加項の係数、λ_A はゼロにならない、M≥1 で厳密 c 復帰)。
  - ノブ動作確認 (`run_0020`): lowMachPrecond=0 + 既定で c' 挙動を再現 (7.5e-8 ≈ 7.1e-8, atomicAdd ノイズ級)。
- `2026-07-19` — **σ 掃引で「c' は σ 弱化で代替できるか」を検証 (ユーザー問い) → 単一マッハ領域では YES と判明・前回の「3.2倍悪い」を訂正**:
  - フル c σ 掃引 (`case/35...run_0021/0022`, `case/09...run_0027/0028`): **σ=0.015 で市松 3.9e-8・KE 1.10% = c' σ=0.05 (7.1e-8, 1.36%) と同等以上の Pareto 点に到達**。前回比較は同一 σ 固定で設計比較としては誤導だった。
  - 市松減衰は σ に非線形 (σ=0.04: 3.1e-10 が σ=0.05: 9.2e-9 より速い**非単調** = 2Δ 定在音響波の過減衰領域)。
  - **c' の残存価値はマッハ混在流のみ**: 「σ を M に応じ下げる」の面ごと自動版が c'。単一領域ならフル c+小 σ で十分、チャンバー M0.05+プルーム M3 同居 (ノズル/ピントル) ではグローバル σ が両立できず c' が優位。
  - **運用指針**: 単一低マッハ領域 = `keepDissCprime: 0` + σ~0.015 / マッハ混在 = `keepDissCprime: 1` + σ=0.05 (既定)。

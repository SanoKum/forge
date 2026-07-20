# case/39.periodic_hills — 周期丘 (ERCOFTAC periodic hills, Re=10595)

剥離・再付着を含む標準検証ケース。SST-DDES ([turbulence-iddes-sst plan](../../plans/active/turbulence-iddes-sst.md))
の grey-area 挙動確認と、将来の SLA (Shear-Layer-Adapted Δ) 導入の判別ケース。
駆動は文献同様「x/z 周期 + 質量流量一定の動的体積力」
([timeint-bodyforce-massflow-control plan](../../plans/active/timeint-bodyforce-massflow-control.md))。

## 幾何・条件

- 丘形状: NASA TMR 配布 `hill-geometry.dat` の 6 区間多項式 (mm 単位, h=28mm)。正本は
  https://tmbwg.github.io/turbmodels/Other_LES_Data/2Dhill_periodic/hill-geometry.dat
- 領域: $9h \times 3.035h \times 4.5h$ (x/z 周期)、山頂が x=0/9h。上壁平坦
- Re=10595 ($U_b h/\nu$, $U_b$=丘頂断面 y∈[h,3.035h] のバルク速度)
- forge 設定: CPG 空気, T_w=300K 等温壁 (体積力仕事の排熱), M≈0.15 (U_b=52 m/s),
  μ=1.617e-4 (定数) → ν=1.374e-4
- 目標体積平均運動量密度 $\langle\rho u_x\rangle_V$ = 44.12 kg/m²/s (=ρ₀U_b·A_crest·L_x/V)
- 初期体積力推定 400 Pa/m (動的制御が上書き)

## 参照データ (裏取り済み・2026-07-20 調査)

- **Fröhlich et al. 2005 JFM 526** — 参照 LES (~5M セル)。x/h=0.05,0.5,2,4,6,8 のプロファイル、再付着 x_r/h≈4.6-4.7
- **Breuer et al. 2009 Comput. Fluids 38** — Re 掃引 DNS/LES 2 コード (13M 格子)
- **ESAIM Proc. 16 (2007) 133-145** — 5 コード DES 比較。**標準 DES 格子 160×100×60 (~1M), Δx⁺,Δz⁺<35, Δy₁⁺≲1** ← 本ケースの格子仕様の根拠
- **ATAAC D3.2-36** (ERCOFTAC KBwiki) — DDES/IDDES @160×160×60。素の DDES/IDDES で参照とよく一致
- データ入手: ERCOFTAC KBwiki UFR 3-30 (https://www.kbwiki.ercoftac.org/w/index.php/UFR_3-30_Description)

## メッシュ

`mesh/make_hill_mesh.py` → `hill_des_160x100x60.msh` (gmsh4.1 直接書き出し, 構造 hex)。

- 160×100×60 = 96 万セル (node 変換で 99.2 万 CV)。x/z 一様 (Δx=0.056h, Δz=0.075h)
- y: 両側幾何ストレッチ (比 1.0961)。y₁=1.5e-3h (平坦部) / 1.0e-3h (丘頂) → y₁⁺≲1 (低 Re 壁解像)
- 品質: **PASS** (AR max=74.6, skew max=0.451; cell 変換で判定)

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_rans_init` | 段階起動スモーク: 定常 SST-RANS (SLAU 1次, lowMachPrecond2, blockDPLUR cfl_pseudo2, 固定 bodyForce 400)。**静止 IC** (u=0, P/T 一様 = 音響過渡ゼロ) から体積力スピンアップ 2000 step | 完走・NaN 無し。スピンアップは健全だが遅い (Ux max 3.4 m/s どまり) → 成形 IC の run_0002 へ | ref (起動方式確認) |
| `run_0002_rans_dev` | RANS 発達本番: run_0001 res_2000 の ω 場を保持しつつ**成形速度 IC** (列ごと質量保存パラボラ u max 78) に張替え、20k step。⚠ cfl_pseudo=10 は step 35 で NaN。**cfl_pseudo=2 で安定** | 完走 exit0・NaN なし・全残差 2-2.6 桁低下 (VERDICT: NOT CONVERGED=継続収束中、DDES の IC としては十分)。固定 β=400 の平衡は M/V=27.4 (目標 44.1 には β~1000 要) と判明。res_20000.h5 | ref (DDES の IC 供給元) |
| `run_0003_ddes_smoke` | **素の DDES + bodyForceCtrl 統合スモーク** (dual-time, SLAU, DESmode1)。IC = run_0002 res_20000 をバルクリスケール (M/V→44.12) + 6% 3D 撹乱。切り分け経過: ① KEEP+dt4e-6 → step7 NaN、② SLAU+dt4e-6 → step5 NaN (rms_roOmega が毎 step ×10 で先行発散 = **サブ反復未収束**、制御器は無罪)、③ SLAU+dt1e-6+nSub20+cfl_pseudo1 → **安定** | **2000 step 完走・NaN 無し**。**bodyForceCtrl 検証合格**: M/Mt=1.00000±7e-6, fx 932±117 に整定 (γ=0.02, deadbeat γ=1 は 1/dt ノイズ増幅で不可)。fd_shield ゾーニング動作 (壁近傍 13% RANS)。※ detectNaN 素通りは **config キー位置の罠** (time.deltaT 配下でのみパース) と判明 → パーサ修正済 (252543d) | ref (制御検証) |
| `run_0004_ddes_keep` | **KEEP スタック再テスト** (dt1e-6/nSub20/cfl_pseudo1、他は run_0003 と同一、IC=run_0003 res_2000)。dt4e-6 の NaN は dt 過大が真因 (SLAU も同 dt で死亡) のため復権判定 | **1000 step 完走・NaN 無し**・fx=1058±42・M/Mt=1.00000。**本番は KEEP 構成で確定** (grey-area 検証に必要な低散逸) | ref (KEEP 可否判定) |
| `run_diag_nsub10` / `run_diag_dt2e6` / `run_diag_dt2e6_nsub10` | コストチューニング診断 (各 500 step、IC は直前 run の末尾場): nSub 20→10 / dt 1→2e-6 / 両方 (4 倍速候補) | nsub10: **PASS** 2.35 step/s (2 倍)。dt2e6: **PASS** 1.18 step/s (物理 2 倍, ω 残差最良 0.105)。複合 (4 倍): NaN 無し完走だが **ω 残差に間欠バースト** (0.1→5.5→回復の繰返し) = 縁辺安定 → **不採用** | 破棄予定 |
| `run_0005_ddes_baseline` | **素の DDES baseline 本番**: KEEP スタック + dt2e-6/nSub20/cfl_pseudo1 (診断で ω 残差最良の 2 倍速構成)、bodyForceCtrl γ=0.02、150k step = 62 FT (spin-up ~15 FT + 統計 ~47 FT、実時間 ~35h)。IC = diag 群で ~1.3 FT 発達済みの場。snapshot 2000 step 毎 (0.8 FT, ~75 個)。**判定目標**: 文献の素 DDES と同じズレ方 (再付着点の下流ずれ・せん断層応力欠損) をするか — 合えば実装健全、その後 SLA (Δ̃_ω) 導入で改善幅を見る | (実行中) | active |

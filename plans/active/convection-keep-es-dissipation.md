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
| 3 | TP 単成分 ($s^0(T)$ 多項式=`thermo_s0_mass` 流用) | **実装+検証済** |
| 4 | **多成分** (混合物性 w/S + 混合エントロピー + カーネル内 Y 正規化) | **完了** (プラトーバグ根治済・下記) |

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
- `2026-07-19` — **node 検証 + L3 (WALE 込み粘性 TGV Re=1600) 完了** (残課題の消化):
  - **node L1 (`case/35...run_0023/0024`)**: null-mode は node でも成立 (step0 全 rms 機械ゼロ)。matrix σ=0.05 で減衰するが **1.5桁/400step と cell (4.2桁) より大幅に遅い**。切り分け: 内部ノード 2.7e-5 vs **周期継ぎ目ノード 3.3e-4 (0.6桁) が律速**。node periodic 合併 DOF 周りの散逸の効き (state mirror タイミング等) は**要別調査** (open item)。
  - **node L2 (`case/09...run_0029`)**: KE cost +1.371% = cell (+1.356%) と一致。**レイヤは node 内部で cell 同等**。
  - **L3 (`case/09...run_0030/0031/0032`)**: WALE+Re1600, 32³, t*≈10。σ=0: ピーク 0.066@t*10.1 (過小・遅め=32³ 解像不足)。**σ=0.05 は層流期 (t*<4) に 5-6% 食い LES には過大**。**σ=0.015 が推奨**: 層流期 0.6%・ピーク 0.080@t*8.4・K/K0(10)=0.60。定性ゲート (t*≈9±1 ピーク・滑らか・NaN なし) は 3 run とも通過。
  - **最終推奨 (用途別 σ)**: 低マッハ市松対策・頑健性優先 = σ0.05 (既定)。**WALE 併用の解像 LES = σ0.015**。市松と LES 精度の同時要求では 0.015-0.03 をケースで較正。
  - **残 open**: node 周期継ぎ目の散逸弱化の根治 / Step 3 (TP) / Step 4 (多成分) / 64³ 以上での DNS 定量比較。
- `2026-07-19` — **★訂正: 「node 周期継ぎ目の散逸弱化」は存在しなかった (測定バイアス)**。根治調査の結果:
  - 境界サブセットの A_cb 射影は sum(par)=2≠0 のため **平均圧×パリティ不均衡の DC バイアス** (予測 3.254e-4) を持ち、これを「残存市松」と誤認していた。**共分散射影 (バイアス除去) では step400 で all/interior/boundary = 6.78e-6 に完全一致 = 減衰は一様・継ぎ目弱化なし**。`plot_checkerboard_decay.py` を共分散射影に修正。
  - cell との桁差 (4.15 vs 2.17 桁/400step) は **node の dt が半分 (4.46e-2 vs 8.92e-2)** なだけ。**単位時間あたりは 0.122 vs 0.116 桁/時間で一致 — 散逸レイヤは node で cell と完全同等**。
  - mirror/gather のコード確認: 両方ともステージごとに実行されており健全 (main.cpp:1187, assembleResidual 末尾)。
  - 新規の軽微 open: **node periodic の dt が合併前 half-CV 体積の CFL で不必要に半分** (正しさでなく効率の問題、2倍コスト)。setDT の周期 group 合併体積使用で回復可能 → median-dual plan 側の改善候補。
- `2026-07-19` — **L3 σ 掛引 (0.02/0.03) + Step 3 (TP 単成分) 実装・検証完了**。
  - **L3 掛引** (`case/09...run_0033/0034`): σ=0.02 → K/K0(10)=0.594・層流期 3.1% / σ=0.03 → **0.557 (DNS帯 0.50-0.57 内)**・4.1%。**LES 実用帯 = σ0.02-0.03** (0.02=層流期重視、0.03=終値重視)。32³ では両指標の同時充足は不可 (解像度律速)。
  - **Step 3 (TP)**: 実装前に `tools/verify_entropy_scaling_tp.py` で TP の $w=\frac1T[g-\frac12|u|^2,\mathbf u,-1]$ (Chalot-Hughes-Shakib)・**S 閉形式 = 音響 $\rho/(2\gamma R)$・エントロピー $\rho/c_p(T)$・せん断 $\rho T$**・エントロピー波 $r_E=\frac12 q+e-c_vT$ を数値検証 (PASS, datum 不変も確認)。実装: KEEP 中心 `Itilde` を TP では保存量 roe ベースに (EOS 整合)、scalar は Δroe=保存量ジャンプ+sonic、matrix は TP w/S (thermo 評価 double・`thermoHrefTemp` 焼き込み係数で自動 datum 整合)。多成分は wrapper でエラー (Step 4)。CPG 経路は分岐分離でビット不変。
  - **L1-TP 検証** (`case/35...run_0025/0026/0027`, N2, P=238kPa/T=705K/u=10m/s): pure で null-mode 不変 (TP でも成立)、matrix σ0.05 で **3.9桁/400step 減衰**、scalar も同等。全 NaN なし。
  - 残: Step 4 (多成分: Chalot w + Gouasmi スケーリング + ゼロ濃度種/datum 対策)。
- `2026-07-19` — **Step 4 (多成分) 部分完了 + matrix×多成分の未解決バグを封印**。
  - 設計: 種輸送は forge の scalar 輸送に任せ (散逸込み massflux で連続式と整合)、**5式側の散逸を混合物性で評価** (w に混合エントロピー −ΣY_kR_k ln X_k を含む; Y ln X→0 でゼロ濃度種正則)。凍結組成では TP 1気体に厳密縮退 (python 検証 PASS)。組成ジャンプ面は「Q が SPD なので無条件散逸的」までで全系厳密 ES は主張しない。
  - **検証結果**: pure=null-mode 不変 (多成分でも成立, `case/35...run_0028`)・**scalar (type1) 多成分=クリーン減衰 1.1e-7** (`run_0033`)・組成半割 smoke=NaN なし場有界 (`run_0030`)。
  - **★未解決バグ (matrix type2 × nSpecies≥2)**: 市松が最初 ~50step 正常減衰後 **2.3e-4 プラトー** (`run_0029`)。bisect で確定した事実: (i) ns=1 は新バイナリでも健全 (`run_0031`)、(ii) **Y=[1,0] 縮退 (物理的に純 N2) でも再現** (`run_0032`) → 実在混合の物理でなく ns≥2 コード経路、(iii) scalar は同条件でクリーン (`run_0033`) → 種輸送機構は無罪。潰した容疑: thermo _mix 関数 (質量重み和で Y=[1,0] 等価)・mixent 符号 (Y=[1,0] で 0)・periodic_d の roY ghost コピー (実装あり)・配列サイズ (nCells_all 確保)・組成ドリフト (厳密 [1,0] 維持)・res_0 状態 (7桁一致)。**残る謎: step0 の rms_roUx のみ 0.07% 差 (rms_ro/rms_roe は一致) = 初手から運動量散逸だけ違う**。プラトー減衰率はエントロピー/せん断波レート (|Un|) に整合 → 基底の運動量成分に ns≥2 でのみ入る差異が疑われる。
  - **処置**: wrapper で `keepDissType==2 && nSpecies>1` をエラー封印 (scalar は許可)。多成分 LES は当面 type1 (KE cost は type2 の 2 倍だが正しく動く)。再開時はこの bisect 事実から。
- `2026-07-19` — **★matrix×多成分プラトーバグ根治 (深掘り完遂)・封印解除**。
  - **切り分け過程**: ①カーネル printf で z/y/x 全法線の面をビット比較 → **KEEP_d 散逸は ns=1/ns=2 で完全同一** (矛盾の深化)。②step0 の rms_roUx 差 0.07% は**圧力フラックス巨大相殺 (p·S~9e3 vs 残差~7e-3) 上の atomicAdd 丸め増幅ノイズ**と判明 (rms_roUy/z は両者厳密ゼロ)。③pure 400step: ns=1≡ns=2 (基盤ソルバ無罪)。④**カーネル内で ns=1 経路を強制する in-place bisect → プラトー消滅** = roY 読み〜混合評価に限定。
  - **真因**: ρY_k と ρ が**別カーネル・別 atomicAdd 順で更新**されるため Σρ Y_k ≠ ρ の**共通モードノイズ (~1e-7)** が生じ、未正規化 Y=ρY_k/ρ 経由でエントロピー変数へ。$s^0$ 項が ×~10 (s⁰/T·Δ換算)、混合エントロピー $\ln X$ が対数微分特異で増幅し、**matrix 散逸が弱くしか舐めないエントロピー/せん断方向へ毎ステップ注入** → 市松の弱減衰成分がプラトー化。scalar (type1) は全モードを音響レートで減衰するため同ノイズを隠蔽していた (「scalar は無罪証明」の誤読に注意)。
  - **修正**: matrix ブランチのカーネル内で **Y を正規化 (Σ_k Y_k=1 強制)** し共通モードを除去 (`bc_cell_Y` と同じ流儀)。ゼロ和フォールバック付き。
  - **検証**: 実混合 N2/O2 0.7/0.3 → A_cb 9.77e-4→**4.79e-8 完全減衰** (`case/35...run_0034`)・単成分回帰不変 1.18e-7 (`run_0035`)・組成半割 smoke NaN なし (`run_0036`)。**wrapper の封印を解除し、多成分でも matrix (type2) を第一候補に昇格**。
  - 教訓: (a) 別々に更新される保存量の**比**をカーネルで使うときは正規化必須 (共通モードノイズが対数/エントロピー項で増幅される)。(b) 「強い散逸で問題が見えない」ことは無罪証明にならない。

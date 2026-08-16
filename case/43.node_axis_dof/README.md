# case/43.node_axis_dof — node 軸対称の軸ノード離散化の検証 (case/41 解析モデル流用)

**目的**: node (median-dual) 軸対称で軸ノードをどう扱うべきかを、生産モデル (case/41 wind-tunnel
Md=4 Euler, `run_0057` の nozzle.msh/config をそのまま複製) と離散演算子テストで切り分ける。
[case/41](../41.wind_tunnel_design/) とは別ディレクトリ・別ブランチ (`feature/node-axis-dof`) で作業。
関連 plan: [architecture-node-centroid-value-position.md](../../plans/active/architecture-node-centroid-value-position.md)
(値の位置=ノード座標 / 軸半径 r̄ の分離)、[axisymmetric/implementation.md](../../methods/axisymmetric/implementation.md)。

## 結論 (2026-08-16)

1. **軸半 CV の O(1) 不整合を演算子テストで定量化** ([optest/](optest/)): ρ=const, u=U₀−2ax, v=ar (連続式残差 ≡ 0) を
   1 step 陽解法で評価。現行 (値=ノード、幾何=双対重心 r̄=Δr/4) は軸 CV の質量残差 e=+0.83 (GG)、
   fx=0.5 or LSQ でも +0.5、第一/第二内点行にも −0.17/+0.03 が漏れ、**h に依らず不変** (メッシュ細分で消えない)。
   `nodeValueAtNode: 1` (値の位置をノード座標に、r̄ は r 重み専用 `rEff`) で全行が丸め床 (≤2e-3) に落ちる。
   1 次風上は元々厳密 (中点平均)。fx=0.5 単独・LSQ 単独では直らない (advisor 指摘どおり)。
2. **ノズル (case/41 モデル, Euler)**: `nodeAxisDirichlet: 0` で軸 DOF を解いても真空崩壊は起きない (case/40 の崩壊は再現せず、
   Pt 1 MPa/Md 4 では出ない)。ただし軸ノード M が滑らかに −0.5〜−0.8% 欠損 (偶関数外挿との差 max 0.033)。
   `nodeValueAtNode: 1` で欠損は ≤0.007 (出口∩軸コーナー 1 ノードを除く) に縮み、近軸場の cell 参照との差も 0.0067→0.0031。
3. **軸量の抽出は軸ノード直読でなく偶関数外挿** (`forge_design/evaluate/axis_extract.py`, `axis_curve_node(mode="evenfit")` 既定化):
   設計指標 ‖ΔM‖∞ は直読で Dirichlet 0.222% / 軸DOF 0.888% / van 0.257% と処理法に敏感だが、偶関数外挿では
   0.147 / 0.112 / 0.165% (cell 参照 0.151%) に収れんする。**生産 run_0057 の 0.224% Md は軸コピーのバイアス分 (~0.07% Md) を含む**。
4. **コーナー欠陥の正体はゴーストの回転体積潰れ → dt 崩壊 (解決)**: 値位置=ノードでは境界ノードが境界面上に乗り
   鏡映ゴーストが同位置に退化する。軸∩境界 (r=0) ではゴーストの回転体積が r 床 (1e-20) に潰れ、`setDT` の
   `dx=vol/|S|` が ~1e-20 → 局所 CFL 1e13〜1e14 → `dt_local` ~1e-22 で**入口軸・出口軸の 2 ノードが IC のまま
   完全凍結**していた (ro がビット一致で実証)。ゴーストの r 重みを所有 CV の r̄ にする修正 (`variables.cpp`) で解消
   (凍結 0・cfl max=目標 4.0)。副次的に残差が 2.4→**2.7 桁** (Dirichlet 1.9 桁) に改善。
5. **NS (粘性) は 2 次再構成が壁スリバーで破綻 (未解決・Euler ゲート)**: case/41 NS 生産モデル (y+1.5,
   `wall_first_frac` 4.5e-5, SST) で `nodeValueAtNode: 1` は step 62–64 で発散。**機構は未特定**で、事実と仮説を分けると:
   - 事実 (A/B): `nodeAxisDirichlet: 1` 併用でも同 step 64 (**軸 DOF は無関係**)、1 次 × cfl 1.0 は完走、
     2 次 × cfl 0.2 は step 154 (**CFL は遅らせるだけ**)、段階起動 (soft→2 次) でも step 62、
     粗い壁格子 (`wall_first_frac` 5e-3) では 2 次でも完走、hoop 面積を離散閉性に置換しても同 step 64。
     **GG → `gradLSQ: 2` だけが効き step 64 → 1693 (26 倍延命)** — ただし最終的には落ちる。
   - 事実 (幾何): 壁 CV の「値の位置→面重心」距離は 4.86e-5 → 9.61e-5 m (**1.98 倍**)、壁法線成分は
     3.89e-7 → 7.84e-7 m (第一セル厚 3.47e-7 の 1.1 → 2.3 倍)。発散起点は第一内点の **P=1.11 MPa > 室圧 1.00 MPa**。
   - → 主因は外挿距離そのものより**勾配・再構成の点値整合**側にある公算 (LSQ の 26 倍延命)。NS 対応には
     (a) 勾配・リミッタの点値整合、(b) `wall_dist` の node 基準再計算、(c) ghost 撤廃 (plan §5 Step 1) が要る。
     config に警告を追加済み。
6. **別件で実在する不整合: r 重みメトリックの自由流閉性が壁 CV で崩れている (float32 由来・生産経路も同じ)**。
   一様圧の整合条件 Σ_f r_f·S_f = (0, A_planar) を実測すると、生産 NS メッシュ (y+1.5) の**壁 CV で最大 59%**
   (x 成分 29%、1e-3 超が 3597/54417 CV)、Euler メッシュで 0.34%。合成メッシュの AR スイープ
   ([optest/closure_ar_sweep.py](optest/closure_ar_sweep.py)) では誤差 ≈ 5 ulp × 打ち消し比 (∝ AR) で**壁テーパ非依存**
   → 幾何の lumping でなく **float32 メトリック丸め** (`geom_float = float`)。`nodeValueAtNode` と無関係で
   **生産 Dirichlet 経路にも同じだけある**。一様圧で壁 CV に 59% の偽半径力が立つ。対策候補: メトリックの FP64 化、
   または hoop 面積を離散閉性 Σ_f r_f S_f に置換 (`axisRFloor` 経路に実装済み)。**NS 発散のトリガーではない** (run_0018)。
   なお**保存則 (telescoping) は崩れていない**: 内部面は 1 本の r_f·S_f を両側で符号反転して使うので Σ_i V_i q_i の
   変化は境界フラックスのみ。体積の r̄ をどう取っても保存は保たれ、かつ V = r̄_area·A_planar = ∫r dA は**厳密**なので
   `rEff` は近似ではない。実際に保存を破っているのは `nodeAxisDirichlet` (状態上書き + 残差 0 化) と修正前の
   凍結コーナー CV の方で、軸 DOF 化は保存を回復する。

7. **自由流の回復は「圧力ゲージ」で足りる (倍精度不要)**: 6. の桁落ちは $p$ を悪条件な metric 和に掛けるから起きる。
   対流流束は既に `space.pRef` で $(p-p_{\rm ref})S$ を使うのに **hoop ソースだけ絶対 $p$ のまま**だったので、
   軸対称 × `pRef>0` は**偽半径力 $p_{\rm ref}A$ (1.25e6 m/s²)** を出すバグ状態だった。ソースを同じゲージに直すと、
   壁 CV の偽半径加速度は AR 1250 で **3.4e4 → 2.1e-3 m/s²**、`hoopAreaFromClosure: 1` 併用で **1.5e-4**
   ([optest/freestream_test.py](optest/freestream_test.py))。生産 Euler では場が 9e-5 動くだけ・収束同等 (run_0021)、
   既定 OFF はビット不変 (run_0022)。**ただし NS の `nodeValueAtNode` 発散はこれでも直らない** (step 64→65, run_0023/0024)
   = 5. の未解決は metric ノイズが原因ではない。plan [axisymmetric-freestream-hoop-gauge](../../plans/active/axisymmetric-freestream-hoop-gauge.md)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| `run_0001_ref41_0057_dirichlet` | case/41 `run_0057_axismach_srcfix_n2000` の入力+最終場のコピー (生産: node, `nodeAxisDirichlet: 1`, Euler cplus n2000)。比較の基準 | 軸−偶関数外挿 max 0.0030、‖ΔM‖∞ 直読 0.222% / 外挿 0.147% Md。`res_12000.h5` | ref |
| `run_0002_axisdof_euler` | 0001 と同一入力で `nodeAxisDirichlet: 0` (軸ノードを DOF として解く、値位置=双対重心) | NaN 0・STEADY・plateau (生産と同型)。真空崩壊なし (ρ_axis/ρ_row1 0.989–1.018)。軸 M 欠損 max 0.033 (0.5–0.8%)、‖ΔM‖∞ 直読 0.888% | active (**軸 DOF 対照**) |
| `run_0003_axisdof_van_euler` | 0002 + `nodeValueAtNode: 1` (試作: 値位置=ノード座標, r̄→rEff) | 軸−外挿 ≤0.0069 (出口コーナー除く)、cell 参照との差 0.0031、‖ΔM‖∞ 直読 0.257% / 外挿 0.165%。出口∩軸コーナー 1 ノード M 3.892 | active (**恒久策試作の根拠**) |
| `run_0004_axisdof_van_lsq` | 0003 + `gradLSQ: 2` | 0003 と 3 桁一致 (コーナー含む) = 勾配法は無関係 | active |
| `run_0010_van_ghostfix` | 0003 + **ゴースト r 重み修正** (境界ゴーストの回転半径を所有 CV の r̄ に) | 凍結ノード 0・`dt_local` 正常 (cfl max 4.0)・残差 2.7 桁 (0003 の 2.4 / Dirichlet 1.9 より良い)。軸−外挿 ≤0.0069、bulk は 0003 と同じ | active (**Euler 生産候補**) |
| `run_0007_van_1st` | 0003 + `convMethod: 0` (コーナー仮説の A/B。結果的にコーナーは凍結由来と判明) | 1 次解に収束。コーナー診断の副産物 | active |
| `run_0008_ns_dirichlet` / `run_0009_ns_van` | **NS スモーク** (case/41 `run_0071_ns_v2` の mesh/config/収束場, 4000 step): 生産 Dirichlet vs van | 0008 完走 (NaN 0)、**0009 は step 64 発散** (`res_nan_64.h5`) | active (**NS 未対応の根拠**) |
| `run_0011_ns_van_axisdir` | van + `nodeAxisDirichlet: 1` (軸 DOF を切り離す) | **同じ step 64 で発散** = 軸 DOF は無関係 | active (切り分け) |
| `run_0012_ns_van_soft` / `run_0013_ns_van_staged` | van soft 起動 (1 次 cfl 0.2, 3000 step) → その場から 2 次 cfl 1.0 | 0012 完走・場は Dirichlet と同等 (P/T/μt 一致水準)、**0013 は step 62 発散** = 段階起動では回避不可 | active (切り分け) |
| `run_0014_ns_van_1st_cfl1` / `run_0015_ns_van_2nd_cfl02` | 1 次 × cfl 1.0 / 2 次 × cfl 0.2 | **0014 完走・0015 step 154 発散** = 原因は 2 次再構成 (CFL は遅らせるだけ) | active (**NS 原因特定**) |
| `run_0016_ns_coarse_van` / `run_0017_ns_coarse_dirichlet` | 粗い壁格子 NS (`run_0069_ns_v1_coarse`, y+~50, 2 次 cfl 1.0, 3000 step) | **両者とも完走・NaN 0** = 破綻は y+≈1 の高 AR 壁スリバー固有 | active (**境界条件の切り分け**) |
| `run_0005_cell_ref` | 同 nozzle.msh を cell 変換 (品質 PASS AR 9.9/skew 0.41)、`interp_field` で 0057 場から段階起動 (soft 1次 cfl0.5 3000 → 本段) | NaN 0・STEADY・plateau 0.8 桁 (cell の既知床)。セル中心列の偶関数外挿 ‖ΔM‖∞ 0.151% Md — node 外挿の独立参照 | active (**参照**) |
| `run_0006_regress_dirichlet` | 0001 と同一入力 (生産 `nodeAxisDirichlet: 1`) を**パッチ後バイナリ**で再実行 (既定経路のビット不変確認) | 最終場は 0001 と最大相対差 4e-6 (ro/roUx/roe/P; node の run 間非決定性レベル) = 既定経路は不変 | active (回帰) |
| `run_0018_ns_van_closure` / `run_0019_ns_dir_closure` / `run_0020_ns_van_lsq` | NS 発散の切り分け続き: van + `axisRFloor 1e-12` (hoop 面積を離散閉性に) / 同 Dirichlet 対照 / van + `gradLSQ: 2` | 0018 **同 step 64 発散** (閉性は無関係)、0019 完走、0020 **step 1693 まで延命** (GG→LSQ で 26 倍) | active (**NS 切り分け**) |
| `run_0021_euler_closure` / `run_0022_regress2` | 生産 Euler + `hoopAreaFromClosure: 1` / 既定 OFF 回帰 | 0021: 場の差 9e-5・残差 2.7 桁 (同等)。0022: 旧 build と 1.8e-6 (node 非決定性レベル) = 既定はビット不変 | active |
| `run_0023_ns_van_pref` / `run_0024_ns_van_pref_closure` / `run_0025_ns_dir_pref` | NS van + `pRef 1e6` / + 閉性面積 / Dirichlet 対照 | **0023・0024 とも step 65 発散** (pRef でも閉性でも NS は直らない)、0025 完走 | active (**NS の原因から metric を除外**) |
| `run_0026_ns_van_trace` | van NS の発散過渡を 8 step 毎に出力 (70 step) | **壁ノード↔第一〜第三内点で壁法線方向の奇偶 (P/Uy) モードが step ~24 から指数成長** (e-fold ~8 step): P が 9.954/9.961/9.954/9.960 → 9.84/1.004e6/9.94、Uy ∓3.7/+1.5。起点は室 (M~1e-3, P 1 MPa) の壁 | active (**NS 発散モードの同定**) |
| `run_0027_ns_van_lmp2` / `run_0028_ns_van_lmp2_lsq` | van + `lowMachPrecond: 2` (/ + LSQ) — 低マッハ市松への常套手段 | **両方 step 54 で発散 (悪化)** = 低マッハ前処理では直らない | active (常套手段の除外) |
| `optest/fs_*` | 自由流テスト (一様圧・静止, AR 0.5–1250) | 偽半径加速度 3.4e4 → 2.1e-3 (pRef) → 1.5e-4 (pRef+閉性) m/s² | active (**ゲージ整合の根拠**) |
| `optest/closure_*` | r 重み閉性の AR スイープ (幾何のみ・solve 不要) | 誤差 ≈ 5 ulp × 打ち消し比、taper 非依存 → float32 丸め。AR 1250 で 2.5% | active |
| `optest/h{0.02,0.01,0.005}_*` | 離散連続式演算子テスト (`optest/run_optest.py`, 直管 1×0.2 m 一様 quad, 1 step 陽解法, `limiter: 0` は線形場再構成のため意図的) | `optest/optest_results.json`。要約は上記 1. | active (**演算子テスト**) |

解析スクリプト: `analyze_axis_row.py` (軸行 vs 第一内点), `compare_axis_runs.py` (→ `axis_runs_compare.png`)。

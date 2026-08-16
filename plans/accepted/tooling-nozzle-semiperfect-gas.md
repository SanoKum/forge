# axis-Mach チェーンの semi-perfect gas 化 (NASA-9/CEA, frozen 組成) と イソブタン風洞スタディ

## メタ

- **area**: `tooling / optimization`
- **status**: `done`  <!-- 2026-08-17 起票・同日 R/L スタディ完了 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md)
  - [`methods/thermophysics.md`](../../methods/thermophysics.md) (forge 側 TP・NASA-9 の正本)
- **related_plans**:
  - 親: [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md)
  - 上位: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) (多目的最適化への接続)
- **created**: `2026-08-17`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的 (ユーザ指示 2026-08-17)

MOC に NASA-CEA の物性を **semi-perfect (thermally perfect, frozen 組成)** として載せ、
入口全圧・全温・組成をインプットに設計できるようにする。案件: **イソブタン燃焼ガス風洞**
($P_t$=5 MPa, $T_t$=1000 K, C₄H₁₀/空気 φ=0.9 生成物, 出口径 1 m, 入口配管径 1 m, $M_d$=4)
で **R (スロート曲率) と L (入口側 $L_U$・出口側 $L_c$) をスタディ**する。
その先は多目的最適化 (軸 M 目標や長さを振る手段としての形状生成)。

## 2. 設計

- **`forge_design/gas/semiperfect.py`**: forge 内蔵 DB (`thermo_d.cu::builtinDB`, CEA
  McBride–Gordon 2002) と**同一の NASA-9 係数**を Python 側に持ち、$\nu(M)$・$A/A^*$・
  $\rho V$・$T(M)$・$\gamma(T)$ を等エントロピー膨張の $T$ パラメトライズで数値積分して
  テーブル化。**MOC は $\nu\leftrightarrow M$ 変換と流束密度だけ差し替えれば thermally
  perfect でも成立** (適合条件 $\theta\pm\nu$ 一定は不変)。
- **MOC カーネル規約**: `pm_nu`/`pm_mach`/`pm_mach_vec`/`_mass_flux_density`/
  `area_ratio_isentropic` は「γ (float) の位置にガスモデルを渡せる」(`_is_gas`)。
  Hall 遷音速級数は定数 γ 前提なので**スロート局所 γ\*** を渡す。
- **forge 側**: frozen 組成を `mixture_pseudo_species` で**単一擬似種 (NASA-9 の質量分率
  加重、厳密)** にまとめ `speciesDBFile` で渡し `thermalMethod: 2` の 1 種 TP で回す
  (多成分 implicit の既知不安定性 [[wys-tp-divergence-is-cold-not-multispecies]] を回避)。
  IC は TP の $e(T)=h(T)-RT$ (絶対基準) で `roe` を構成。
- **probdef**: `gas.model: semiperfect` + `gas.species {名: 質量分率}`。

## 3. 実測・知見 (2026-08-17)

**(a) 物性**: N2 のみで CPG(1.4) に $\nu$ 2e-4 (M=2)〜4e-3 (M=3) — 差は NASA-9 の 200 K
未満外挿 (forge 本体も同じ) で真の差。一定 cp 擬似種では 1e-7 で機械精度退化 (実装は正しい)。
イソブタン φ=0.9 生成物 (CO₂ 16.8 / H₂O 8.6 / O₂ 2.2 / N₂ 72.4 wt%): $R$=292 J/kgK、
**γ は 1.31 (スロート, T\*=868 K) → 1.38 (M4, 266 K)** と変化。$A/A^*(M4)$=13.105。
**一定 γ=γ\* の CPG は出口半径を +8% 誤る** (15.29 vs 13.11) — semi-perfect が必須。

**(b) 設計チェーン**: 同じ Hall アンカー・law で semi-perfect MOC は 0.2 s、$r_F$ が
1D 理論と −0.19% (CPG 系列と同精度)、$\kappa_0 R$=0.996。$r_t$=0.138 m で出口径 0.998 m。

**(c) 設計スイープ (36 点, 0.4 s/点)**: R∈{1.5,2,3,5}×L_U∈{4,6,9}×L_c∈{max,8,6}:
- **R=5 は全滅** (壁非単調/リンギング — $M''_A$ 過大で単一 quintic では壁にならない)
- **R=1.5–3 は L_c=max/8 が有効、L_c=6 は R=1.5 のみ** (R=2,3 はリンギングゲート =
  MOC 壁の要求曲率が急、内部衝撃の前兆)
- **L_U=9 は R=1.5 で μ>20 (非単調収縮)**。R≥2 なら可
- 全長 3.4–4.8 m。$\theta_{max}$ 19.0°(R3) → 19.8°(R1.5)

**(d) CFD 投入で踏んだ罠 (3 つ、全て解決)**:
1. `interp_field.py` が res の primitives から `roe=P/(γ−1)+½ρu²` (CPG) で再構成 → TP では
   $T$ が 3000 K に跳ぶ。**res に保存量 `roe` があればそのまま使う**よう修正 (CPG は不変)。
2. `nodeAxisDirichlet=1` は forge が TP で弾く → `axisymMethod: 1` に逃げたら本段 step ~1000
   で**出口軸コーナーから発散**。**CPG でも method 1 は同じ発散** (case/41 run_0076 で切り分け、
   TP は無罪) → 軸対称高度化 (別セッション) の対象。回避 = method 0 + Dirichlet 0
   (+ 必要なら `axisRFloor`)。
3. 出口背圧: 1000 Pa は 5 MPa 系では出口静圧 28 kPa と桁違い → 整合圧に。
4. TP は本段 cfl4 で爆発 → mid 段 (2次 cfl1) を挟み cfl2 (`evaluate.cfl_main`/`mid_stage`)。
5. **TP × node 軸対称 (method 0, Dirichlet 0, axisRFloor 2.8 mm) でも soft 段 step 271 で
   膨張部中域 (x≈5–6) の軸近傍から発散**。**同一壁・同一メッシュ・同一段階起動を CPG(γ\*) で
   回すと完走・良好** (run_0002: ‖ΔM‖∞ 0.284% $M_d$、出口 ε_M 0.013%/ε_θ 0.0096°、
   品質 PASS、軸 M は 4000 step 間 8.8e-6 で凍結)。→ **TP × node 軸対称は forge 側の穴**
   (`nodeAxisDirichlet` が TP で使えず、代替 [`axisymMethod:1` は CPG でも出口軸コーナーで
   発散 / `axisRFloor` は真空化は抑えるが中域軸で破綻] が無い)。**別セッションの軸対称
   高度化への申し送り事項**。本計画では **設計 = semi-perfect、CFD 検証 = CPG(γ\*)** で
   進める (`evaluate.cfd_gas: cpg`)。設計点間の相対比較 (衝撃・一様性の傾向) には十分。
   TP CFD はソルバ対応後の再検証項目。
6. **`cfd_gas: cpg` は設計も CPG に落とす**: 最初「設計 semi-perfect・CFD だけ CPG」にしたら
   壁 (A/A*=13.1) と CFD (CPG 1.309 → 15.3) が食い違い、9 点全部で出口コア M 3.86・
   ΔM 4.5% (run_0003–0011 初回、破棄)。**設計と CFD の熱力学は必ず一致**。CPG 検証では
   出口径が 1.078 m になる (semi-perfect 正本 0.998 m) ことを明記。

## 4. 完了条件

- [x] ガスモデル + MOC 配線 + 単体テスト (`run_gas_tests.py` 18 項目 ALL PASS)
- [x] forge TP 用擬似種 + IC + config 注入
- [x] イソブタン基準 run の CFD 完走 (CPG(γ*) 経路, run_0002)・TP は forge 側穴として申し送り
- [x] R/L スタディの CFD (M4.2, 8 点, run_0012–0019) + 結果ページ
- [x] methods 更新・README・commit

## 3.5 R/L スタディ結果 (M4.2、ユーザ訂正後の正本)

**設計スイープ (semi-perfect, 36 点)**: R=1.5–3 × L_c∈{max, 8} が有効、L_c=6 は全 R で棄却
(M4.2 は膨張が大きく短ノズルは急曲率)、R=5 は全滅、R1.5×L_U9 は μ>20。全長 3.65–4.73 m、
$r_t$=0.1258 m (スロート径 252 mm)。

**CFD (CPG γ\*=1.309、設計・CFD 熱力学一致、8 点。全 run 品質 PASS・NaN 0・軸 M 凍結 ~1e-5)**:

| R | L_U | L_c | ‖ΔM‖∞ [% M_d] | 出口 ε_M | ε_θ |
| --- | --- | --- | --- | --- | --- |
| 1.5 | 6 | max | 0.328 | 0.013% | 0.018° |
| 2 | 6 | max | 0.306 | 0.017% | 0.013° |
| **3** | 6 | max | **0.256** | 0.025% | **0.003°** |
| 1.5 | 6 | 8 | 0.436 | 0.021% | 0.021° |
| 2 | 6 | 8 | 0.523 ✗ | 0.020% | 0.021° |
| 3 | 6 | 8 | 0.681 ✗ | 0.023% | 0.024° |
| 2 | 4 | max | 0.296 | **0.061%** | 0.013° |
| 2 | 9 | max | 0.302 | **0.010%** | 0.018° |

- **R**: L_c=max では R↑で誤差↓ (0.33→0.26%)。L_c=8 では逆転 (R↑で悪化、R≥2 でゲート外) —
  短ノズルでは大 R の小さい $M''_A$ が尾部減速を急峻化。R=5 は設計不能。
- **L_c**: max (≈10.8 r_t) が推奨。8 は全長 −0.35 m で誤差 2 倍。M4.2 は M4 より短縮余地が狭い。
- **L_U**: 軸 M に効かず**出口一様性に効く** (L_U=4 で ε_M 3–6 倍)。R=1.5 では μ≤20 から L_U≤8.9。
- **推奨ベースライン R=3, L_U=6, L_c=max**: 全長 4.35 m、スロート 252 mm、出口 1.0 m。
  MOO の探索域目安 R∈[2,3]・L_U∈[6,9]・L_c∈[8,max]。
- 結果ページ: https://claude.ai/code/artifact/07b52b7f-54b1-4cc8-97a5-4ee9af380b2d

## 5. 変更ログ

- `2026-08-17` — 起票・実装・設計スイープ完了。CFD は罠を潰して再投入。
- `2026-08-17` — M_d を 4.0→4.2 に訂正 (ユーザ)、M4.2 で設計スイープ・CFD 8 点完了、推奨ベースライン確定。**完了**。
- `2026-08-17` — **新 binary (`ea08bcbe`, nodeAxisDirichlet 撤去) で TP CFD を再挑戦 → 再び発散**
  (run_0020: soft 段 step 190、x≈7 の軸上で P が床へ)。旧 binary の症状 (step 271, x≈5–6) と同性質で、
  **TP × node 軸対称の発散は nodeAxisDirichlet と無関係、forge の TP 経路 (EOS 温度反転か軸
  ソース項の TP 版) 固有**と確定。同一メッシュ・IC の CPG は新旧 binary とも完走。
  申し送りを更新: 「TP 単一種 (thermalMethod 2, species [MIX]) + node 軸対称 Euler が膨張部の
  軸で発散する。CPG は同条件で完走。case/42 run_0001 (旧) / run_0020 (新) が再現ケース」。
- `2026-08-17` — **TP × node 軸対称の発散原因を特定・解消**。切り分け (run_0021–0025、soft 段設定で 1 要素ずつ):
  一定 cp 種 ○ / 200 K 未満 cp 凍結 ✗ / 陽解法 ○ / 軸ホップ Jacobian の一般 EOS 化 ✗ /
  **`thermoHrefTemp: 298.15` ○**。→ **原因 = 陰解法 Jacobian の $\chi_{\rm eos}=c^2-\kappa h$ に
  絶対基準の $h$ (生成エンタルピー込み) が入り桁違いになること**。sensible datum で解消。
  runner の TP 経路に `thermoHrefTemp` を既定化 (`evaluate.thermo_href_temp`, 既定 298.15) し
  IC も同 datum。**フル TP CFD (run_0026, R3/L_U6/L_c=max): ‖ΔM‖∞ 0.116% $M_d$・ε_M 0.014%・
  ε_θ 0.0036°・出口径 0.999 m** — 設計と CFD の熱力学が完全整合し CPG(γ\*) 検証 (0.256%) の
  半分以下。継続 (run_0027, +24000) で軸 M の変動は 1.3e-3→1.0e-3→3.7e-4 と単調減少 (収束途上、
  残差床 3e-2 は TP EOS 反転の float 床と推定)。軸ホップ Jacobian の一般 EOS 化はソルバに残置
  (CPG は解不変、差 2e-5 = χ_eos の float 丸め分)。申し送りは「解消済み、thermoHrefTemp 必須」に更新。
- `2026-08-17` — **TP 残差床の正体 = TP 入口 BC の緩和振動**: run_0026/0027 は発散せず (NaN 0、24000 step) 内部場も
  凍る (軸 M 変動 4e-4) が、`check_convergence` は FAIL (rms_ro 3e-2 プラトー、CPG は 1e-9)。入口面で
  **P/Pt=1.03 (全圧超過)** が snapshot 間 1.3e-3 で振動 = TP `inlet_Pressure` (静圧参照/速度参照 rf=0.5
  緩和、コード内コメントの「純静圧は振動・純速度は P>Pt 反射」の病理) が大径低速入口 (M≈0.03) で顕在化。
  `inlet_Pressure_dir` (run_0028) は step 9059 で**発散**。**forge 側 TP 入口 BC への申し送り** (CPG 経路
  との非対称)。「収束した」とは報告しない。ノズル内部の設計評価値 (0.116% Md 等) は場が凍っているので
  相対比較には有効だが VERDICT は NOT CONVERGED と併記する。
  配管長 0.5→2 r_t (run_0029) でも不変 → 入口 BC 自体の性質。設計側の手は尽くした。
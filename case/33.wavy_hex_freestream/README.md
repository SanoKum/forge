# 33. wavy hex free-stream 検証

正弦波で歪ませた非直交 hex (10^3) で forge の free-stream(一様場)保存を検証する。
非構造/非直交メッシュで forge が発散する原因が「free-stream 保存崩れ(metric誤差×流束)」か
「別機構(エネルギー式等)」かを切り分ける。

- mesh/wavy.h5         : forge 入力メッシュ(閉性 1.28e-9・体積厳密に検証済み, 純幾何は正しい。
  品質: AR max 1.5 / skew max 0.412 = PASS)
- mesh/wavy_fluent.h5  : 合成 Fluent 形式の元データ
- 生成: /tmp/make_synth_fluent.py の gen_wavy (内部ノードを amp*h*sin で歪ませた構造hex。
  スクリプトは散逸済み。境界面は平坦で、xmin/xmax は (y,z) 面重心が厳密一致=periodic 化可能)
- mesh/wavy_node.msh   : **node (median-dual) 用の gmsh 4.1 復元メッシュ**。wavy.h5 と同一
  節点・同一 hex 結線 (生成: `mesh/make_wavy_node_msh.py`。xmin/xmax の (y,z) 節点も厳密一致
  = node periodic 可)。各 run 内で `convertGmshToForge` (solverConfig `discretization: node`)
  で median-dual HDF5 化する (closure 5.4e-7 正規化)

検証方針(別LLM助言 + 自前): 一様場の forge 残差を 1ステップ逆算し、
F(U_c)·r_c (r_c=Σ dir·sv = セル幾何閉じ誤差) と成分ごと・符号付きで突き合わせる
(詳細は [VERIFICATION.md](VERIFICATION.md))。

- U=0: 一致(運動量=p·r, 質量=0, エネルギー=0) → 圧力項が主原因 → `space.pRef` で根治
- U∞≠0: 質量/運動量/エネルギー残差比が u・H に一致 → 移流項 (ρu, ρu², ρuH)·r_c が主原因
  → `space.roRef`+`space.uRef` (KEEP 移流基準差分) で根治
  ([plan §8](../../plans/active/convection-freestream-preserving-flux.md))

## 計算 run 一覧

U=0 系は全面 slip、U∞≠0 系 (run_0007 以降の xperiodic) は x 面 periodic + 側面 slip。
IC はいずれも一様場 (ρ=1.177, P=101325; M0.1 は u_x=34.716 m/s)。

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0003_pref_long` | SLAU + pRef=101325 一様静止 長時間 | 機械精度で 300step 安定 (pRef 根治の SLAU 検証) | ref (SLAU pRef 検証) |
| `run_0004_keep_pref0` | **KEEP** (pRef なし) 一様静止 | 運動量残差 7.2e-5 → **step11 発散** = float32 free-stream 崩壊を KEEP でも実証 | ref (KEEP 症状記録) |
| `run_0005_keep_pref` | **KEEP + pRef=101325** | 運動量残差 **3.96e-12 (機械精度)**・300step 安定。SLAU と同水準の根治 | ref (KEEP pRef 検証) |
| `run_0006_keep_dissmat_pref` | KEEP + **matrix ES 散逸 (keepDissType=2)** + pRef | 4.07e-12・安定。一様場で Δw=0 → 散逸レイヤは free-stream を厳密に保つ | ref (散逸レイヤ無害確認) |
| `run_0007_keep_pref_u0_xperiodic` | KEEP+pRef, U=0, **x-periodic 化対照** | run_0005 と同一挙動 (roUx 4.06e-12/300step) = periodic 化は無害 | ref (periodic 対照) |
| `run_0008_keep_pref_uxM01_xperiodic` | KEEP+pRef, **U∞=M0.1** (移流桁落ち実測) | step0 rms ro/roUx/roe = 3.7e-8/1.3e-6/**1.3e-2** (残差比=u,H で移流項閉じ誤差と定量一致) → **step14 発散** | ref (移流桁落ち実測) |
| `run_0009_keep_dissmat_uxM01_xperiodic` | 同 + matrix ES σ=0.05 | **step13 発散** = 散逸レイヤでは救えない | ref (散逸で不可の記録) |
| `run_0010_slau_pref_uxM01_xperiodic` | SLAU+pRef, U∞=M0.1 | **step6 発散** = upwind でも救えない (全スキーム共通のボトルネック) | ref (SLAU 対照) |
| `run_0011_keep_advgauge_uxM01_xperiodic` | **KEEP + 移流基準差分 (roRef=1.177, uRef=[34.716,0,0])** U∞=M0.1 | step0 rms 1.7e-14/8.1e-12/5.6e-9 = **U=0 と同水準の機械精度**・300step 安定 (エネルギー 6.5 桁改善) | ref (advGauge 根治検証) |
| `run_0012_keep_u0_regr_newbin` | 回帰: advGauge 実装後バイナリで run_0007 再実行 (roRef=0) | run_0007 と atomicAdd ノイズ床 (4e-12) 内で一致 = **roRef=0 経路に回帰なし** | ref (回帰エビデンス) |
| `run_0013_keep_advgauge_dissmat_uxM01` | advGauge + matrix ES σ=0.05 (LES 想定構成) | step0 rms 1.4e-12/2.5e-9/3.3e-7・300step 安定 = ゲージと散逸レイヤは併用可 | ref (併用確認) |
| `run_0014_gaugecurve_sweep/` | **非一様場でのゲージ効果カーブ** (M0.1+速度擾乱 a=1e-4..1e-1 × ゲージ有無の step0 残差; サブ run 群) + 発散真因の切り分け dbg 群 | ①擾乱下の step0 残差はゲージ有無で一致 (物理∝a が支配、OFF ノイズ床 1.3e-2 との交差は a≈5e-7)。図 `gauge_decay_curve.png`。②dbg 群で **run_0008-0010 の「発散」は unsteady:0 (定常局所dt) との交絡**と判明 (下注) | ref (ゲージ効果カーブ+真因切り分け) |
| `run_0015_advgauge_urefm20_uxM01` | **uRef 20% ミスマッチ耐性**: uRef=27.77 (真値 34.72) の一様 M0.1 | step0 roUx 4.5e-7・300step 安定 (定常モードでも)。一様場なら偏差も空間定数 → telescoping が生き、ノイズは偏差比例 (~0.2×) に留まる | ref (ミスマッチ耐性) |
| `run_0016_advgauge_urefm20_dissmat` | 同 + matrix ES σ=0.05 | 同様に安定 | ref (ミスマッチ+散逸) |
| `run_0017_node_advgauge_uxM01_xperiodic` | **node (median-dual)** + advGauge U∞=M0.1 (境界半割面ゲージ plan §8.5 の検証) | step0 rms 4.5e-14/7.8e-12/1.9e-9 = **cell run_0011 と同水準の機械精度**・300step 安定・最終場は厳密一様 (min=max) | ref (node advGauge 根治検証) |
| `run_0018_node_nogauge_uxM01_xperiodic` | node ゲージ無し対照 (roRef 無し) | step0 roe rms **9.9e-3** (cell run_0008 と同水準のノイズ床)・定常局所 dt で step15 発散 = node でも移流桁落ちを実証 | ref (node ノイズ床対照) |
| `run_0019_cell_advgauge_regr_newbin` | 回帰: node 対応後バイナリで run_0011 (cell) 再実行 | step0 rms 1.7e-14/8.2e-12/5.6e-9・末尾 1.37e-6 = run_0011 と一致 (cell 経路無影響) | ref (回帰エビデンス) |
| `run_0020_node_advgauge_dissmat_uxM01` | node advGauge + matrix ES σ=0.05 (LES 想定構成) | step0 roe 1.5e-7・300step 安定 = cell run_0013 と同水準 (散逸レイヤ併用可) | ref (node 併用確認) |

**重要な訂正 (2026-07-19, run_0014 の dbg 群で判明)**: run_0008-0010 の「step6-14 発散」は
**`unsteady: 0` (定常・局所 dt) モードとの交絡**。同じゲージ無し M0.1 一様流を `unsteady: 1`
(物理 dt) で回すと発散せず、残差はノイズ床 (roe 1.3e-2) に恒久プラトーで 300step 安定
(`run_0014_gaugecurve_sweep/dbg_uniform_nogauge_unsteady`)。さらに物理擾乱 a=1e-2 を与えると
**ノイズと無関係に** 定常モードは全スキーム (SLAU CFL0.1 でも)・全 BC (全面 periodic でも)・
**Cartesian (case/35 run_0038 steady_variant) でも** step5-11 で発散する = 「閉じた周期系の
音響過渡×定常局所 dt」というソルバ一般の限界であり、wavy メッシュ・periodic 実装・ゲージとは
無関係。**step0 残差 (ノイズ種) の測定と帰属・ゲージによる機械精度化はこの訂正の影響を受けない**。
ゲージの実利は (a) ノイズ種を消して定常局所 dt モードでも U∞≠0 を回せるようにする
(run_0011/0013 は定常モードで安定)、(b) 非定常 LES では物理の下に埋まる恒久ランダム強制
(roe 残差 1.3e-2 相当) を機械精度まで除去する、の 2 点。

(run_0001/0002 は SLAU 初期切り分け run。結果は [VERIFICATION.md](VERIFICATION.md) に集約済みでディレクトリは削除済み)

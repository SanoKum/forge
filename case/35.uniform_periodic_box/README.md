# 35. uniform periodic box (Taylor-Green)

全 6 面 Cartesian periodic の 3D 構造化 hex box (32³, `mesh/tg.msh`)。一様流/Taylor-Green 初期値で
周期境界・3D 経路の最小検証に使う。3D median-dual (M4, [`discretization-median-dual-3d`](../../plans/active/discretization-median-dual-3d.md))
の双対メッシュ生成 (体積総和・closure) 検証ケースでもある。

## 計算 run 一覧

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_uniform` | cell モード一様流ベースライン (convMethod=0) | `res_*.h5`, `residual_history.csv` | ref |
| `run_0002_node_dualcheck` | node モード (`discretization: node`) で 3D median-dual 生成の幾何検証 | 体積 dual=248.05=primal (relErr 3.5e-8)・closure normalized 2.5e-6・負体積 0 (`convert.log`)。solver run は未実施 | active |
| `run_0003_node_periodic_smoke` | node 3D で全面 periodic の forge smoke (20 step, §4.5 実装前) | **step 4 で DIVERGED**。半割面を主面接続する旧挙動は幾何が誤りで発散 → §4.5 実装の動機 | 破棄予定 |
| `run_0004_node_periodic_freestream` | **§4.5 実装後**: 構造化 periodic 立方体(8³ hex, 全方向 periodic)+ 一様 +x 流で free-stream 保存検証 | **PASS**: 全 50 step で残差**機械ゼロ**。union-find `217 slave merged`=729−512(triperiodic 一意 CV)で完全一致。`conv main planes periodic 0`(移流から periodic 除外) | ref (検証エビデンス) |
| `run_0005_node_periodic_tg` | node 全面 periodic + Taylor-Green IC (M0.1) explicit RK3 | step4 DIVERGED。union-find `3169 merged`=35937−32768 正。発散は periodic 起因でなく node+低マッハ+explicit 既知不安定(run_0006 で実証) | 破棄予定 |
| `run_0006_node_slip_tg_isolation` | **切り分け**: run_0005 の periodic を slip に置換(periodic 無し)、他同一 | step4 で**同様に DIVERGED** → TG 発散は periodic 実装でなく node+低マッハ+explicit RK3 の不安定と確定(periodic 無罪) | ref (切り分けエビデンス) |
| `run_0007_node_periodic_tg_lmp2` | run_0005 + lowMachPrecond=2 (explicit) | step1 DIVERGED。explicit では低マッハ安定化せず(node 低マッハは implicit 必須)。非自明 periodic は implicit(§4.5.7)待ち | 破棄予定 |
| `run_0008_node_periodic_tg_implicit` | node 全面 periodic + TG + implicit block-DPLUR(lowMachPrecond=2, cfl_pseudo=1, nStepInner=5) | step64 で **NaN 発散**(explicit step4 より大幅に粘るが発散) | 破棄予定 |
| `run_0009_node_slip_tg_implicit_iso` | **切り分け**: run_0008 の periodic を slip に置換、他同一 | step99 まで **NaN なし**(plateau/roe rising だが発散せず)。**implicit では periodic だけ NaN 発散** → §4.5.7 Jacobian 行 fold が implicit periodic 安定化に必要と確定(RHS gather だけでは LHS 不整合) | ref (§4.5.7 必要性の根拠) |
| `run_0010_node_periodic_tg_implicit_mirror` | **§4.5.7 実装後**: run_0008 + sweep ごと dq ミラー(slave=master) | **step100 完走・NaN なし・場健全**(ρ 0.83–1.02, P/T/roe 有界, 渦発達)。run_0009(slip)と同等プラトーになり periodic 固有発散消失。**implicit periodic 安定化達成** | ref (§4.5.7 検証エビデンス) |
| `run_0011_node_periodic_tg_gradgather` | **勾配 periodic gather 実装後**: run_0010 同設定で勾配 gather 有効 | step100 安定・場健全(ρ 0.825–1.017)で run_0010 と実質不変、seam viscous が片側→両側へ微修正(回帰なし)。free-stream も機械ゼロ維持 | ref (勾配 gather 検証) |
| `run_0012_gradcheck` | **勾配 gather 定量検証**: TG 解析勾配 ∂ux/∂x=M0cosx cosy cosz と照合 | 境界(合併)勾配 mean\|err\|=5.3e-4(\|g\|max 0.398 vs 解析0.4)で内部5.8e-4と同等、master==slave一致(2.9e-7)。periodic境界勾配が内部同等に正しいと確認 | ref (勾配検証エビデンス) |
| `run_0014_cell_keep_cbd_pure` | **L1a/b 市松診断** (ラダー L1, [survey §10](../../notes/investigations/convection-central-scheme-oscillation-control.md)): cell 32³, 一様 M0.1 流+圧力のみ市松摂動 ε=1e-3 (`make_checkerboard_ic.py`), pure KEEP (keepDissType=0) | **null-mode 実証**: 全 rms 残差が厳密ゼロ (中心平均が市松を相殺し原理的に見えない)・A_cb=1e-3 が 400step 後も**不変** | ref (L1a null-mode エビデンス) |
| `run_0015_cell_keep_cbd_diss` | 同 IC + **keepDissType=1 (σ=0.15, lowMachPrecond=1)** | rms_roe のみ非ゼロ (散逸が Δp 経由でのみ市松を見る=設計通り)。**A_cb 1e-3→8.0e-10 (6桁減衰/400step)** | ref (L1b 減衰エビデンス) |
| `run_0016_cell_slau_cbd_ref` | 同 IC + SLAU (対照) | SLAU は市松を ro/roUx/roe 全部で見る (rms 即非ゼロ)。A_cb 1e-3→1.2e-8 | ref (対照) |
| `run_0017_cell_keep_cbd_diss005` | 同 IC + keepDissType=1 **σ=0.05** (較正) | A_cb 1e-3→1.3e-7 (~4桁減衰/400step)。TGV KE cost 2.7% (case/09 run_0024) とのバランスで **σ=0.05 を既定に採用** | ref (σ 較正) |
| `run_0018_cell_keep_cbd_dissmat005` | 同 IC + **keepDissType=2 (matrix ES, σ=0.05)** | A_cb 1e-3→**7.1e-8** (scalar σ=0.05 の 1.3e-7 より良) — 音響のみ c' 込み固有値の選択的減衰が市松に有効 | ref (Step 2 matrix L1b) |
| `run_0019_cell_keep_cbd_dissmat005_fullc` | matrix σ=0.05 を**フル c** (lowMachPrecond=0, 旧結線) で | A_cb 1e-3→**9.2e-9** (最速)。ただし TGV KE cost 4.35% (case/09 run_0026) と3.2倍悪化 → c' の価値の実測 | ref (c' vs フル c 切り分け) |
| `run_0020_cell_keep_cbd_dissmat005_knobchk` | **keepDissCprime ノブ分離後**の動作確認 (lowMachPrecond=0 + 既定 keepDissCprime=1) | A_cb→7.5e-8 ≈ run_0018 (c' 挙動を lowMachPrecond と独立に再現、差は atomicAdd ノイズ級) | ref (ノブ分離検証) |
| `run_0021_cell_keep_cbd_dissmat_fullc_s0015` | フル c で σ 弱化の代替可否: σ=0.015 | A_cb→**3.9e-8** = c' σ=0.05 (7.1e-8) と同等以上。**単一マッハ領域なら σ 弱化で c' 代替可** (KE も 1.10% = case/09 run_0027) | ref (σ掃引・c'代替検証) |
| `run_0022_cell_keep_cbd_dissmat_fullc_s004` | 同 σ=0.04 | A_cb→3.1e-10。**σ=0.05 (9.2e-9) より速い非単調** = 2Δ 定在音響波の過減衰領域を観測 | ref (σ非単調の記録) |
| `run_0023_node_keep_cbd_pure` | **node L1**: node 変換 (median-dual 35937CV) + 市松 IC (COORD パリティ・17969/17968)、pure KEEP | **node でも null-mode 実証** (step0 全 rms ~1e-10 = 機械ゼロ、A_cb 400step 不変)。explicit RK3 でも発散せず | ref (node L1a) |
| `run_0024_node_keep_cbd_dissmat005` | node + matrix σ=0.05 (既定) | A_cb 1.03e-3→3.5e-5 (**1.5桁**)。**cell (4.2桁) より大幅に遅い**: 内部 2.7e-5 に対し**周期継ぎ目ノードが 3.3e-4 (0.6桁) で律速**、継ぎ目から染み戻る。node periodic の合併 DOF 周りの散逸の効きが弱い (state mirror タイミング等、要別調査) | ref (node L1b・既知の弱点記録) |

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

# 35. uniform periodic box (Taylor-Green)

全 6 面 Cartesian periodic の 3D 構造化 hex box (32³, `mesh/tg.msh`)。一様流/Taylor-Green 初期値で
周期境界・3D 経路の最小検証に使う。3D median-dual (M4, [`discretization-median-dual-3d`](../../plans/active/discretization-median-dual-3d.md))
の双対メッシュ生成 (体積総和・closure) 検証ケースでもある。

## 計算 run 一覧

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_uniform` | cell モード一様流ベースライン (convMethod=0) | `res_*.h5`, `residual_history.csv` | ref |
| `run_0002_node_dualcheck` | node モード (`discretization: node`) で 3D median-dual 生成の幾何検証 | 体積 dual=248.05=primal (relErr 3.5e-8)・closure normalized 2.5e-6・負体積 0 (`convert.log`)。solver run は未実施 | active |

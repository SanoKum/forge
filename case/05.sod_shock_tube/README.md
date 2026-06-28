# 05. sod shock tube

## 3D periodic 化検証 (median-dual M4 §4.5)

2D/1D sod を 3D box ([0,1]×0.1×0.1, 200×8×8 hex) に押し出し、transverse (y,z) を Cartesian periodic にして node periodic を検証。IC=sod (x<0.5 左状態 ro1/P1, 右 ro0.125/P0.1)、inviscid、explicit RK3。

| run_* | 目的・設定 | 結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_node3d_periodic_sod` | node + y,z periodic, RK3, SLAU/1次 | **node explicit が step16-24 で発散**。periodic は動作 (union-find 3417 merge=16281−12864, periodic 除外)。slip/z-only periodic でも同 step16-18 発散 → periodic 非依存の **node-explicit-shock 一般不安定** (TG/backstep と同類)。4-way overlap エッジに float closure 起因の微小異常 (Ux 7e-7) はあるが論理バグでない (free-stream 8-way が機械ゼロ=多重加算無し) | 破棄予定 (node-explicit 限界) |
| `run_sod3d_cell_ref` | **cell** + y,z periodic (参照) | **完走・ro 0.125..1.0 正常・x=0.6 で ro std=0.0 = 完全 spanwise 均一** → periodic setup・メッシュ健全。cell periodic は正しく機能 | ref (periodic setup 検証) |

**結論**: 3D periodic の setup・幾何 gather は正しい (cell 完走+均一、free-stream 機械ゼロ)。node periodic も無罪。sod node explicit は node-explicit-shock 限界で未走 (要 implicit dual-time 等)。

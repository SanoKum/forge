# 05. sod shock tube

## 3D periodic 化検証 (median-dual M4 §4.5)

2D/1D sod を 3D box ([0,1]×0.1×0.1, 200×8×8 hex) に押し出し、transverse (y,z) を Cartesian periodic にして node periodic を transient 衝撃で検証。**物性・IC は run_0002_slau_cpg と同一** (cpg 物理単位: 左 P1e6/T2000, 右 P1e5/T400, roe=P/(γ-1), cp1038.8/γ1.4)。**教訓**: 非次元 `initial:"sod"` (P=1) は forge thermo と非互換 (T≈0.003K で roe 破壊→P 均一) なので使わない。cfl は run_0002 同様 0.2。

| run_* | 目的・設定 | 結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_node3d_periodic_sod` | **node + y,z periodic**, run_0002 同一 config (cpg 物理単位 PL1e6/TL2000・PR1e5/TR400, cfl0.2, SLAU, RK3) | **完走 (200 step NaN なし)・物理 sod 構造 (ro 0.84..1.92, P 1e5..1e6, 衝撃右伝播)・spanwise std=0〜1e-6 完全均一** = node periodic が transient 衝撃で機能。node は cell より前面がクリーン (決定的) | **検証完結** |
| `run_sod3d_cell_ref` | cell + y,z periodic (同 config, 参照) | 完走・物理 sod。滑らか域 spanwise std=2e-7 均一、衝撃前面のみ std=0.5 波打ち (cell atomicAdd ノイズ、periodic とは別) | ref |
| `run_node3d_slip_lowcfl` | **IC バグ調査の記録**: 当初 `initial:"sod"` (非次元 P=1) は forge CPG thermo (cp1038.8,R296.8) と非互換で T≈0.003K→roe 破壊→**P 均一の異常状態**。これに node checkerboard が乗り発散していた | 破棄予定 (IC バグ記録) |

**結論**: node sod 3D periodic が **安定・物理的・spanwise 完全均一**で動作 (検証完結)。当初「node-explicit-shock 限界」とした発散は**誤診**で、真因は壊れた IC (非次元→P均一) + 過大 CFL だった (user 指摘で判明)。periodic 実装は transient 衝撃でも健全。

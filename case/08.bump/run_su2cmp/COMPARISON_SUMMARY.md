# bump 三者比較 (forge node / cell vs SU2) — y=0.25 ライン

- 形状: x∈[0,3], y∈[0,1], 下壁 bump (peak 0.1 @x=1.5)、非粘性・slip 上下壁、低マッハ (Pt/Ps=1.19, M~0.5)。
- サンプリングライン: y=0.25 水平 (下端から25%、全 x で流体内)。
- SU2: 2D Euler ROE/Venkat、mesh/bump2d.su2 (32004 quad)、`-> 収束` (rms[Rho]≈-6.25 plateau)。
- forge-node: run_tmp_hlle_node/res_150000.h5、HLLE、`check_convergence -> PASS (converged)`。
- forge-cell: run_tmp_hlle_cell/res_40000.h5、HLLE、`check_convergence -> NOT CONVERGED (plateau ~6e-4, atomicAdd floor)`。

## y=0.25 ライン SU2 基準の相対差
| 量 | node max/mean | cell max/mean |
|----|---------------|---------------|
| P  | 0.04% / 0.01% | 0.25% / 0.12% |
| T  | 0.01% / 0.00% | 0.07% / 0.03% |
| |U|| 0.06% / 0.03% | 0.49% / 0.30% |

## 結論
- **forge-node は SU2 と <0.1% で一致** (P/T/|U|)。node 系統の物理は SU2 検証ソルバと整合。
- **forge-cell は ~0.25-0.5% ずれる**。これは cell の **atomicAdd 非決定性で残差が ~6e-4 で頭打ち=未収束**の刻印
  (物理誤差ではなく収束度の差)。cell は機械収束しないため SU2 比較では node より劣る。
- 図: bump_line_compare.png

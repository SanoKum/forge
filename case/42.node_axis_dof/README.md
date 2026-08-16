# case/42.node_axis_dof — node 軸対称の軸ノード離散化の検証 (case/41 解析モデル流用)

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
4. 未解決: `nodeValueAtNode: 1` の出口∩軸コーナー 1 ノード (M −0.11, P +4.5%; GG/LSQ 同一)、入口∩軸 −1%。
   NS (粘性) は境界半割面 dcc 退化のため本フラグ未対応 (ghostless 化が前提, plan §5 Step 1)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| `run_0001_ref41_0057_dirichlet` | case/41 `run_0057_axismach_srcfix_n2000` の入力+最終場のコピー (生産: node, `nodeAxisDirichlet: 1`, Euler cplus n2000)。比較の基準 | 軸−偶関数外挿 max 0.0030、‖ΔM‖∞ 直読 0.222% / 外挿 0.147% Md。`res_12000.h5` | ref |
| `run_0002_axisdof_euler` | 0001 と同一入力で `nodeAxisDirichlet: 0` (軸ノードを DOF として解く、値位置=双対重心) | NaN 0・STEADY・plateau (生産と同型)。真空崩壊なし (ρ_axis/ρ_row1 0.989–1.018)。軸 M 欠損 max 0.033 (0.5–0.8%)、‖ΔM‖∞ 直読 0.888% | active (**軸 DOF 対照**) |
| `run_0003_axisdof_van_euler` | 0002 + `nodeValueAtNode: 1` (試作: 値位置=ノード座標, r̄→rEff) | 軸−外挿 ≤0.0069 (出口コーナー除く)、cell 参照との差 0.0031、‖ΔM‖∞ 直読 0.257% / 外挿 0.165%。出口∩軸コーナー 1 ノード M 3.892 | active (**恒久策試作の根拠**) |
| `run_0004_axisdof_van_lsq` | 0003 + `gradLSQ: 2` | 0003 と 3 桁一致 (コーナー含む) = 勾配法は無関係 | active |
| `run_0005_cell_ref` | 同 nozzle.msh を cell 変換 (品質 PASS AR 9.9/skew 0.41)、`interp_field` で 0057 場から段階起動 (soft 1次 cfl0.5 3000 → 本段) | NaN 0・STEADY・plateau 0.8 桁 (cell の既知床)。セル中心列の偶関数外挿 ‖ΔM‖∞ 0.151% Md — node 外挿の独立参照 | active (**参照**) |
| `run_0006_regress_dirichlet` | 0001 と同一入力 (生産 `nodeAxisDirichlet: 1`) を**パッチ後バイナリ**で再実行 (既定経路のビット不変確認) | 最終場は 0001 と最大相対差 4e-6 (ro/roUx/roe/P; node の run 間非決定性レベル) = 既定経路は不変 | active (回帰) |
| `optest/h{0.02,0.01,0.005}_*` | 離散連続式演算子テスト (`optest/run_optest.py`, 直管 1×0.2 m 一様 quad, 1 step 陽解法, `limiter: 0` は線形場再構成のため意図的) | `optest/optest_results.json`。要約は上記 1. | active (**演算子テスト**) |

解析スクリプト: `analyze_axis_row.py` (軸行 vs 第一内点), `compare_axis_runs.py` (→ `axis_runs_compare.png`)。

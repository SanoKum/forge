# レビュー依頼プロンプト: node (median-dual) 軸対称離散化の 3 方式と今後の方針

以下を Codex への最初のプロンプトとして使う。

---

forge (自作 圧縮性 FVM ソルバ, CUDA/float32, cell 中心と node 中心 median-dual の 2 離散化、今後は node 主体) の**軸対称解析の離散化方針について、忖度なしの厳しいレビュー**をしてほしい。実装の正しさ・数学的健全性・運用判断のすべてを疑ってかかり、主張は必ずコードと実測で検証すること。あなたの結論が「現方針は誤り」であっても構わない。

## 背景と経緯 (2026-08-03〜05 に確定した事実)

対象: 軸対称ベルノズル (case/40, Pt 4 MPa / Tt 1500 K / ε9, SLAU+SST, 構造化 quad 221×65 を node=median-dual で解く)。

1. **基底欠陥の発見**: forge の軸対称は「r 重み幾何」(B 流儀: 体積・面積に r を乗じる。`methods/axisymmetric/implementation.md`)。node では軸上ノードが半 CV の中心になり、その体積・面積が r̄∝Δr で消えるため、**軸行が radial 圧力平衡からデカップルして真空まで過膨張** (ρ が隣接行の 1/10〜1/50、有効圧力負 → EOS 床 P=1 Pa ピン)。laminar でも発生 (=SST 以前の欠陥)。崩壊軸行と第一内点行の偽せん断が SST の k 生産シート (周囲の ~10³ 倍) を作り、**陰解法 (block-DPLUR, 点 Jacobi 型 defect-correction) では遅発性発散** (roK が e-fold ~2000 step で成長し step ~8000–10600 で爆発)。cell 離散化は同メッシュで健全。経緯全記録: `plans/active/boundary-node-nozzle-wall-outlet-stability.md` §2.6。

2. **対策を 3 方式実装済み** (いずれも case/40 で 12000 step 検証済み、全域 2 次 + block-DPLUR 陰解法 cfl_pseudo 4):

| 方式 | 中身 | 代表 run | rms_ro | row1 k 帯 | η_CF |
| --- | --- | --- | --- | --- | --- |
| `nodeAxisDirichlet: 1` (現 runner 既定) | 軸ノードを解かず、radial 代表隣接ノードから毎ステージ状態コピー (∂q/∂r=0 の 1 次, u_r=0)。全残差 0 化 + DPLUR 5 行単位行化 | `case/40.../run_0022_node_runner_nobfo` (旗艦 ALL PASS は run_0016) | **1e-9 床 (PASS 可)** | 3.4 | 0.9896 |
| `axisymMethod: 1` (SU2 流) | planar 幾何 + 1/y ソース (非粘性 4 式+解析 Jacobian、粘性 AuxVar μv/y 系の GG 勾配込み、SST 1/y 移流拡散)。軸ノードは通常 DOF として解く。SU2 `.external/su2-src` から移植 | `case/40.../run_0026_node_su2axis` | **implicit ~1e-4 プラトー** (explicit RK3 なら 3e-7) | **1.35 (最清浄)** | 0.9904 |
| `axisRFloor: 3.0e-4` (r 床) | r<床 の面・セルは r 重みを床値で打切り (幾何が消えない=環状流路近似)。床帯は hoop ソース・Jacobian・u_r/r 不課。**要点: ソース面積を解析 A_planar でなく床適用後の離散閉性 Σ_f S_f (A_closure_x/y) に置換** (素朴な面床は混在 CV の圧力閉性が壊れ step 19 で発散した) | `case/40.../run_0027_node_rfloor` (+24000 step soak で床帯縁の k 帯 ~230 が有界平衡と確認) | **~1e-7 プラトー** | 230 (有界) | 0.9906 |

3. **SU2 式 implicit プラトーの切り分け済み** (`plans/active/axisymmetric-su2-source-formulation.md` 変更ログ): 5×5 solve の double 化=不変 / 非粘性ソース Jacobian OFF=平均流微改善だが SST 悪化 / 粘性ソース OFF=不変 / SST ソース OFF=不変。ソース CFL ~0.1 で形式的に stiff でない。explicit は 3e-7 到達 = 空間スキーム健全。→ 結論「対角 Jacobian の工夫では解けない。FVS 近似 Jacobian と実残差の不整合が大 Δτ で露呈するもの」。次の打ち手候補として **radial line-implicit (原義 DPLUR の線緩和、median-dual では line detection が前準備)** を挙げている。

4. **別系統の未解決欠陥**: 壁ノード列の**等圧エントロピー市松** (P 滑らか 0.6% / T·ρ 逆相 ~12%、ベル部 |ΔΔT|~206 K、**定常解自体が市松を持つ** — IC を平滑化してもビット同一の市松に再収束、laminar で増悪 426 K、cell は 32 K)。壁ノードは u=0 で移流がなく断熱壁は T のアンカーを持たない。切り分け: 伝導は活きている (κ×10→振幅×0.42)、粘性勾配補正項は無罪、μ(T) は一部寄与 (×0.65)。項別計装が未実施。`plans/active/boundary-node-nozzle-wall-outlet-stability.md` §2.8。

5. **その他の既知事項**: node−cell の η_CF 差 +1.1〜1.2% は未帰属 (katoLaunder でも軸修正でもない、node 2 次化の離散差。Rao 照合は未実施)。`isAxisymmetric`/`axisymMethod` は mesh ブロックへ移設済 (physProp は deprecated 読み)。壁市松対策で warm start 時の壁列 T 平滑化が runner に入っている (定常解は不変なので実質 IC 清浄化のみ)。

## レビューしてほしい点 (優先順)

1. **SU2 移植の検証** (`solver_density_cuda/cuda_forge/axisymmetricSource_d.cu` の `axisymmetricSourceSU2_d` 系 / `ransSource_d.cu` の axisym 節 / `timeIntegration_d.cu` の `isAxisymmetric==2` Jacobian) を `.external/su2-src` (CSourceAxisymmetric_Flow, ResidualDiffusion, ResidualAxisymmetricConvectionDiffusion, BC_Sym_Plane) と突き合わせ、符号・項・Jacobian 写像の誤りを探せ。意図的スキップ 2 点 (エネルギー粘性ソースの `+ρk` turb_ke 項 = 次元不整合の疑い / エネルギー Jacobian の `1/2` 整数除算 = SU2 側バグの疑いで数学的に正しい形を実装) の判断も是非を問う。
2. **axisRFloor の数学的健全性**: 環状近似 (軸に半径 r_floor の slip 円柱を置くのと等価) の保存性・整合次数・床値のメッシュ依存性 (軸行重心 < 床 < 第一内点行重心 という条件が崩れるメッシュでの挙動)。離散閉性面積 A_closure の考え方は正しいか。床帯縁の有界 k 帯 (~230) の物理的意味とリスク。
3. **nodeAxisDirichlet の精度リスク**: 軸行を 1 次コピーにすることの系統誤差 (軸 Mach・軸上衝撃・推力積分への影響)。「一番深く収束するから既定」という運用判断は妥当か。
4. **implicit プラトー診断の妥当性と次の一手**: 切り分け (T1–T4) から「対角では解けない」と結論したのは正しいか。radial line-implicit は妥当か、それとも近軸限定 Δτ 制限・CFL ランプ・Newton-Krylov・単純に「implicit で場を作り explicit で磨く」二段運用のほうが費用対効果が良いか。
5. **生産既定の選択**: 現在は nodeAxisDirichlet (最深収束) を既定にしている。設計最適化ループ (η_CF 評価, 12000 step ≈ 20 秒) の用途で、この選択は正しいか。η の 3 方式差は 0.1% 内だが、これを「方式によらず η は信頼できる」と読んでよいか。
6. **壁エントロピー市松**: §2.8 の診断と「項別計装が次ステップ」という方針の妥当性。もっと直接的な対処 (壁ノードのエネルギー扱いの変更等) はあるか。
7. **全体**: 我々が見落としている構造的リスク (float32、median-dual 境界弱形式、SLAU との相性等) を挙げよ。

## 検証に使えるもの

- 文書: `plans/active/axisymmetric-su2-source-formulation.md` (SU2 式と r 床の全記録) / `plans/active/boundary-node-nozzle-wall-outlet-stability.md` (§2.6–2.8) / `methods/axisymmetric/implementation.md` / `case/40.nozzle_design_tool/README.md` (run 一覧と 3 方式比較表)
- コード: `solver_density_cuda/cuda_forge/axisymmetricSource_d.cu`, `timeIntegration_d.cu`, `ransSource_d.cu`, `variables.cpp` (r 重み・床・閉性面積), `mesh/mesh.cpp` (axis_flag/axis_rep), `main.cpp` (assembleResidual の呼び順)
- SU2 原典: `.external/su2-src/SU2_CFD/src/numerics/flow/flow_sources.cpp`, `include/numerics/turbulent/turb_sources.hpp`, `include/solvers/CFVMFlowSolverBase.inl`
- 実測: `case/40.nozzle_design_tool/run_0015〜0027` (各 `CONVERGENCE_VERDICT.txt`, `residual_history.csv`, `res_12000.h5`)。診断 env ゲート `FORGE_DIAG_SU2JAC_OFF` / `FORGE_DIAG_SU2VISC_OFF` / `FORGE_DIAG_SU2SST_OFF` で再実験可 (実行は `solver_density_cuda/tools/run_case.sh <run_dir 絶対パス>`)

## 出力形式

- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (ファイル:行 または run の数値) と対案をセットで。
- 最後に「軸対称解析の推奨構成 (短期・中期)」を 1 つに絞って提言すること。両論併記で逃げないこと。

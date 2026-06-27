# node-centered SST 壁関数の近壁カバレッジ不足による過剰乱流 (調査)

- **日付**: 2026-06-28
- **対象**: `case/18.backstep` (2D 平面 backstep, RANS-SST, `wallTreatmentSST=1` automatic)
- **きっかけ**: user 指摘「node の乱流分布が明らかにおかしい・wall distance が怪しい」。node の vis_turb が段差凹コーナーで vis_t/l=6789 (cell=77) と異常。
- **関連**: [discretization-node-wall-implicit-dirichlet.md](../../plans/accepted/discretization-node-wall-implicit-dirichlet.md), [methods/turbulence/theory.md](../../methods/turbulence/theory.md) §6.5, memory [[node-mode-periodic-and-backstep-status]]

## 結論 (データ確定)

**node-centered SST の壁関数近壁処理 (wf_pk 生産置換 + ω ピン) が「壁ノード1点」しか覆わず、cell が覆う「第一内層」を取りこぼす**のが過剰乱流の真因。生産項 Pk が支配 (ω でも decouple でも factor でも strain でもない)。

### k 収支実測 (段差垂直壁の角, x≈2.1)

診断フィールド `wf_pk`(≥0=壁関数処理/-1=未処理)・`Pk_diag`(確定生産)・`src_jac_k`(=β*ω, Dk=src_jac_k·ρk) を出力に追加して計測 (run_0060 node / run_0061 cell, 収束場から 200 step)。

| 位置 | wf_pk | Pk(生産) | Dk(消散) | k | ω | vt/l |
| --- | --- | --- | --- | --- | --- | --- |
| NODE y=0 (壁ノード) | 0.0 (処理) | 0.0 | 104 | 51.9 | 18.8 | 3293 |
| NODE y=0.062 (第一内層) | **-1 (未処理)** | **166** | 99 | 52.0 | 17.8 | 2221 |
| NODE y=0.125 | -1 | **190** | 80 | 52.2 | 14.2 | 1957 |
| CELL y=0.031 (第一セル) | 0.03 (処理) | 0.0 | 25 | 3.1 | 75.9 | 49 |
| CELL y=0.094 | 0.05 (処理) | 0.1 | 9 | 3.8 | 20.9 | 73 |
| CELL y=0.281 | 0.64 (処理) | 0.6 | 15 | 4.6 | 31.1 | 77 |

- **cell**: 段差垂直壁の第一セル列 (x≈2.07, y=0.03/0.09/0.28) が**全て wf_pk≥0**。淀みコーナーでは wf_pk≈0 (u_τ≈0) なので生産がほぼ消え、k=3-5 に抑制。
- **node**: **壁ノード (y=0) のみ wf_pk≥0**。第一内層以降 (y=0.06/0.125/0.25) は **wf_pk=-1=未処理**で標準 SST 生産 Pk=μt·S²=**166〜228** が働く。μt が高い(=k/ω)ため生産が暴走し k=52 に到達 → vis_t/l=3293。

### 棄却した仮説 (記録)

1. **「omega 行 decouple 欠落」** → 誤り。`applySSTPointImplicit_d` に omega decouple は実装済みで機能 ([update_d.cu](../../solver_density_cuda/cuda_forge/update_d.cu) `decouple_wall_omega`)。壁ノード ω はピン値 (≈6ν/β₁y_eff²) を保持している。
2. **「壁 ω 値が SU2 より低い (factor 6 vs 60)」** → 誤り。automatic wall は係数 **6** が正 (解析漸近, [theory.md](../../methods/turbulence/theory.md) §6.5(b))。60 は低Re 壁解像 BC の過大指定で automatic には高すぎる (user 指摘で確認)。SU2 `BC_HeatFlux_Wall`(係数60) は今回の MARKER_WALL_FUNCTIONS run と別経路で、比較対象違いだった。
3. **「凹コーナーで strain S 偽増幅」** → 誤り。node の S はむしろ cell より低い (node 2.8 vs cell 19.4)。
4. 真因は上記の通り **wf_pk 生産置換の近壁カバレッジ不足** (および同根で ω ピンも壁ノードのみ)。

## 【解決・2026-06-28】修正完了

修正は plan [`turbulence-node-wall-function-coverage.md`](../../plans/accepted/turbulence-node-wall-function-coverage.md) に集約。要点:

- **正しい設計**: **生産 P_k 置換 (`wf_pk`) のみ壁ノード+第一内層ノードの両方に被覆。ω ピン/残差ゼロ化/decouple は壁ノードのまま据え置く。**
- **失敗した当初案** (wf_pk と ω ピンの両方を第一内層 irep へ移設) は: ① 凹コーナーで底壁/段差面が同じ irep を共有し ω ピン値が race、② 壁ノードの ω アンカが外れ ω 崩壊 → k がさらに暴走 (52→246-350)。**ω ピンは壁ノードに残すのが鍵**。
- **結果** (`case/18.backstep` run_0065/0067): 段差コーナー vis_t/l 6789→77・k 52→3 (cell 同等)、場平均 424→207 (cell 198)、x_R 6.71→7.63 (cell 7.95/SU2 7.89 整合)。cell ビット不変 (run_0066)。
- **残課題**: node 壁ノードは cell 第一セルより壁から遠く ω ピンが ~1/4 低いため、再付着せん断層で μt ピーク残 (vt/l~5800 vs cell 990、局所 ~27 node)。場平均・x_R は整合。完全一致には node near-wall ω 解像 or k 上限 (SU2 SetTurbVars_WF の k Dirichlet) の追加対策が要る。

## 修正方針 (要 plan) — 当初の検討 (記録)

cell では「壁隣接セル=第一内層」が壁関数処理を受ける。node でこれに整合させるには、**壁関数処理 (wf_pk 生産置換 + ω ピン) を壁ノードだけでなく第一内層ノードにも適用**する必要がある。具体策の候補:
- `computeWallFrictionSST_d` / `rans_wall_scalar_boundary_d` が選ぶ代表内部点 (SU2 Normal_Neighbor, 既存) を**生産置換と ω ピンの適用先**にも使う (現状 u_τ 計算だけに使い、wf_pk/ω は壁ノード ic に書いている)。
- すなわち wf_pk と ω ピンを「壁ノード ic」でなく「第一内層ノード irep」に書く (cell の wall-adjacent cell に対応)。壁ノード自体は no-slip Dirichlet (速度0) + ω は ic にも保持。
- decouple (point-implicit dω=0) も適用先ノードに合わせる。

cell 側の「セルをピンする是非」(user 別途議論希望) と設計が絡むため、plan 化時に整合させる。

## 成果物 / run

- 診断 run: `case/18.backstep/run_0060_node_kbudget` (node) / `run_0061_cell_kbudget` (cell)。収束場 (run_0059/run_0057) から 200 step、`wf_pk`/`Pk_diag`/`src_jac_k` 出力。
- コード: `variables.hpp` に `Pk_diag` 追加 + 診断 4 フィールドを output_cellValNames へ、`ransSource_d.cu` に `Pk_diag[ic]=Pk` 書き込み (一時診断, 本調査用)。

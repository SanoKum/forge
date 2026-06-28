# 作業プロンプト: backstep 細密メッシュで node+KEEP+WALE LES を回す

> 使い捨ての作業プロンプト (notes/sessions)。終わったら破棄/更新してよい。
> ブランチ: `feature/median-dual-3d`。ビルド: `solver_density_cuda/build` で `ninja forge`、
> 実行は必ず `solver_density_cuda/tools/run_case.sh <run_dir>`。
> `export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial:$LD_LIBRARY_PATH`。

## ゴール
3D backstep の **node-centered LES** を**細密メッシュ**で回し、**解像された乱流**(剥離せん断層の roll-up・
再付着・spanwise 3D 乱れ)を得る。旧メッシュ(spanwise たった 8 セル)は WALE が過剰散逸して定常化していた。

## なぜ新メッシュか (前セッションの結論)
- 旧 `backstep_3d.h5` は 228×48×**8** (Δz=0.5H)。LES には粗すぎ、WALE の μ_sgs∝Δ²=(体積)^(2/3) が
  **ν_lam の 40–380 倍** (RANS 級) になり解像乱流を全部潰して**定常化**していた。
- 同じ旧メッシュでも **SLAU ILES (LESorRANS=0, WALE 無し)** は活発な乱流が出ていた
  (late 連続スナップショットで Ux \|Δ\|~38 m/s)。→ 犯人は WALE の過大渦粘性。
- 対策 = **spanwise を細かく** (Δz を小さく→μ_sgs を LES 相応 O(1–10)×ν に下げる)。

## 新メッシュ (作成済み)
- `.geo`: [`case/18.backstep/mesh/backstep_les.geo`](../../case/18.backstep/mesh/backstep_les.geo)
  — **230(x) × 50(y) × 80(span z)**。段差 H=1, x=2、x[0,30] y[0,3] z[0,4] (span 4H, **Δz=0.05H**)。
- `.msh`: `case/18.backstep/mesh/backstep_les.msh` (gmsh -3 で生成済、**933,606 nodes / 989,002 elements**, 98MB)。
- bcond physID: 1=inlet, 2=outlet, 3=top wall, 4=bot wall, 5/6=periodic(side2/side1, Cartesian dz=∓4)。
- 再生成: `cd case/18.backstep/mesh && gmsh -3 backstep_les.geo -o backstep_les.msh -format msh4`。

### HDF5 変換
新しい run_dir に `.msh` と `solverConfig.yaml`(node モード)・`bcondConfig.yaml` を置き、その中で:
```
cd <run_dir>
../../../solver_density_cuda/build/convertGmshToForge backstep_les.msh backstep_les.h5
```
- `meshFileName`/`valueFileName` を `backstep_les.h5` に合わせる。出力 ~303MB。
- **変換に ~50 秒かかる** (wall_dist が brute-force O(N×N_wall)=3.6e10; kdtree 未リンク・converter が -fopenmp 未適用)。
  一回きりなので許容。速くしたいなら nanoflann 導入 or converter を OpenMP ビルド (別タスク, 前セッションで分析済)。
- 変換後 **`tools/check_mesh_quality.py backstep_les.h5`** で AR≤1000/skew≤0.9 を確認 (構造化 hex なので通る見込み、Δz0.05/Δx0.13/Δy0.06 で AR~2.6)。

## ソルバ設定 (node + KEEP + WALE LES)
`solverConfig.yaml` の要点 (旧 `run_node3d_keep_pure_wale_long/solverConfig.yaml` を雛形に):
```yaml
mesh: { discretization: "node", nodeWallDirichlet: 1, meshFileName/valueFileName: backstep_les.h5 }
solver: "KEEP"               # KEEP 中心流束 (+Roe 散逸 = keepDissipation:1)
physProp: { thermalMethod:0(CPG), visc:0.001, cp:1004.5, gamma:1.4 }  # Re_H~40000
time:
  unsteady:1, dualTime:0
  deltaT: { control:1(cfl), cfl: 1.0, dt_max: 0.001 }   # ← CFL 上限は 1.0 (下記)
  timeIntegration:3 (explicit RK3), nStepOuter: 大きく (例 150000)
  outStepInterval: 1000〜2000
space: { convMethod:2(MUSCL), limiter:2(venkata), keepDissipation:1 }
turbulence: { LESorRANS:1, LESmodel:1(WALE) }
initial: "backstep"
bcondConfig: 旧 backstep の inlet_uniformVelocity(Ux≈49 等)/outlet_statPress/wall/periodic を流用。
```

### 重要な既定事実 (前セッションで確定)
1. **CFL 上限 = 1.0** (explicit RK3)。sweep 実測: cfl 0.5/1.0 安定、**1.5 で発散** (P 非物理 1–227000 Pa)、2.0/2.5 は step~10 で NaN。→ **cfl 1.0 を使う**。細密化で許容 cfl が下がる可能性があるので、発散したら 0.5 へ。
2. **periodic state-mirror は実装済** (commit 70af723, `periodicMirrorNSState`)。周期 master/slave の保存量を初期化時+各 RK stage で同期する。これが無いと seam に圧力欠陥 (~150 Pa) が出るが、**もう自動で効くので不要な心配**。ただし **seed の摂動は z 周期的に** (z=0=z=4) 与えること (mirror が拾うが物理的整合のため)。
3. **pure KEEP も可**: `keepDissipation:0` (Roe 散逸無し) も state-mirror 後は spanwise checkerboard が出ず綺麗だった (mirror 前の市松は seam desync が主因)。低散逸 LES を狙うなら pure KEEP + WALE も選択肢。まずは安全な `keepDissipation:1` (KEEP+Roe) 推奨。
4. **k/omega は LES では vestigial**: LESorRANS=1 では roK/roOmega 輸送は全 gate で no-op (凍結)。出力の k=roK/ro が ro 変動で動いて見えるだけ。無害。気になるなら `setInitial` で roK/roOmega=0 にしてもよい (WALE は k 不使用)。
5. **dt は global 一様** (unsteady=1 && dualTime=0 → `setDTlocal_uniform_cell_d`)。local time stepping ではない (時間精度あり)。

## LES 品質チェック (必須)
- **解像乱流が出ているか**: late 連続スナップショット (例 res_(N-1000) vs res_N) で **Ux \|Δ\| が大きい**(=非定常)ことを確認。
  小さい (例 \|Δ\|<0.2 mean) なら laminarize=失敗 (WALE 過大 or メッシュまだ粗い)。SLAU ILES は \|Δ\| max~38 だった。
- **vis_turb/ν_lam を確認**: LES なら O(1–10)。まだ O(100) なら WALE 過大 → さらに細密化 or **ILES に切替**
  (`LESorRANS:0` + `keepDissipation:1`(KEEP+Roe 散逸を implicit SGS に); 旧メッシュでも ILES は乱流が出た)。
- 派生量 (再付着長 x_R, せん断層, 時間平均) を報告するときは **`tools/check_quasisteady.py` / `check_convergence.py` の VERDICT を必ず引用** (AGENTS 必須)。過渡ピークを定常値と誤認しない。非定常 LES は「収束」でなく統計的定常を時間平均で見る。

## 実行・運用 (AGENTS 準拠)
- 必ず複製した新 `run_*` で実行、`run_case.sh` 経由 (VERDICT/residual_history.png 自動生成)。
- 序盤で NaN/Inf 確認 (`residual_history.csv` の rms 列)。case README (`case/18.backstep/README.md`) の run 一覧表を同期。
- `case/` の run 成果物 (res_*.h5 等) は commit しない。ソース/docs/plan のみ commit、feature ブランチへ push。

## 参考 (実装の所在)
- KEEP: `cuda_forge/convection/convectiveFlux_keep_d.inc.cuh` (commit 79b4e67)。
- node WALE: `cuda_forge/turbulent_viscosity_d.cu` の `dimGrid_normalcell` 起動。
- periodic state-mirror: `cuda_forge/periodicNode_d.cu` `periodicMirrorNSState_d_wrapper` (commit 70af723)、`main.cpp` の init + RK stage 呼び出し。
- 設計: [`plans/active/discretization-median-dual-3d.md`](../../plans/active/discretization-median-dual-3d.md) §4.5.9、
  [`plans/active/convection-keep-revive-node.md`](../../plans/active/convection-keep-revive-node.md)。
- 旧 run (参照): `case/18.backstep/run_node3d_iles_unsteady` (SLAU ILES, 乱流出た)、
  `run_node3d_keep_pure_wale_long` (旧粗メッシュ WALE, 定常化した反例)。

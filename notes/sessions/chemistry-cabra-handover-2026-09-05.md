# 引継ぎ: 有限速度化学 (H₂ 燃焼) と Cabra 検証 — 2026-09-05 時点

Claude CLI で続きを行うための引継ぎ文書。まず [`AGENTS.md`](../../AGENTS.md) のルール (run ディレクトリ命名・README の run 表・
`check_convergence.py` の VERDICT 必須・commit 運用) に従うこと。

## 1. 作業場所 (最重要)

| 項目 | 値 |
| --- | --- |
| worktree | `/home/sano/work/forge-chem` (**メインの `/home/sano/work/forge` ではない**。メインは別作業中のため触らない) |
| branch | `feature/chemistry-finite-rate` (origin に push 済み, HEAD `42568a09`) |
| ビルド | `solver_density_cuda/build` (→ `.build-native/relwithdebinfo` の symlink), CUDA arch **86 必須**。`cd solver_density_cuda/build && cmake --build . -j 8` |
| 実行 | `bash solver_density_cuda/tools/run_case.sh <run_dir>` (別バイナリは `tools/run_case_bin.sh <bin> <run_dir>`)。`build/forge` を直接呼ぶと hook に止められる |
| Python | Cantera 入り: `/home/sano/work/forge/.venv-chem/bin/python` (setup/compare スクリプト用)。CEA2: `/home/sano/work/forge/.venv-cea/` |
| 計算結果 | すべて `forge-chem/case/4[5-8].*/run_*` (commit 対象外)。tmp には置いていない |

`solverConfig.hpp` を変えたら**必ずフルビルド**相当になっているか確認 (差分ビルドの CUDA obj 取りこぼしで step0 NaN 凍結の前例あり)。

## 2. 何ができているか (機能)

計画: [`plans/active/chemistry-finite-rate-h2.md`](../../plans/active/chemistry-finite-rate-h2.md) (§9 変更ログが時系列の正本)、
仕様: [`methods/chemistry.md`](../../methods/chemistry.md)、設定: [`procedures/solver-settings.md`](../../procedures/solver-settings.md)
(`physProp.chemistry`, `lowMachPrecond`, `speciesPrecondDt` の節)。

- Arrhenius + 三体 + 逆反応 (Kc) + Lindemann/Troe、解析 Jacobian、種ブロック point-implicit (`jacobianMode 2`)、
  反応熱の陰的注入 (案C 予測子)、Jacobian 凍結 `jacobianInterval` (**推奨 5**)、Strang 分離 (非定常 RK 用)、PaSR (`tci 1`, `tciCmix`, `tciTauChem 1`)。
- 機構: `solver_density_cuda/tools/mechanisms/` (Jachimowski 13sp33r / 9sp20r, GRI H2/O2 NASA-9)。参照解: `tools/chem_reference_cantera.py`。
- **定常陰解法の反応流は `speciesImplicitCoupling: 2` + `jacobianMode: 2` 必須**、`cfl_pseudo` 上限 ≈ 4 (ノズル) / 1 (Cabra)。

### 検証済み
- 0-D 着火 (case/35) Cantera 一致、Q1D ノズル再結合 (case/46 run_0008) Y_OH 一致。
- Burrows–Kurkov (case/47, best `run_0025_v2_delta10_kx3`): 出口組成一致、全温ピーク 1.08 vs 1.18、着火位置 5–9 cm vs 実験 18–25 cm (未解決)。

## 3. 2026-09-05 に直したバグ 3 件 (Cabra 立ち上げで発覚)

すべて commit `876ec8b4`。詳細は [plan lowmach 変更ログ 2026-09-05](../../plans/accepted/time_integration-lowmach-preconditioning.md)、
[plan outlet §2.12](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)。

1. **node × `lowMachPrecond>=2` の境界ノード凍結** — `implicit_defect_correction_block_precond_d` (timeIntegration_d.cu) に node 用行処理
   (境界半割面の粘性対角スキップ・軸/壁行 decouple) を追加。これ以前の node precond run (case/39, 43) は境界が凍結していた可能性。
2. **TP 亜音速 `outlet_statPress` の γ 混用** (boundaryCond_d.cu) — セル局所 γ_mix で特性構成を統一。低マッハ出口の一様逆流 (−65 m/s) の指紋。
   BK は壁近傍の亜音速出口列が変わり残差 0.5 % 変化 (再現ノイズ 1e-5)。BK の出口プロファイル再確認は **`run_0028_exit_recheck` で完了 (2026-09-05, 結論維持)** — case/47 README の run 表を参照。
3. **多成分 × 定常擬似時間 × precond の前線速度不整合** — `time.deltaT.speciesPrecondDt` (**既定 1** = 化学種は前処理拡大前 Δτ を使う。0 は旧挙動で A/B 用)。
   setDT が `dt_local_sp` を用意し化学種 DPLUR が読む。k/ω は同構造の潜在問題 (未対応)。

## 4. Cabra (case/48) の現状

[`case/48.cabra_h2n2/README.md`](../../case/48.cabra_h2n2/README.md) の run 表が正本。結果ページ (図クリックで拡大):
https://claude.ai/code/artifact/27dc951c-fd98-4185-a64e-96afe231bfbd

- レシピ: `setup_cabra_case.py run_dir --chem 0 --nstep N --out M --conv 1 --relax 0.5 --cfl 1.0 --iccol 1 [--restart res.h5]`
  (node 軸対称, precond 2, ε 0.15, 2 次, SST 壁関数, TP 9 種)。**cfl_pseudo 2 は step 932 で NaN** (`run_0066`)。
- 混合場 (化学 OFF): `run_0065` → `run_0067` 計 26000 step、有界、NOT CONVERGED (発達途上)。最新 `run_0067/res_20000.h5`。
- 反応 ON (`--chem 1 --jac 2 --ji 5 [--tci 1 --cmix C]`): `run_0068`→`run_0069` (tci 0), `run_0070` (PaSR C_mix 1), `run_0071` (PaSR C_mix 4)。
  **全て x/d≈20 で自着火 → 火炎基部が上流伝播 → リップ付着** (実験 H/d≈10)。PaSR の C_mix は推移を変えない。
  下流 (z/d 26) の半径 T は実験と整合、中心軸 T は z/d≤15 一致・z/d 20 で早すぎる立ち上がり。全 run NOT CONVERGED (過渡)。
- 切り分け用オプション (`--single/--jetcof/--jetn2/--iccol 0/--planar/--precond/--coupling/--far`) は README 参照。`run_0001`–`run_0064` は
  バグ切り分けのエビデンスで破棄予定。

### ユーザ判断待ち (次の一手)
1. 自着火火炎向け TCI (EDC / 輸送 PDF / 着火遅れ時間ベース) を実装して浮き上がり位置まで合わせに行く、か
2. 燃焼加熱器の検証は下流の温度・組成 + BK 出口組成で足りるとして加熱器の設計計算へ進む、か。

どちらでも、入口 k/ω (ジェット乱流強度) と SST 混合速度の感度は未実施 (中心軸 H₂ 減衰がやや速い)。

## 5. 未了・注意

- ~~BK 出口プロファイルの再確認 (上記 3-2)~~ → `run_0028_exit_recheck` で完了 (結論維持)。
- BK の早着火 (5–9 cm vs 18–25 cm) は未解決 (run_0028 でも Y_OH>1e-4 位置が 4.36→3.93 cm と緩慢な上流ドリフト継続)。
- `chemistry_source_d` の占有率最適化 (3 KB スタック/126 レジスタ) は未着手。
- メイン `feature/median-dual-3d` に化学の初期 3 commit が混入している (ユーザ了承済み・そのまま)。force-push 禁止。
- `pkill -f <run名>` は自分のシェルを殺す (パターンが自分のコマンド行に一致) — 使わない。

## 6. 参照

- 論文 (untracked): `papers/combustion/` (Jachimowski, Burrows–Kurkov, Cabra NASA/CR-2004-212887)。
- 調査メモ: [`notes/investigations/chemistry-finite-rate-h2-survey.md`](../investigations/chemistry-finite-rate-h2-survey.md)。
- メモリ (Claude 側): `chemistry-finite-rate-direction`, `node-precond-boundary-freeze-and-tp-outlet-gamma`, `node-inlet-wall-corner-conflict`。

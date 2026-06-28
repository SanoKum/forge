# backstep の P・ρ 振動の正体と打ち手 — 低マッハ odd-even チェッカーボードと `lowMachPrecond=2`

> 調査メモ (2026-06-28)。`case/18.backstep` で node/cell を問わず観測される P・ρ の振動の正体を
> データで切り分け、有効な打ち手を A/B で実証した記録。発端プロンプト:
> [`notes/sessions/backstep-pressure-oscillation-session-prompt.md`](../sessions/backstep-pressure-oscillation-session-prompt.md)。
> 関連 plan: [`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md)、
> [`convection-slau2-lowmach.md`](../../plans/accepted/convection-slau2-lowmach.md)。

## 結論 (要約)

1. **振動の正体 = 低マッハ再循環域に局在する P・ρ の odd-even チェッカーボード** (collocated 圧力-速度デカップリング)。
   段差リップ (x=2.0, y=1.0) を1セル過ぎた縦断面 (x≈2.13) の y 分布で、**P が節点ごとに ±50 Pa 鋸歯振動**
   (2次差分 rms = 全振幅の 48–54%)。**node・cell 両方**に出る。**速度 (Ux) は滑らか**。
2. **正体でないもの**: 時間リミットサイクルの大振幅ではない (プローブ時系列 <10 ppm)・速度の市松ではない・
   リミッタチャタリングではない (`limiter_P` は安定)・**float32 桁落ちではない** (振幅 5e-4 ≫ float32 ε≈1e-7、
   かつ float32 のままスキーム散逸変更で消えるため精度は無関係)。
3. **打ち手 = `lowMachPrecond=2`** (Weiss–Smith 完全前処理 block-DPLUR, plan Phase 4)。
   **P 市松振幅を ~76–77% 削減** (cell 78.4→18.4 Pa, node 96→22 Pa) し、かつ **残差プラトーを打破**:
   同一 restart からの cell A/B で control(lmp0) は rms_roUx≈**0.29 で停滞** (roUy 上昇) なのに対し lmp2 は
   **5.5e-4 まで収束** (~500× 低、全列 falling)、rms_roOmega も control 1.13→lmp2 **8.9e-4**。**プラトーの正体が
   低マッハ市松リミットサイクル**だったと確定。quasisteady は pmax/machmax とも STEADY。
4. **`lowMachPrecond=1` (流束のみ前処理) は不可**: cfl_pseudo=1 で発散 (NaN)、0.3 で高残差停滞。
   RHS のみ c' 化し LHS Jacobian が c_hat のままで不整合。**完全前処理の =2 のみが整合収束**。

## 背景

`case/18.backstep` は低マッハ (M≈0.11、Pt=100000/Ps≈99118、U≈38 m/s) 後退段 RANS-SST を block-DPLUR
陰解法で回す。直近まで残差プラトー (`NOT CONVERGED (stalled/plateau — needs scheme change)`) が続き、
P・ρ の「振動」が node/cell 共通で観測されていた。代表 run: `run_0067` (node)・`run_0057` (cell)・`run_0058` (SU2)。
メッシュ `mesh/backstep_planar.h5` (10997 node / 10720 cell 平面2D)。段差は **x=2.0** (入口 x∈[0,2] が y∈[1,3]、
下流 x∈[2,30] が y∈[0,3]、段差高 H=1)。

## 切り分けの手順とデータ

### (1) 空間 odd-even — グローバル指標では薄まる

`res_*.h5` 最終場で近傍平均偏差 `dev = φ - mean_nb(φ)` の市松相関 `ρ = <dev·mean_nb(dev)>/<dev²>` を計算
(滑らか→ρ>0, 市松→ρ<0)。全域では `P` ρ≈−0.03〜−0.08 と弱く、市松は**局在しているため薄まる**。

### (2) 角部1セル下流の y 分布 — 市松が顕在化 (user 指摘で判明)

段差角を1セル過ぎた x≈2.13 (cell は x≈2.20) の縦列で y 方向に並べ、2次差分 (odd-even) rms / 全振幅:

| 場 | P odd-even [Pa] | P rel | ρ odd-even | ρ rel |
| --- | --- | --- | --- | --- |
| node baseline (`run_0067` 100k) | 96.3 | 0.538 | 1.83e-3 | 0.266 |
| cell baseline (`run_0057` 100k) | 78.4 | 0.476 | 1.41e-3 | 0.245 |

リップ下 (y<1) の再循環域で P が節点ごと ±50 Pa の鋸歯。Ux は滑らか (rel≈0.11)。図:
[`run_0067.../corner_oddeven_baseline.png`](../../case/18.backstep/run_0067_node_wfcov_long/corner_oddeven_baseline.png)。

### (3) 時間性 — リミットサイクルは微小

収束後 (step≥50000) の固定点プローブの snapshot 間変動: P std 0.5–7 ppm、ρ ~1 ppm。`limiter_P` は 0.97–0.99 で安定
(チャタリングなし)。→ **大振幅の時間振動でもリミッタ起因でもない**。残差プラトーは微小振幅リミットサイクル。

### (4) 精度 — float32 桁落ちは棄却 (倍精度ビルド 2×2 で決定的に確認)

市松振幅 ±50 Pa / 1e5 ≈ **5e-4 (相対)** は float32 ε≈1.2e-7 の数千倍 → round-off ではなく**スキームが減衰
させない実在の数値モード**。**`build-double` (全 float64) を現行 HEAD で再ビルドし、precond × precision の 2×2 を
同一 restart・10k step で実施** (cell):

| case | build | precond | P odd-even [Pa] | rms_roUx | rms_roOmega |
| --- | --- | --- | --- | --- | --- |
| B0 (`run_0074`) | float | lmp0 | 69.9 | 0.291 | 1.13 |
| B1 (`run_0080`) | **double** | lmp0 | **69.9** | 0.291 | 1.13 |
| B2 (`run_0076`) | float | lmp2 | 18.5 | 0.086 | 0.20 |
| B3 (`run_0081`) | **double** | lmp2 | **18.5** | 0.086 | 0.20 |

**B1=B0・B3=B2 が完全一致** → **倍精度は市松にも残差にも全く効かない。改善は 100% `lowMachPrecond=2` の前処理物理**。
`lowMachPrecond=2` 経路は内部で Jacobian 組立・1/β 分母を部分倍精度化 (timeIntegration_d.cu:153/1055) するが、
**全 double 化 (B3) で float lmp2 (B2) と一致**するため内部 double も勝因でない。`mixed-precision-axisym-refuted` と
同じ「精度に見えてスキーム」のパターンを実データで確認。
注: `build-double` は double のレジスタ消費増で既定 `blocksize=512` が起動上限超過 → `FORGE_CUDA_BLOCKSIZE=128` で実行。

#### (4b) c' 前処理「物理」 vs double Jacobian「構築」の決定的切り分け

`lowMachPrecond=2` は (a) RHS 流束を c' 散逸化 + (b) LHS 前処理 Jacobian を **double 構築**
(`build_jacobian_split`/`eos_split_jacobian_general_closed`、`lowMachPrecond=0/1` は float の `accumulate_split_jacobian_cf`)。
「市松を直すのは (a) 物理か (b) double 構築か」を **2 つの既存ノブで分離** (cell, 同一 restart, 10k):

| run | c' 前処理 (RHS) | double Jacobian (LHS) | P odd-even [Pa] |
| --- | :---: | :---: | --- |
| `run_0074` lmp0 float | ✗ | ✗ | 69.9 |
| `run_0082` lmp0 + `implicitSolvePrecision=1` | ✗ | **✓** | **69.9** (= baseline・bit 一致、残差も 0074 と同値) |
| `run_0075` lmp1 (cfl0.3) | **✓** | ✗ (float) | **21.1** |
| `run_0076` lmp2 | ✓ | ✓ | 18.5 |

→ **double Jacobian 単独 (c' なし) は市松に全く効かない** (0082=0074)。**c' 前処理は float Jacobian のまま 78→21 Pa と
大半を落とす** (0075, `lowMachPrecond=1` は float の `accumulate_split_jacobian_cf` を使う)。lmp2 の double Jacobian は
21→18.5 の僅かな上積みで、本質的役割は**前処理系 (Γ⁻¹A の ~1/β 悪条件) を安定に solve すること** (lmp1 は停滞、lmp2 は収束)。
**結論: 市松を直すのは c' 前処理の物理。double は前処理系の数値条件付け (収束安定化) のためで、市松除去の主因ではない。**

### 機構

SLAU の質量流束の圧力散逸項 `χ/c_hat·ΔP` が音速 `c_hat` でスケールし、M→0 で圧力-速度カップリングが O(M) に
縮退 → odd-even を減衰できない。これは plan
[`time_integration-lowmach-preconditioning.md`](../../plans/accepted/time_integration-lowmach-preconditioning.md) §9 が
ノズル/チャンバーで特定済みの機構で、backstep 再循環の市松も同一原因と確認された。

## 打ち手の A/B (cell、`run_0057` の収束場から restart、cfl_pseudo=1)

| run | 設定 | 結果 (10k step) |
| --- | --- | --- |
| `run_0074` | control (lmp0) 10k | rms_roUx **0.29 で停滞** (roUy 上昇)・roOmega 1.13。P odd-even ~70 Pa 維持。A/B 基準 |
| `run_0073` | `lowMachPrecond=1` cfl1 | **DIVERGED (NaN)**。step0 rms_ro 5e-3 に跳ね倍々成長 |
| `run_0075` | `lowMachPrecond=1` cfl0.3 | NaN回避だが高残差停滞 (rms_roUx 0.41)。P odd-even 21 Pa |
| `run_0076`/`run_0077`/`run_0079` | **`lowMachPrecond=2` cfl1** (`run_0079`=完全収束 200k 相当) | **NaN無・全残差 falling**。P odd-even **18.4 Pa (−76%)**、rms_roUx **0.29→5.5e-4**・roOmega **1.13→8.9e-4**、quasisteady STEADY |

`run_0077` (cell lmp2 50k) の角部 P 分布は baseline の ±50 Pa 鋸歯が ±10 Pa 級の小波へ平滑化 (ρ も同様、Ux 不変)。
図: [`run_0077.../corner_oddeven_lmp2_vs_baseline.png`](../../case/18.backstep/run_0077_cell_lmp2_long/corner_oddeven_lmp2_vs_baseline.png)。
場は物理的 (vt/l 平均 196 vs baseline 198、P 範囲同等、NaN 無)。市松指標は step10k=30k=50k で同値 (18.4) = 形状定常。

### node 側 (`run_0078`, `run_0067` 収束場から restart)

`lowMachPrecond=1` は **DIVERGED** (cfl1)。`lowMachPrecond=2` は **発散せず** (旧 `run_0042` の =2 step0 NaN は explicit
時代の交絡)、P odd-even **96→22 Pa (−77%)**、平均流残差 −2.5〜2.8桁。ただし **rms_roOmega は ~9.5 で依然プラトー**
= node の omega プラトーは低マッハ市松とは**別問題** (README 既知の壁近傍 ω ソース挙動)。

## 推奨 / 残課題

- **低マッハ剥離・再循環ケースは `lowMachPrecond=2` + `blockDPLUR=1` を既定候補**とする (cfl_pseudo=1 で整合収束)。
  `lowMachPrecond=1` は LHS 非前処理で不整合のため使わない (発散/停滞)。`precondEps=0.15`。
  注: `implicitSolvePrecision=1` (double solve) は `lowMachPrecond=2` と非対応 (config ガードで停止)。
- **残課題**: ① 残留 ~18 Pa の小市松 (リップ近傍は物理せん断層構造を含む) のさらなる低減、
  ② node の rms_roOmega プラトー (別 SST 問題、本市松とは独立)、③ `run_0079` で cell lmp2 を完全収束まで回し最終 VERDICT 確定。

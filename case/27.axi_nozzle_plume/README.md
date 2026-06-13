# 27.axi_nozzle_plume — Rao Mach-4 ノズル + プルーム (TP gas / CEA 照合)

`case/23.axi_nozzle/run_0097_slau_sst_dilAB_2nd_rebuilt` の設定を母体に、ノズル+遠方場プルームを
一から計算する独立ケース。SLAU・軸対称・Rao ノズル (設計 $M_e=4$)、入口 $P_t=4\,\mathrm{MPa},\ T_t=1500\,\mathrm{K}$、
背圧/周囲 40 kPa。高 $T_t$ なので thermally-perfect gas の $c_p(T)$ 効果が大きく、CEA 照合に好適。

## メッシュ・初期条件

- メッシュ: `mesh/axi_nozzle_2d_publicrao_plume.geo` (+ `.msh`)。単一流体領域 (regionId=8)、
  境界 physID = inlet(1)/wall(3)/axis(4)/outer1(5)/outer2(6)/outer3(7)。27499 セル。
- IC: `convertGmshToForge` が `solverConfig.yaml` の `initial: "uniform_p101325_u10"` を適用して
  一様 ($P=101325,\ u=10$) のベース IC を書く。これがノズル内部 IC の実体 (スクリプト相当)。
  プルームは単一領域なので領域パッチ (`setInitial_plume_outer.py`, 旧 2 領域メッシュ用) は不要。
- 起動: 一様 IC は $P_t=4\,\mathrm{MPa}$ に対し桁違いに低いため、**1 次精度で起動** → 2 次へ移行
  (確立済みワークフロー)。

## 手法決定 (パラメータスタディ)

run_0001 (1 次・非粘性で確立) から warm-start し 2 次精度で比較 (全て NaN/発散チェック済):

| run | 手法 | 結果 |
| --- | --- | --- |
| `run_0002_exp_2nd_cfl03` | explicit RK3, cfl 0.3 | 安定 |
| `run_0003_exp_2nd_cfl06` | explicit RK3, cfl 0.6 | 安定・最速 |
| `run_0006_exp_2nd_cfl08` | explicit RK3, cfl 0.8 | **NaN (step 333)** |
| `run_0007_exp_2nd_cfl10` | explicit RK3, cfl 1.0 | **NaN (step 17)** |
| `run_0004_imp_2nd_cflp5` | implicit blockDPLUR, cfl_pseudo 5 | **NaN (step 2170)** — MUSCL+陰解法不安定 |
| `run_0005_exp_2nd_precond` | explicit + lowMachPrecond=2 | **NaN (step 2)** — 低Mach前処理は超音速核で不適 |

→ **採用手法: explicit RK3 / 2 次 (convMethod=1) / cfl 0.6**。preconditioning は超音速ノズルでは
発散 (M≪1 用)、陰解法は MUSCL と高 cfl_pseudo で発散。残差は $\sim2\times10^{-3}$ で頭打ちだが
これはプルーム (超音速ジェットの衝撃ダイヤモンド) が本質的に非定常なため。**ノズル核は定常**。

## TP 効果 / CEA 照合

採用手法で N2 を計算 (高温で $c_p(T)$ が顕著):
- `run_0008_cpgN2_2nd` — CPG-N2 (定数 $c_p=1038.8,\ \gamma=1.4$)
- `run_0009_tpN2_2nd`  — TP-N2 (`thermalMethod=2, species:[N2]`, NASA-9)

`compare_cea.py` がノズル内部 ($x\le0.0675$) 中心線の $T(P)$ を、$P_t,T_t$ からの
**NASA-N2 $c_p(T)$ 1 次元等エントロピー解 (= 非反応フローズン CEA)** と比較。

**結果** (`cea_comparison.png`):
- **TP-N2 は NASA $c_p(T)$ 等エントロピーに 0.04% で一致** (CEA を再現)。
- **CPG-N2 は 5.99% ずれる** (定数 $c_p$ で $c_p(T)$ 効果を取りこぼす)。
- 膨張を通して TP-N2 は CPG-N2 より高温 (N2 の $c_p$ が高温で大きく、運動エネルギーへの変換が
  少ない実在気体効果)。Mach はほぼ一致。

## 粘性 / RANS + TP の検証 (M3 輸送係数・RANS×TP)

非粘性 CEA 照合に加え、TP-N2 を段階的に粘性化して未検証だった経路を確認 (全て NaN チェック済):

- `run_0010_tpN2_visc_kintheory` — TP-N2 + `viscMethod=2` (kinetic theory 混合粘性・熱伝導, laminar)。
  安定。**`vis_lam` = 1.28e-5〜5.51e-5 Pa·s が T[194,1551]K に対し物理的な N2 粘性** (M3 を実粘性 TP 流れで検証)。
- `run_0011_tpN2_rans_sst` — 上記 + RANS SST (dilatation 補正 2)。安定 (残差 1.08e-3、乱流粘性が
  プルームせん断層を平滑化)。`vis_turb` [0,0.57]、`k` [0,9.5e4] が物理的 (RANS×TP のエンドツーエンド検証)。

> M4 の乱流化学種拡散 ($\mu_t/\mathrm{Sc}_t$) は単成分 N2 では dormant のため未検証 (多成分 RANS が必要。plan §10)。

## 再現手順

```bash
# 1) run_0001: メッシュ変換 + 1次起動
cd run_0001_cpg_1st_establish && convertGmshToForge axi_nozzle.msh axi_nozzle.h5 && forge
# 2) 採用手法で N2 (run_0001/res_5000.h5 から warm-start; TP は roe を NASA-N2 基準へ書換)
#    -> run_0008_cpgN2_2nd, run_0009_tpN2_2nd
# 3) CEA 照合
python3 compare_cea.py run_0008_cpgN2_2nd/res_8000.h5 run_0009_tpN2_2nd/res_8000.h5
```

## 計算 run 一覧 (抜粋)

| run_* | 目的・設定 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_0004_imp_2nd_cflp5` | 軸対称 SLAU 2nd block-DPLUR 陰解 (cfl0.5) | **float では step3000 付近で NaN** (近軸不安定) | active |
| `run_regr_cf` | 回帰: 閉形式 FVS 既定 float。run_0004 同条件 | float は同様に NaN (≈step4000) = 閉形式由来でない | active |
| `run_regr_cf_double` | `implicitSolvePrecision=1` (double solve)。run_0004 同条件 | **完走・NaN なし → 近軸 double solve が発散を安定化**。plan [precision-mixed-axisym.md](../../.github/plans/precision-mixed-axisym.md) | active |

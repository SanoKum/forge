# case/46.nozzle_h2o2_kinetics — 燃焼室 2788 K の化学非平衡ノズル (有限速度化学 Phase 2 検証)

有限速度化学 ([`plans/active/chemistry-finite-rate-h2.md`](../../plans/active/chemistry-finite-rate-h2.md),
[`methods/chemistry.md`](../../methods/chemistry.md)) の**定常陰解法**検証。形状は case/44 va3 ノズル
(run_0091 のメッシュ `nozzle.h5`, node 軸対称, 23725 ノード, $r_t$=0.2057 m, 出口 $A/A^*$≈15.1)。
入口は CEA 平衡 (`cea_A3000_eq_chamber.json`: 反応物 3000 K 指定 → 解離後 $T_c$=2788 K, $p_c$=11.39 bar, $X_{OH}$ 1.37 %, $X_{NO}$ 3.1 %)。
機構は Jachimowski 1988 13 種 33 反応 (`mech13.yaml`) + 不活性 AR, CO2 (CO は CO2 に合算, N₂O は落とす) の 15 種。

- `setup_kinetics_case.py run_dir [--chem 0|1] [--jac 1|2] [--coupling 0|1|2] [--cfl] [--nstep]`: IC (run_0091 収束場を $T_c/1161$ でスケール、組成は入口組成で一様)・BC・config を生成。
- `q1d_kinetics_ref.py nozzle.h5 out.csv [--frozen]`: 準 1 次元参照解 (Cantera `mech15.yaml`、スロート CEA 平衡状態から A/A*=1.02 で開始)。
  `q1d_finite_rate.csv` / `q1d_frozen.csv`。
- `compare_axis.py run_dir[,run_dir...] step out.png`: 軸線の T / M / Y_OH / Y_NO を Q1D と重ねる。

Q1D と forge 軸線の T/M は 1D 面平均 vs 軸値の差 (出口 M で 3.7–3.85 vs 3.9–3.95) があるので、**化学の検証は $Y_{OH}$ の減衰と
frozen との差 (ΔT, ΔM)** で見る (`--start-from-forge` で forge の M=1.3 軸点から始めた Q1D が同一組成始点の比較)。
参照値 (Q1D, 出口 x=4.986 m): finite-rate T=939 K, M=3.846, $Y_{OH}$=6.9e-5; frozen (スロート組成凍結) T=899 K, M=3.890, $Y_{OH}$=5.1e-3。
CEA: 平衡 974 K / M 3.80、スロート凍結 891 K / M 3.89。

## 計算 run 一覧

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_frozen` | 化学 OFF (frozen, 15 種輸送のみ), implicit cfl_pseudo 2, 12000 step | 完走。残差プラトー (`rms_roe` 1.0 絶対, run_0091 系と同じ warm 床) = `check_convergence` NOT CONVERGED (STALLED)。frozen 基準 | ref |
| `run_0002_kinetics_jac2` | 化学 ON, jacobianMode 2, coupling 2 (反応熱は陽的注入) | **step 5 で NaN 発散**: スロート付近の再結合熱 $\dot Q$~7e11 W/m³ が 1 擬似 step で内部エネルギーの 30 % | 破棄予定 (発散記録) |
| `run_0003_kinetics_jac1` | 同 jacobianMode 1 | step 3 で NaN | 破棄予定 |
| `run_0004`–`run_0007_dbg_*` | 切り分け: coupling 1/0, cfl_pseudo 0.2 | 全て NaN (cfl 0.2 で step 60 まで延命) → 陰解結合方式でなく陽的反応熱が原因と確定 | 破棄予定 |
| `run_0008_kinetics_jac2_heatimpl` | **反応熱の陰的注入** (案C 予測子の線形化, `chemistry_heat_inject_d`) + jac 2 + coupling 2, cfl_pseudo 2 (frozen と同じ), 12000 step | **安定・完走**: `rms_roe` 1.77e6→86.5 (4.3 桁), `rms_ro` 3.8e-7, `check_convergence` は frozen 基準と同じプラトー (rms_roUy STALLED) = NOT CONVERGED だが出口量は 4000→12000 step で T 896.8→897.3 K, $Y_{OH}$ 6.365e-5→6.385e-5 (準定常)。**軸線 $Y_{OH}$ は forge スロート状態始点の Q1D と全域一致** (出口 6.39e-5 vs 6.44e-5)、NO は凍結 (3.07e-2, Q1D 同値)。frozen 比 出口 T +58 K・M −0.045 (Q1D: 平衡スロート始点で +40 K・−0.044)。図 `axis_vs_q1d_from_forge.png`, `axis_vs_q1d_cea_throat.png` | ref (Phase 2a 検証エビデンス) |
| `run_0009_kinetics_cfl4` | run_0008 と同設定で `cfl_pseudo` 4, 4000 step | 安定・NaN なし。step 4000 の残差は cfl 2 と同水準 (`rms_ro` 4.2e-7, `rms_roe` 171 vs 104) → 収束加速はしないが上限内 | ref (CFL 上限) |
| `run_0010_kinetics_cfl8` | 同 `cfl_pseudo` 8 | NaN にはならないが `rms_roUx` 29.5・`rms_roe` 2.9e4 と 2 桁悪化 (リミットサイクル) → **反応 ON の実用上限は cfl_pseudo ≈ 4** (frozen 系の既知の単桁上限と同等) | ref (CFL 上限) |
| `run_0011_regr_pasr_off`, `run_0012_regr_pasr_off_rerun` | PaSR 実装後の回帰 (既定 `tci: 0`, run_0008 と同設定 300 step) と同一設定の再実行 | step 299 の残差は run_0008 と 5 桁目で差、再実行同士も同程度の差 (`rms_roe` 2429.7 / 2430.6 / 2426.2) = 化学種移流 atomicAdd の run 間ノイズ内で不変 | 破棄予定 (回帰記録) |

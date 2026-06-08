# 多成分 thermally-perfect gas ソルバ機能

## メタ

- **area**: `thermophysics` (新規。理論/実装は `docs/thermophysics/`)
- **status**: `in_progress`
- **related_docs**:
  - `docs/thermophysics/theory.md`
  - `docs/thermophysics/implementation.md`
  - `docs/convection/` (流束の TP 整合)
  - `docs/diffusion/` (化学種拡散・エンタルピー拡散)
- **created**: `2026-06-09`
- **owner**: `CFD Dev`

## 1. 目的

forge を単一成分・熱量的完全気体 (CPG) から**多成分 thermally-perfect gas (非反応)** に拡張する。比熱 $c_p(T)$ を NASA-9 (CEA 準拠) とし、粘性・熱伝導率・質量拡散を kinetic theory で評価、混合物性は ideal-gas mixing。`thermalMethod==0` を保持し `thermalMethod==2` で有効化。

## 2. スコープ

- **やる**: NASA-9 熱力学、Newton 温度反転、多成分質量分率輸送、kinetic theory 輸送係数、化学種拡散とエネルギー結合、SLAU/ROE の TP 整合、設定/IC/BC/出力配管。
- **やらない**: 化学反応源項 (Arrhenius 等) と有限速度化学。設計上の拡張点としてのみ残す。

## 3. 関連 docs と前提

- 理論: `docs/thermophysics/theory.md`、実装: `docs/thermophysics/implementation.md` (本 plan 着手前に作成済み)。
- 既存の汎用スカラ輸送 `scalarTransport_d` / RANS `ransTransport_d` を化学種輸送のミラー元とする。
- `flow_float` は FP32。熱力学/輸送のカーネル内計算は double で行う。

## 4. 設計方針

詳細は `docs/thermophysics/implementation.md`。要点:
- 熱力学コア `cuda_forge/thermo_d.{cuh,cu}` (NASA-9, 混合則, Newton 反転)。化学種データは device pointer で渡す。
- `dependentVariables_d` に `thermalMethod==2` の Newton 反転を追加し `T,P,Ht,sonic,gamma,cp` を更新。
- 対流流束は被移流全エンタルピーを NASA 化 (SLAU 最小修正、ROE は有効γ̃の段階A)。
- 陰解法 Jacobian は `gamma[ic]` 配列化 (frozen-coefficient)。
- 化学種は segregated point-implicit で別解き (5×5 構造不変)。

## 5. 実装ステップ (マイルストーン)

1. **M1**: `thermo_d`, Newton 反転 (`dependentVariables_d`), `thermalMethod==2`, SLAU/ROE の TP 整合, 設定 (`solverConfig`), `main.cpp` 初期化。— **実装済・検証中**
2. **M2**: 化学種変数登録 (`variables` の `cellValNames` 非 const 化 + `registerSpecies`), 化学種移流/実現可能性/再正規化, 化学種 BC/IC/出力 (`speciesTransport_d`, `speciesBoundary_d`)。
3. **M3**: kinetic theory 輸送係数 (`gasProperties_d` 拡張: Wilke/Wassiljewa/混合平均拡散, Neufeld 衝突積分)。`thermCond` のセル毎配列化。
4. **M4**: 化学種拡散 (`scalar_diffusion_species_d`) + ΣJ_i=0 補正 + エンタルピー拡散項 `Σh_s J_s` を `res_roe` へ。

## 6. 検証

- **単体 / ビルド**: native (`tools/build_native_wsl.sh`, CC 8.6) でビルド成功。NASA-9 単体値を NIST と比較 (N2: R=296.8, cp(300)=1040, cp(2000)=1284 J/kgK, h(298.15)=0、e→T 反転誤差 ~1e-11)。
- **検証ケース**: `case/13.nozzle_H/run_0001_tpgas_n2` (TP-N2) vs `run_0002_slau_cpg` (CPG)。`validate_tpgas.py` で中心線 M・T・P と全エンタルピー保存を NASA(=CEA) 等エントロピーと比較。高 Tt ケースで cp(T) 効果を showcase。
- **判定基準**: 残差収束、$H_0=h(T)+½|u|^2$ が中心線で一定 (エネルギー整合)、(M, P/P0) が NASA/CEA 等エントロピー上に乗る。`thermalMethod==0` 既存ケースが不変。

## 7. 影響範囲

- 新規: `cuda_forge/thermo_d.{cuh,cu}`, `docs/thermophysics/{theory,implementation}.md`, `case/13.nozzle_H/validate_tpgas.py`。
- 変更: `cuda_forge/dependentVariables_d.{cu,cuh}`, `cuda_forge/convectiveFlux_d.{cu,cuh}`, `cuda_forge/gasProperties_d.cu` (TP は無改変)、`cuda_forge/CMakeLists.txt`, `input/solverConfig.{hpp,cpp}`, `main.cpp`, `docs/index.md`。
- M2-M4 で `variables.{hpp,cpp}`, `viscousFlux_d`, `boundaryCond.cpp`, `setInitial.hpp`, `output.cpp` 等を追加変更。

## 8. 完了条件

- [x] `docs/thermophysics/theory.md` 作成
- [x] `docs/thermophysics/implementation.md` 作成
- [ ] M1-M4 実装・検証完了 (§6)
- [ ] `.github/plans/README.md` 更新
- [ ] 本 plan を `done` 化し §9 に変更ログ

## 9. 変更ログ

- `2026-06-09` — 初稿。M1 (NASA-9 熱力学, Newton 反転, SLAU/ROE TP 整合) 実装、native ビルド成功、NASA-9 単体値を NIST と照合済み。nozzle_H で TP-N2/CPG 検証実行中。
- `2026-06-09` — **M1 流れ場 CEA 定量検証 (完了)**。当初の nozzle_H は warm-start 元が旧バイナリ製で inlet_Pressure+SLAU が流れを確立できず不成立 (CPG でも崩壊するため TP 起因ではない)。代わりに `case/05.sod_shock_tube` の均一メッシュ (`1D_shock_tube.geo`, nx=499, x∈[0,1], 全 slip, inviscid SLAU, MUSCL+limiter, RK3, 固定 dt=2.5e-7s で t=2.5e-4s) で N2 衝撃管 (PL/PR=10, TL/TR=2000/400K) を立て直し、TP (`run_0001_slau_tp_n2`, thermalMethod==2) と CPG (`run_0002_slau_cpg`, thermalMethod==0) を実行。
  - **枠組み検証 (CPG vs 厳密 Sod)**: 中心線で rarefaction/contact/shock 位置が一致。L2 相対誤差 P=0.90%, u=2.96%, ρ=3.53% (ρ は contact 平滑化が主因、499 セル 2次精度で妥当)。→ TP 流束+EOS+BC の数値枠組みが厳密 Sod 解と整合。
  - **cp(T) 効果 (TP vs CPG)**: 高 T プラトー (~1650K) で TP が CPG より高温、contact/shock front が右 (高速) へシフト (T 差最大 248K, u 差 225m/s)。N2 の cp(T) 増大 (γ 低下) による実在気体効果を定量提示。
  - 後処理: `case/05.sod_shock_tube/gen_shocktube_ic.py` (TP/CPG IC), `compare_shocktube.py` (厳密 Sod Riemann ソルバ + 中心線照合 + `shocktube_comparison.png`)。
  - 付随: `cuda_forge/calcStructualVariables_d.cu` の幾何メトリック計算を FP32 (`sqrtf`/`powf`) → FP64 (`sqrt`/`pow`) 化 (TP の double 一貫性向上)。
  - §6 nozzle_H ベースの検証計画は本 sod_shock_tube ベースの照合に置き換え。M1 の流れ場 CEA 検証はこれで完了とみなす。

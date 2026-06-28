# EOS 正値化フロアの config 化 (pMin / roMin / tMin)

## メタ

- **area**: `thermophysics`
- **status**: `done`
- **related_docs**:
  - [`methods/thermophysics.md`](../../methods/thermophysics.md) §実装 2「従属変数と温度反転」内「EOS 正値化フロア (config 化)」
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 1. 目的

`dependentVariables_d` の密度・圧力・温度フロアがハードコード (`P_min=1.0` Pa, `ro_min=1e-4`, `T_min=1e-4`) であり、無次元・低圧ケース (Taylor-Green は $P_0=1/\gamma\approx0.714$ Pa) では**初期場全体をフロアがクランプして解を破壊**する。フロア値を `solverConfig.yaml` で指定可能にし、ケースの圧力/密度スケールに合わせて下げられるようにする。

## 2. スコープ

- **やる**: `physProp.{pMin,roMin,tMin}` を追加 (任意キー)。`dependentVariables_d` の全 EOS 経路 (単相 CPG / 二相 CPG / TP gas) のハードコードフロアを config 値に置換。既定値は従来値と同一でビット不変。
- **やらない**: フロア発火時の `roe` 整合再構成ロジック自体の変更 (フロアが発火しない通常運転では恒等)。CPU フォールバック (`dependentVariables.cpp`) は元々フロア無しのため変更不要。

## 3. 関連 docs と前提

- 現象の発見元: `case/09.Taylor-Green` で node/cell pure KEEP の保存性検証中、res_0 で $P\equiv1.0$ (本来 $[0.65,0.78]$) を観測。質量・全エネルギーが ~3e-3 ドリフト、エントロピーが nan 化。
- 既定フロアは大気スケール (`P_min=1.0` Pa = 大気圧の 1/10⁵) 想定。

## 4. 設計方針

- `solverConfig` に `flow_float pMin=1.0, roMin=1e-4, tMin=1e-4` を追加 ([`input/solverConfig.hpp`](../../solver_density_cuda/input/solverConfig.hpp))。`getOptionalValidatedValue` で `physProp` から読む ([`input/solverConfig.cpp`](../../solver_density_cuda/input/solverConfig.cpp))。
- `dependentVariables_d` カーネル引数に `pMin,roMin,tMin` を追加し、置換:
  - `ro_temp = max(ro, roMin)` (旧 `1e-4f`)
  - 単相 CPG: `T_temp = max(intE/(cp/γ), tMin)`、`P_temp = max((γ-1)(roe-ρ·ek), pMin)` (旧 `1e-4f`/`1.0f`)
  - 二相 CPG: `Pn = max(Pn, pMin)`、Tguess の下限 `tMin`
  - TP gas: `Pnew = max(Pnew, pMin)`
- wrapper で `cfg.pMin/roMin/tMin` を渡す。

## 5. 実装ステップ

1. `solverConfig.{hpp,cpp}` にキー追加。 ✅
2. `dependentVariables_d.cu` のカーネル signature・使用箇所・wrapper を更新。 ✅
3. `methods/thermophysics.md` に節を追加。 ✅
4. Taylor-Green で `pMin=1e-6` を指定し、初期場が $P\in[0.65,0.78]$ に復元されることを確認。 ✅

## 変更ログ

- `2026-06-28` — 実装・検証完了。`case/09.Taylor-Green/run_0004,0005,0006` で `pMin=1e-6` を指定し、res_0 の圧力場が解析初期値 ($P\in[0.654,0.774]$, mean 0.7155) に復元されることを確認。既定値据え置きで既存ケースはビット不変。node 双対メッシュ非粘性 pure KEEP で全エネルギー保存が ~1e-4 (フロア破壊時の ~1e-3 から改善)、KE ±1.6% / エントロピー ±0.6% 保存を確認 (pure KEEP の保存性実証)。

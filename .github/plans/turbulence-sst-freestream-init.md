# SST 初期乱流の freestream 初期化 (kInit / omegaInit)

## メタ

- **area**: `turbulence`
- **status**: `done`
- **related_docs**:
  - `docs/turbulence/implementation.md`
- **related_plans**:
  - `.github/plans/turbulence-node-sst-wallfunction.md`
  - `.github/plans/diffusion-node-wall-viscous-distance.md`
- **created**: `2026-06-26`
- **owner**: `CFD Dev`

## 1. 目的

IC に `roK`/`roOmega` が無いとき forge は両者を **0 初期化**していた。`ω=0` は `μ_t=k/ω` が
ill-posed で SST が cold start から発散する (特に node `wallTreatmentSST=1`)。config で
freestream 初期乱流 (`kInit`/`omegaInit`) を与えて 0 初期化を避けられるようにする。

## 2. スコープ

- **やる**: `solverConfig.turbulence.kInit`/`omegaInit` (既定 0) を追加。
  `variables::readValueHDF5` で IC に roK/roOmega が無いとき `roK=ρ·kInit`, `roOmega=ρ·omegaInit`。
- **やらない**: 既存 SST 壁処理 (wt=0/1) の挙動・SST 生産項・wall function 本体は不変。
  cell の SST 過小 mut・SU2 自由流設定は別問題。

## 3. 背景 (バグ実証)

NACA2412 (20.naca_ml, y+~50) の node SST `wt=1` が cold IC (naca.h5, roK/roOmega 無し) から
step~1 で ω 発散。発散源は翼面壁でなく **inlet∩top コーナー** で、step1 に bulk ω≈0 / inlet
ω≈14 (BC 1000 のはず) / コーナー ω≈1e8。原因は IC の roK/roOmega 欠落 → 0 初期化
([variables.cpp] readValueHDF5)。当初「node 壁関数バグ」と疑ったが、真因は SST cold-init。

## 4. 設計

- `solverConfig.hpp`: `flow_float kInit=0.0; flow_float omegaInit=0.0;`
- `solverConfig.cpp` turbulence ブロック: `getOptionalValidatedValue<flow_float>("kInit"/"omegaInit", 0.0)`, `>=0` 検証。
- `variables.hpp`/`variables.cpp`: `readValueHDF5(fname, msh, kInit=0, omegaInit=0)`、欠落時 `ro*kInit`/`ro*omegaInit`。
- `main.cpp`: `cfg.kInit, cfg.omegaInit` を渡す。
- 既定 0 で `ρ·0=0` ＝ 従来動作・ビット不変。

## 5. 検証

case: `case/20.naca_ml/001.test`、メッシュ y+~50 (壁関数領域)。

| run | 初期化 | 結果 |
| --- | --- | --- |
| node SST wt=1, cold IC (旧) | roK=roOmega=0 | step~1 で ω 発散 |
| node SST wt=1, cold IC + kInit=1/omegaInit=1000 (新) | ρ·1 / ρ·1000 | **SURVIVED**, bulk k≈1.25 |
| `run_cmp_node_sst_wt1_fsinit` (120k) | freestream | **CL 0.782, CD 0.0074, mut/mul max508/mean123** = 付着乱流の妥当値 (~0.78-0.82) |

- 比較: node wt=0 (CL 0.429, mut 2199 過剰) / cell wt=1 (0.794, mut 42 過小) / SU2 (0.602, 自由流過大) の中で
  **node wt=1+fix が最も物理的**。
- 既定 (kInit=omegaInit=0) は `ro*0=0` で従来とビット一致 (回帰なし)。

## 変更ログ

- 2026-06-26: 実装・検証完了。`kInit`/`omegaInit` 追加、IC 欠落時 freestream 初期化。
  node SST wt=1 が y+~50 メッシュで cold IC から安定し CL 0.78 (妥当) を出すことを確認。
  これにより node が適切な壁関数 (wt=1) を使えるようになり、「node SST が wt=0 を強制され
  乱流暴走」問題が解消。

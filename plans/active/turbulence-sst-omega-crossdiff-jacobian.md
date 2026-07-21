# SST ω 交差拡散の point-implicit Jacobian 追加

## メタ

- **area**: `turbulence / time_integration`
- **status**: `in_progress`
- **related_docs**:
  - `methods/turbulence/implementation.md` (SST ソース項・point-implicit)
- **related_plans**:
  - [`time_integration-explicit-pointimplicit-sst.md`](../accepted/time_integration-explicit-pointimplicit-sst.md) — 既存 point-implicit の親
  - [`turbulence-iddes-sst.md`](turbulence-iddes-sst.md) — 動機 (case/39 dual-time で ω サブ反復が無収束)
- **created**: `2026-07-21`
- **owner**: `CFD Dev`

## 1. 目的

case/39 周期丘 DDES dual-time で、サブ反復内の残差低下が運動量 ~19x に対し **ω は 0.9x (無収束)**。
ω 点陰対角に唯一欠けているソース Jacobian = 交差拡散項を追加し、ω のサブ反復収束を立てる。

## 2. スコープ

- **やる**: `src_jac_omega += max(CD_ω, 0)/(ρω)` (CD_ω>0 のときのみ = 対角を強める安全側のみ)。
  config `sstCrossDiffJac` (既定 0 = ビット不変)。ransSource_d.cu の 1 カーネル + config 配管。
- **やらない**: 生産項対角化 (過去に試験済・効果なしの記録がコード内にあり)、off-diagonal (k-ω) 連成、
  block 5x5 への k/ω 統合。

## 3. 設計

- SST ω 方程式ソース: $S_\omega = P_\omega - \beta\rho\omega^2 + CD_\omega$,
  $CD_\omega = 2(1-F_1)\rho\sigma_{w2}\nabla k\cdot\nabla\omega/\omega$。
- $\partial CD_\omega/\partial(\rho\omega) = -CD_\omega/(\rho\omega)$。$CD_\omega>0$ (通常の境界層外縁〜wake) で
  負 → point-implicit の対角 $D_\omega$ へ正の寄与 $+CD_\omega/(\rho\omega)$ [1/s]。$CD_\omega<0$ は対角を
  弱めるので除外 (SU2 と同じ)。
- 収束解は不変 (Jacobian のみの変更)。

## 4. 検証

- case/39 500 step 診断 (run_diag_*): サブ反復内 rms_roOmega drop の A/B (off/on)、
  あわせて cfl_pseudo 1→2 の複合も測る。NaN 無し + drop 改善で採用、production へ反映。

## 変更ログ

- 2026-07-21: 起票・実装着手。

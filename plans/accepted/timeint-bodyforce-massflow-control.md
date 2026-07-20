# 体積力の質量流量一定制御 (bodyForceCtrl) — 周期丘駆動

## メタ

- **area**: `time_integration`
- **status**: `done`
- **related_docs**:
  - `methods/time_integration/implementation.md` (§一様体積力と質量流量一定制御)
- **related_plans**:
  - [`turbulence-iddes-sst.md`](turbulence-iddes-sst.md) — 周期丘 DDES 検証 (本機能の最初の利用者)
- **created**: `2026-07-21`
- **owner**: `CFD Dev`

## 1. 目的

周期丘 (case/39, Re=10595) を文献 (Fröhlich 2005 / Breuer 2009 / ESAIM 2007 DES) と同一の駆動方式
「x/z 周期 + 質量流量一定の動的体積力 (Benocci & Pinelli 1990)」で回せるようにする。
既存 `bodyForce` は config 固定値のみで、断面が x 依存し抗力が先験的に不明な周期丘には不十分。

## 2. スコープ

- **やる**: `bodyForceCtrl` (0/1)、`bodyForceCtrlTarget` (体積平均 ρu_x 目標)、`bodyForceCtrlRelax` (緩和 γ)。
  x 成分のみ制御。毎物理ステップ 1 回の deadbeat 更新 (`f_x^n = f_x^{n-1} + γ(M_t−2M^n+M^{n−1})/(VΔt)`)。
  `bodyforce_history.csv` 出力。`unsteady:1` 必須 (それ以外 throw)。
- **やらない**: y/z 成分の制御、断面流量の直接積分 (体積平均で等価)、仕事項の Jacobian 対角。

## 3. 関連 docs と前提

- 更新則・等価性の理論は `methods/time_integration/implementation.md` に記載済み (本 plan と同時作成)。
- 既存 `bodyForce` 実装 (`bodyForce_d.cu`, 運動量+エネルギー仕事項) はそのまま流用し、
  `cfg.bodyForceX` をステップ境界で書き換えるだけの疎結合とする。
- 既定 `bodyForceCtrl: 0` で既存ケースとビット不変。

## 4. 設計方針

- 制御量 $M=\int\rho u_x dV$ = `thrust::inner_product(roUx, volume)` (1 スカラー/step の D2H、コスト無視可)。
- フックは `advanceOneStep` 冒頭 (物理ステップ境界)。dual-time サブ反復中は $f_x$ 固定。
- 状態 ($V$, $M^{n-1}$, CSV stream) は `bodyForce_d.cu` 内 file-static (単一ソルバインスタンス前提)。
- 陰解法は既存どおり明示ソース (ステップ内定数につき Jacobian 不要)。

## 5. 実装ステップ

1. `solverConfig.hpp/cpp`: キー 3 つ追加 + `bodyForceCtrl=1 && unsteady=0` の拒否。
2. `bodyForce_d.cu/.cuh`: `bodyForceCtrlUpdate_wrapper` (reduction + 更新則 + CSV)。
3. `main.cpp`: `advanceOneStep` 冒頭にフック。
4. 検証: ①既定 off ビット不変 (ビルドのみ)、②case/39 定常 RANS→dual-time 引き継ぎで
   $M/M_t \to 1$ に整定し $f_x$ が抗力と釣り合う値へ収束すること。

## 6. 検証

- case/39.periodic_hills スモーク: `bodyforce_history.csv` で M/M_t の整定 (数 FT 以内、振動発散しない) を確認。
- γ=1.0 で不安定なら 0.5 に落とす (履歴 CSV で判定)。

## 変更ログ

- 2026-07-21: plan 作成、methods 更新と同時。実装着手。
- 2026-07-21: 実装・検証完了 (`case/39.periodic_hills/run_0003_ddes_smoke`)。SST-DDES dual-time 2000 step で
  **M/Mt = 1.00000 ± 7e-6、fx は 932±117 Pa/m に整定** (deadbeat γ=1 はノイズ増幅 1/dt で NaN 級の暴れ →
  **γ=0.02 で運用**、ζ≈0.07 の弱減衰リンギングは 600 step で ±6500→±117 に減衰)。非有限 M で更新を
  スキップする NaN ガードを追加。教訓: 制御量測定に対する deadbeat の実効ゲインは 1/(V·dt) で dt に
  反比例するため、音響スケール dt では γ≪1 が必須。ゲイン分離 (γP/γD 独立) は必要になったら (別 plan)。

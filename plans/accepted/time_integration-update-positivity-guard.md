# 陰的更新の正値性ガード (commit 時の局所 under-relax)

## メタ

- **area**: `time_integration`
- **status**: `done` (負の結果 — opt-in 残置)
- **related_docs**:
  - [`methods/time_integration/implementation.md`](../../methods/time_integration/implementation.md) (「陰的更新の正値性ガード」節)
- **related_plans**:
  - [`time_integration-lowmach-preconditioning.md`](../accepted/time_integration-lowmach-preconditioning.md) (lowMachPrecond=2 との併用を検証)
- **created**: `2026-09-02`
- **owner**: ユーザ指示 (2026-09-02)

## 1. 目的

NS 陰解法 (block-DPLUR) の cfl_pseudo 上限は、case/45 スイープ (run_0015_cflsweep) で
「**下流の厚い低圧 BL で陰的更新が P をアンダーシュート → pMin EOS 床で洗浄 → NaN**」が
律速と特定された。現行の回避策 `implicitRelax` は全セル一律の減速で、床に当たるのは
~121k 中 20-60 ノードに過ぎない。**commit 時にセル単位で更新をスケールする正値性ガード**
(FUN3D/VULCAN 系の solution update limiting と同系) を入れ、全域減速なしで CFL 上限を
引き上げる。lowMachPrecond=2 (壁近傍低速セルの Δτ' を拡大し発散を加速する) との併用も
ガードで救えるかを検証する。

## 2. スコープ

- **やる**: `applyBlockImplicitCorrection_d` (定常 commit) と同 InPlace (dual-time commit)
  にガードを追加。config キー `updateGuardAlpha` (time.deltaT, 既定 0.0=OFF、>0 で有効)。
  既定 OFF でビット同一 (分岐でガード径路を完全に迂回)。
- **やらない**: SST point-implicit 更新のガード (ω 爆発は症状と特定済み) / scalar-DPLUR
  (blockDPLUR=0) 経路 / RHS・Jacobian の変更。

## 3. 設計方針

commit `q ← qN + Δq` の直前に、セルごとに縮小率 $s\in\{1,\tfrac12,\dots,\tfrac1{32}\}$ を

$$\rho(q_N+s\Delta q) \ge \alpha\,\rho(q_N),\qquad
  e_i(q_N+s\Delta q) \ge \alpha\,e_i(q_N),\quad
  e_i \equiv \rho e - \tfrac12|\rho\mathbf u|^2/\rho$$

を満たす最大の半減列値として選び、$q \leftarrow q_N + s\Delta q$ で commit する
(α = `updateGuardAlpha`、推奨 0.5 = 「1 step で密度・内部エネルギーを半分未満にしない」)。

- **P でなく (ρ, e_i) でガードする理由**: CPG では $P=(\gamma-1)e_i$ で等価、TP では
  EOS 反転 (T Newton) を避けつつ P 崩落 (=e_i 崩落) を捕まえられる。EOS 非依存で両ガス共通。
- 5 回半減しても満たせない場合は $s=1/32$ で commit する (EOS 床は最後の防波堤として残す)。
- 既に $\rho\le0$ or $e_i\le0$ のセルはガード対象外 ($s=1$) — 病的状態の隠蔽をしない。
- Δq 配列は書き換えない (診断・他読者への影響なし)。α=0 (既定) は分岐で完全迂回し
  ビット同一。

## 4. 検証

1. 回帰: α=0 で既存 probe (cpg_cfl2 相当) と残差ビット同一。
2. case/45 スイープ再走: ガード ON (α=0.5) × {lowMachPrecond 0, 2} × cfl {4, 8, 16, 32}
   warm restart 2000 step。素の上限 (TP2/CPG3) と implicitRelax 系 (積≈6-8) に対し
   どこまで伸びるか。
3. 収束レース (coarse IC 4000 step): ガード ON 最速設定 vs cfl8+relax0.7 (現最速)。

## 5. 結果 (2026-09-02) — **負の結果: ガードは CFL 上限を上げない**

run_0015_cflsweep (probes5/6, 13 本):

- **guard 単独 (α=0.5, lmp0)**: cfl 4/8/16/32 全て発散 (無ガード比で発散は遅延 97→126 step)。
  NaN ダンプで P は依然 pMin 床に到達 — α=0.5 は「毎 step 半減」を許すため、持続的な
  押し下げには床までの歩行を遅らせる効果しかない。
- **lowMachPrecond=2 + guard**: cfl 4/8/16 全て 9-17 step で発散。**P は 1237 Pa 止まりで
  床に一度も触れず**、ω が 2.8e22 まで爆発 (roK 成長 1e17) — lmp2 の破綻は正値性経路ではなく
  **ガード対象外の SST 隔離更新チャネル**で起きる。併用案は不成立。
- **機構の改訳**: cfl_pseudo 上限の真因は「陰的 defect-correction (近似 Jacobian + SST 隔離
  更新) の反復不安定」であり、EOS 床 NaN はその終端症状。症状のクリップ (本ガード) では
  成長そのものは止まらない。**全域 `implicitRelax` が効くのは反復のスペクトル半径を縮める
  正しいノブだから** (cfl×relax ≈ 6-8 飽和 = 反復の安定境界)。
- **回帰**: guard OFF (既定) の残差軌道差は同一バイナリ 2 回走の再現ノイズ (3-5.6% 相対、
  warm プラトー上の軌道ばらつき) と同水準 = 挙動無変化。
- **処置**: コードは opt-in (`updateGuardAlpha`, 既定 0=OFF) のまま残置 (無害・発散遅延の
  診断的価値はある)。生産推奨は据え置き: **cfl 8 + implicitRelax 0.7 (CPG) / cfl 4-6 +
  0.7 (TP), lowMachPrecond 0**。真の上限引き上げは LHS の整合性改善 (2次/リミッタ項の
  Jacobian 反映、壁法線 line-implicit 化、SST 連成陰化) が本丸 (future work)。

## 変更ログ

- 2026-09-02: 起票。
- 2026-09-02: 実装 (commit カーネル 2 種 + `updateGuardAlpha`)・検証 13 probe。負の結果で
  done。lowMachPrecond=2 併用も不成立 (SST チャネル)。

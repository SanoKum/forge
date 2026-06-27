# SST 生産項 Kato–Launder 補正 (`katoLaunder`)

## メタ

- **area**: `その他 (turbulence)`
- **status**: `done`
- **related_docs**:
  - `methods/turbulence/theory.md` (§7.5)
  - `methods/turbulence/implementation.md` (§5.3.1)
- **related_plans**: なし
- **created**: `2026-06-12`
- **owner**: `CFD Dev`

## 1. 目的

strain-based SST 生産 $P_k=\mu_t S^2$ は、非回転の強加速場 (ノズル喉中心線・よどみ点) で
偽の乱流を生む (stagnation/round-jet anomaly)。これは `dilatationCorrection` (等方分のみ除去)
では消えない別系統の欠陥。Kato–Launder $P_k=\mu_t S\Omega$ を opt-in で追加し、中心線の偽生産を
除去しつつ壁 BL の本物の乱流を保つ。

## 2. スコープ

- **やる**: `rans_sst_source_d` に渦度 $\Omega=|\boldsymbol\omega|$ を組み、`katoLaunder==1` で
  $P_k, P_\omega$ の $S^2$ を $S\Omega$ に置換。config フラグ `katoLaunder` (既定 0) 追加。
- **やらない**: realizability (Durbin) 制限、Kato–Launder 以外の anomaly 対策、LES への適用。

## 3. 関連 docs と前提

- 理論: [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §7.5。
- 実装: [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §5.3.1。
- 前提検証: case 29 conical で軸中心 $k$ 過大を確認 ($S/\Omega\approx24$)。**3D 90° セクタでも
  再現**することを先に確認し、軸対称定式化由来でなくモデル由来と切り分け済み
  ([`case/29.bell_vs_conical/`](../../case/29.bell_vs_conical/README.md))。

## 4. 設計方針

- 渦度大きさ二乗 (既存速度勾配から):
  `Om_sq = (dUxdy-dUydx)^2 + (dUxdz-dUzdx)^2 + (dUydz-dUzdy)^2` $= 2\Omega_{ij}\Omega_{ij}$。
  軸対称(無旋回)でフープは渦度を持たないため補正なし (ひずみ $S^2$ の `+2 S_tt^2` と非対称)。
- `dilatationCorrection` 適用後の `S_sq` を使い、`katoLaunder==1` のとき生産係数を
  `S_prod = sqrt(max(S_sq,0)) * sqrt(Om_sq)` に置換。`Pk = mu_t_eff * S_prod`,
  `Pw = alpha*rho*S_prod`。`==0` では従来の `S_sq` をそのまま使い**ビット一致**。
- リミタ `min(Pk, 10 β* ρ k ω)` と (B) 等方項は後段で従来どおり。

## 5. 実装ステップ

1. `input/solverConfig.{hpp,cpp}`: `int katoLaunder = 0;` 追加 + `turbulence.katoLaunder`
   パース (0/1 検証)。
2. `cuda_forge/ransSource_d.cu`: kernel に `int katoLaunder` 引数追加、`Om_sq` 計算、
   生産係数置換。wrapper で `cfg.katoLaunder` を渡す。
3. ビルド (`build-verify`)。

## 6. 検証

- **ビルド**: `cmake --build build-verify`。
- **検証ケース**: `case/29.bell_vs_conical` の conical (軸対称 2D `run_0017` 相当 + 3D
  `run_3d_conical_rans`)。
- **判定基準**:
  - `katoLaunder:0` 再現 (軸中心線 μt/μlam~38) を基準に、`1` で**中心線 μt/μlam が大幅低下**
    (目標 <5)、**壁 BL の μt ピーク (~26) と推力 (mdot 一致・λ) が不変**。
  - `katoLaunder:0` が既存 RANS 結果とビット一致 (回帰なし)。

## 7. 影響範囲

- ファイル: `input/solverConfig.{hpp,cpp}`, `cuda_forge/ransSource_d.{cu,cuh}`。
- 既存ケース: `katoLaunder` 未指定で従来挙動 (影響なし)。
- docs: `methods/turbulence/{theory,implementation}.md` 更新済み、`methods/index.md` は既存項目内。

## 8. 完了条件

- [x] `methods/turbulence/theory.md` 更新済み (§7.5)
- [x] `methods/turbulence/implementation.md` 更新済み (§5.3.1)
- [x] 実装・検証完了 (§6)
- [x] `.github/plans/README.md` を `done` に更新
- [x] 本 plan の `status` を `done` に変更し §9 に変更ログ

## 9. 変更ログ

- `2026-06-12` — 初稿。case 29 で軸中心 $k$ 過大を診断、3D 検証中。docs §7.5 / §5.3.1 追記。
- `2026-06-12` — **3D 90° セクタ検証完了**: 軸中心の鋭いスパイクは **3D では再現せず** (第1セル/核
  ~0.8, 径方向平坦)。2D 軸対称の中心線スパイクは軸面積 $2\pi r\to0$ の封じ込め由来の増幅と判明。
  ただし**核全体の strain-dominated な $\mu_t$ 過大 (アノマリー本体) は 2D/3D 共通=幾何非依存**
  ($S/\Omega\sim100$, 渦度は物理的に小)。→ Kato–Launder は依然有効な opt-in 修正として採用。
- `2026-06-12` — **訂正 (重要)**: その後の grid-refinement 検証で、case 29 の軸中心スパイクの
  **主因は本アノマリーではなく軸対称 FV の不整合 (planar 面積 vs $r$ 重み体積) の数値特異点**と判明
  (軸第1セル $k$ が細分で 17→836→1956 と発散、背景/BL は収束、3D で非再現)。KL は症状を緩和するが
  根治ではない。根治は [`architecture-axisym-faceweight.md`](architecture-axisym-faceweight.md)。
  KL 自体は一般 SST 欠陥への有効な opt-in 修正として残す。
- `2026-06-12` — **実装・検証完了**: `katoLaunder` (既定 0) 追加
  (`solverConfig.{hpp,cpp}`, `ransSource_d.cu`)。`Om_sq=|ω|²` を速度勾配から組み
  `S_prod=√(S²·Ω²)` で生産置換。**検証** (`run_0017` off vs `run_0018_conical_kl` on):
  軸 $k$ 17.0→1.93, 核 $k$ 5.0→0.47, **壁 BL μt/μlam 13.0→12.9・推力 1953.0→1953.1 (mdot/λ 不変)**。
  `katoLaunder:0` は `S_prod=S_sq` で従来とビット一致 (回帰なし)。

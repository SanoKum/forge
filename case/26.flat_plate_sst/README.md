# 26. 乱流平板 (SST 壁法則検証)

Menter SST k-ω モデルが乱流平板の**壁法則 (law of the wall)** を再現できるかを確認するケース。
陰解法 (block DPLUR, `timeIntegration: 11`) + SLAU + MUSCL で定常計算する。

## 流れ条件

- 入口: 全圧 `Pt=100000 Pa`, 全温 `Tt=288.15 K` (`inlet_Pressure`)
- 出口: 背圧 `Ps=97250 Pa` (`outlet_statPress`)、`M≈0.2` (U∞≈67.8 m/s)
- 空気 `mu=1.8e-5 Pa·s`, `Re/m≈4.46e6` → 平板長 1 m で `Re_L≈4.5e6` (乱流)
- 流入乱流: `k=0.3, omega=300` (freestream `mu_t/mu≈60`, naca_ml 相当)

## 領域・メッシュ (`mesh/flat_plate.geo`)

- `x∈[-0.1,1.0]`, `y∈[0,0.2]`, z 厚さ 0.01 (1 層押し出し)
- 底面: `x∈[-0.1,0]` slip (対称), `x∈[0,1.0]` no-slip wall (平板)
- 上面 slip, 入口/出口は上記 BC、z 両面 slip
- 壁垂直: 90 点, 第一セル高さ **4.14e-6 m**, 等比 1.10 → **y⁺₁≈0.35〜0.4 (<1)**
- streamwise: 前縁 (x=0) へ上流側・平板側の両方を寄せ、前縁でセルサイズ連続 (≈0.4 mm)

## 実行手順 (段階法)

**重要**: 定常計算の実効 CFL は `cfl` ではなく `cfl_pseudo`
(`.github/forge-solver-settings.md` の「CFL の定義」を参照)。
本ケースは壁セルが極めて細く陰解法でも CFL 制限が低いため、段階的に進める。

1. **stage A (`run_0005_slau_spinup`)**: 一様初期値からの cold start。
   `convMethod: 0` (1次風上) + `cfl_pseudo: 2.0` で 25000 step 回し、平滑な発達場を作る。
2. **stage B (`run_0006_slau_muscl`)**: stage A の `res_25000.h5` から restart。
   `convMethod: 2` (MUSCL) + `cfl_pseudo: 3.0` で 40000 step 回し、2次精度で収束させる。

> cold start で `cfl_pseudo≳5` や、いきなり MUSCL + 高 `cfl_pseudo` にすると発散する。
> 発達場からの restart でも本ケースの安定上限は `cfl_pseudo≈3` 程度。

メッシュ生成・変換・実行は `.github/forge-calculation-workflow.md` の標準手順に従う
(Docker, `convertGmshToForge`, `forge`)。

## ポスト処理

```bash
python3 tools/postprocess_wall_law.py run_0006_slau_muscl 0.3 0.6 0.89
```

`wall_law.png` (u⁺-y⁺ と Cf-Re_x) を出力する。

## 収束について

陰解法 `blockDPLUR` は減衰付き point-Jacobi の defect correction で、近似ヤコビアン+高アスペクト比
壁セル+陽的 SST 生産項のため安定 `cfl_pseudo` は MUSCL で ≈5(cold start は 2〜3)が上限。
そのため収束が遅く、`run_0006`(40000 step)時点では外層がまだ発達途中(Cf が 10000 step あたり
3〜4% ドリフト)だった。

- **stage C (`run_0007_slau_muscl_innersweep`)**: `run_0006/res_40000.h5` から restart し、
  `cfl_pseudo: 5.0`, `nStepInner: 20`, MUSCL で **120000 step** 継続。残差 3.8e-7、
  最後の 20000 step での Cf ドリフト 0.1% 以下まで収束。**最終結果はこの run を使う。**

> 安定 CFL が低い構造的原因と改良案(SST 生産項の陰的化, 真の粘性ヤコビアン)は
> `.github/plans/implicit-acceleration-session-prompt.md` に整理(別セッションで実装予定)。

## 結果 (run_0007, 収束済み MUSCL)

| station | Re_x | y⁺₁ | u_τ | Cf | Cf/Schlichting |
|---|---|---|---|---|---|
| x=0.30 | 1.34×10⁶ | 0.36 | 2.707 | 3.12×10⁻³ | 0.885 |
| x=0.60 | 2.66×10⁶ | 0.34 | 2.580 | 2.82×10⁻³ | 0.918 |
| x=0.89 | 3.99×10⁶ | 0.33 | 2.511 | 2.67×10⁻³ | 0.944 |

- u⁺ vs y⁺: 3 ステーションが普遍プロファイルに重なる。粘性低層 `u⁺=y⁺` (y⁺<5)、
  対数則 `u⁺=(1/0.41)ln y⁺+5.0` (y⁺=30〜300)、外層 wake までよく再現。
- Cf vs Re_x: 発達領域で Schlichting 1/5 乗則と 1/7 乗則の間に乗る。1/5 乗則より数% 低いのは、
  前縁からの層流→遷移区間ぶんの virtual-origin オフセット(Re_x を前縁起点で測るため)による系統差。
- y⁺₁≈0.33〜0.36 (<1) を全域で満たす。

## 計算 run 一覧

| run_* | 目的・設定 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_rans_implicit` | 標準: SLAU + SST + block-DPLUR 陰解 (2D) | rms_ro 2.08e-9 収束、壁法則再現 | active |
| `run_regr_cf` | 回帰: 閉形式 FVS 既定化 (`implicitSolvePrecision=0`) の確認。run_0001 と同条件 | 残差収束が legacy と一致 (rms_ro 2.0e-9)、NaN なし → **閉形式は 2D 陰解に無害**。plan [precision-mixed-axisym.md](../../.github/plans/precision-mixed-axisym.md) | active |

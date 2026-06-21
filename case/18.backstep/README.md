# case/18.backstep — 後退段 (Backward-Facing Step)

ステップ後方の剥離・再循環・再付着を扱う基本剥離流ケース。RANS/LES の速度評価
(`SPEEDUP_CAMPAIGN.md`) と **SST-DES (DDES) の機能検証 (T1-B)** に用いる。

## メッシュ

| ファイル | 次元 | 要素数 | 用途 |
| --- | --- | --- | --- |
| `mesh/backstep_2d.msh` | 2D | 32,712 | RANS・速度評価 (スパンなし) |
| `mesh/backstep.msh` | 3D (スパン 4H, 80 層) | 923,200 (→ 857,600 cells) | **DDES**: スパン方向解像が必要 |

- 物理面: `inlet`(1) / `outlet`(2) / `top`(3, wall) / `bot`(4, wall) / `side2`(5) / `side1`(6)。
- DDES では `side1`/`side2` を **周期 BC** にしてスパン方向の解像乱流を許す (2D・slip では RANS 同然)。
- 3D メッシュ品質: `check_mesh_quality.py` VERDICT **PASS** (AR max 2.3, skew 0.333)。

## 流れ条件

低マッハ (M≈0.1 級) の非圧縮的後退段。NASA/TMR 2DBFS (Driver & Seegmiller 1985) を
定量参照とする (Re_H≈36,000, 再付着 x_R/H=6.26±0.10)。定量比較は Phase 1.5 以降。

## 計算 run 一覧

> 速度最適化キャンペーン (`run_0001`〜`run_0034`) の詳細は [`SPEEDUP_CAMPAIGN.md`](SPEEDUP_CAMPAIGN.md)。
> ここでは代表 run と DDES 検証 run を示す。

| run_* | 目的・主要設定 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_slau_rans_implicit` | 2D RANS SST + block-DPLUR 陰解 (`backstep_2d.h5`) | 収束 RANS 場 (再循環泡)。DDES restart 元 | active |
| `run_0002_slau_rans_rk3` | 2D RANS RK3 陽解 | 速度評価基準 | ref |
| `run_0007_slau_2ndup_baseline` 〜 `run_0034_slau_K5_final` | 速度最適化キャンペーン (2D) | カーネル最適化の段階計測 ([SPEEDUP_CAMPAIGN.md](SPEEDUP_CAMPAIGN.md)) | ref |
| `run_0035_des3d_ddes` | **SST-DES T1-B 機能検証 (3D)**: `backstep.h5`(857k, 周期スパン)、`DESmode:1`、2D RANS 場から `interp_field` cross-mesh restart、500 step implicit | **NaN なし**。f_d: 付着近壁 BL≈0.02 (frac<0.1=0.96)、剥離せん断層 (x∈[3,9))≈0.67=LES 活性、再循環泡内部 (Ux<0,高νt)≈0=RANS、再付着~6H。Δmax≈0.13、rd∈[5e-4,10]。`residual_history.png` | active |

> **SST-DES (DDES) T1-B 機能検証**: f_d が剥離域で立ち付着 BL で 0、NaN なし、を満たせば Phase 1
> 合格。**収束 VERDICT は NOT CONVERGED が想定どおり** (DDES は本質的に非定常で残差は下げ止まらない)。
> 解像乱流の定量化 (RMS・スペクトル −5/3・x_R/H) は Phase 1.5 (DES 用低散逸 flux) 後の仕事。
> 計画: [turbulence-iddes-sst.md](../../.github/plans/turbulence-iddes-sst.md)。

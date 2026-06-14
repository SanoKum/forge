# case/16.nozzle_wys — Wyslouzil 超音速ノズル

Wyslouzil et al., *J. Chem. Phys.* **113**, 7317 (2000)
("Binary condensation in a supersonic nozzle") の超音速 Laval ノズル。
凝縮実験用の矩形断面ノズルで、forge による 2D / 3D・非粘性 / 粘性 / 乱流の計算を行う。

## ジオメトリ (`mesh/nozzle_H.geo`, 単位 cm → ScalingFactor 0.01 で m)

論文との寸法照合 (おおむね一致):

| 量 | `.geo` 由来 | 論文 |
| --- | --- | --- |
| 入口全高 | 2×1.236 cm ≈ **24.7 mm** | flow straightener 25.4 mm |
| 断面幅 (z 押し出し) | 1.259 cm = **12.59 mm** | 12.7 mm |
| スロート半高さ | ≈ 0.225 cm = **2.25 mm** | — |
| 出口/スロート面積比 | **1.69** → 等エントロピー M≈2.0 | 設計 M≈1.95 (γ=1.4) |

矩形断面のため 2D 平面近似が妥当 (ユーザ確認済)。

## 条件

- 入口: `inlet_Pressure`, Pt=101325 Pa, Tt=293.15 K (亜音速圧力入口)。
- 出口: `outlet_statPress`, Ps=2000 Pa (超音速流出のため実質外挿)、逆流用 Pt/Tt 設定。
- 気体: γ=1.4, cp=1039 J/kgK (N2 相当)。時間積分は陰解法 (block-DPLUR, `timeIntegration:11`)。
- 初期場: `initial: nozzle_wys` (一様 M≈0.3, Pt/Tt そろえ。`setInitial.hpp` に追加)。
- 段階起動: 非粘性(slip) → 粘性層流(no-slip) → 乱流(SST) を引き継ぎ計算で段階的に。

### 2D / 3D メッシュの壁の扱い

- **2D** (`nozzle_2d*.geo`, 1層スラブ): Euler は全面 slip。粘性/乱流は輪郭壁=no-slip、
  front/back=slip (対称面) とするため `nozzle_2d_visc.geo` で front/back を別 physID(4) に分離。
- **3D** (`nozzle_3d.geo`, 12層): 矩形ダクトの4壁すべてが実壁。Euler は slip、粘性/乱流は全壁 no-slip。
  ※ **乱流(SST)では mesh を no-slip 壁付き bcond で `convertGmshToForge` する**こと
  (壁が slip だと `wall_dist=0` になり SST の ω が発散する)。

## 計算 run 一覧

すべて SLAU + 陰解法 (block-DPLUR)。中心線 exit Mach は等エントロピー値 2.00 を基準に比較。
詳細・後処理図は各 run の `centerline_compare.png` / `residual_history.png`。

| run | 物理 | 壁 | exit M (中心線) | exit Ps | 状態・備考 |
| --- | --- | --- | --- | --- | --- |
| `run_0001_slau_2d_imp` | 2D 非粘性 Euler | slip | **1.989** | 12.9 kPa | active. 等エントロピー 2.00 と一致 (0.6%) |
| `run_0002_slau_3d_imp` | 3D 非粘性 Euler | slip | **1.990** | 12.9 kPa | active. 2D と完全一致 (検証) |
| `run_0003_slau_2d_visc` | 2D 粘性 層流 | 輪郭no-slip / fb slip | **1.940** | 13.9 kPa | active. BL 変位で M 低下 |
| `run_0005_slau_3d_visc` | 3D 粘性 層流 | 全壁 no-slip | **1.937** | 14.0 kPa | active. 2D 層流とほぼ同等 |
| `run_0004_slau_2d_sst` | 2D 乱流 SST | 輪郭no-slip / fb slip | **1.883** | 15.3 kPa | active. 乱流 BL でさらに M 低下。**乱流熱伝導修正後 T≤Tt** |
| `run_0006_slau_3d_sst` | 3D 乱流 SST | 全壁 no-slip | **1.534** | 25.4 kPa | active. 4壁 no-slip の閉塞で M 大幅低下。T≤Tt |

exit Mach の単調序列 (非粘性 1.99 > 層流 1.94 > 2D乱流 1.88 > 3D乱流 1.53) は
境界層変位による有効面積減少として物理的に整合。中心線静温はいずれも全温 293 K 以下。

## 既知の課題

- 近壁メッシュは Euler デモ用で y+~1 の壁解像ではない。乱流 BL を定量評価するには
  壁直交方向を細分した専用メッシュが必要 (中心線コア量は妥当)。
- 本ケース検証中に **RANS エネルギー方程式の乱流熱伝導欠落バグ**を発見・修正
  (`.github/plans/diffusion-turbulent-thermal-conductivity.md`)。修正前は近壁静温が
  449 K (全温 293 K 超過)、修正後 293 K に収束。

## 後処理スクリプト

- `postproc_centerline.py <run> <label>` — 中心線 Mach / 静圧を面積比等エントロピーと比較。
- `plot_residuals.py <run> <label>` — 全 rms 残差の片対数プロット。
- `build_restart.py <src_res> <src_mesh> <dst_mesh> [--rok K --roomega W]` — 引き継ぎ初期場生成。

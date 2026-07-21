# 29.bell_vs_conical — Rao TOP ベルノズル vs 等軸長コニカル (軸対称・推力比較)

RRS [*Making Correct Parabolic Nozzles*](https://rrs.org/2023/01/28/making-correct-parabolic-nozzles/) (2023-01-28)
の作図法に従って**滑らかな Rao 推力最適化放物線 (TOP) ベルノズル**を生成し、
**同じ軸長**のコニカルノズルと推力を比較する軸対称ケース。

23/27 番は「最小長さ (MOC) ノズル」でスロート直後に中心膨張扇 → 反射で**圧縮波**が出ていた。
本ケースは接線連続な Bézier 壁で発散部を構成し、これを避ける。

## ノズル作図法のまとめ (RRS 記事)

記事の核心は **「放物線は *canted* (傾けて) 描かねばならない」**。
軸並行の放物線 $y^2 = 4px$ は、ベル発散部に課す 6 条件
(始点 $N$ の $x,y$・終点 $E$ の $x,y$・$N$ での接線 $\theta_n$・$E$ での接線 $\theta_e$)
に対し自由度が 5 しかなく**過拘束で破綻**する。記事は放物線を回転させた
パラメトリック式 (式 7.1-7.2) を 6 元連立で数値的に解く。

実装上これと**等価で頑健**なのが **2 次 Bézier 曲線**である
(2 次 Bézier はちょうど「傾いた放物線」)。本ケースはこちらを採用:

### 手順 (上流 → 下流)

1. **収縮部**: 半径 $R_c$ のチャンバー直管 → 収縮半角 $\beta$ (既定 30°) の直線収縮。
2. **スロート上流円弧**: 半径 $1.5\,R_t$。収縮直線とスロート (垂直) に接する。
3. **スロート下流円弧**: 半径 $0.382\,R_t$。スロートから角 $\theta_n$ まで回し、
   放物線始点 $N$ を与える。
   $$N=\big(0.382R_t\sin\theta_n,\; R_t+0.382R_t(1-\cos\theta_n)\big)$$
4. **放物線 (2 次 Bézier)** $N\to E$。終点 $E=(L_n, R_e)$。
   制御点 $Q$ は「$N$ から $\theta_n$ の直線」と「$E$ から $\theta_e$ の直線」の交点:
   $$P(t)=(1-t)^2 N+2(1-t)t\,Q+t^2 E,\quad t\in[0,1]$$
   始点接線 $=\tan\theta_n$、終点接線 $=\tan\theta_e$ が自動的に満たされ、壁は接線連続。

### 設計パラメータの決め方

- **膨張比** $\varepsilon=A_e/A_t=(R_e/R_t)^2$。設計出口 Mach から等エントロピー面積比で決定。
- **$\theta_n,\theta_e$**: Rao(1960) のチャート (記事/ Sutton に再掲) を $\varepsilon$ で読む。
  本実装は 80% ベル基準の代表値を digitize し $\log\varepsilon$ 線形補間 (`--thetan/--thetae` で上書き可)。
- **ベル長** (15° コニカル比 $K\%$):
  $$L_n=K\cdot\frac{(\sqrt{\varepsilon}-1)R_t}{\tan 15°}$$
  80% ベル = 同膨張比の 15° コニカルの 0.8 倍長さ。

## 推力比較の設計 (本ケース)

「**同じ軸長**でベル vs コニカルが推力でどれだけ変わるか」を見るため、両者で
**スロート・チャンバー・出口面積 ($R_e$)・軸長 ($L_n$) をすべて一致**させる。
設計手順 (ユーザ確定方針):

1. **良い C1 コニカルを先に決める**: スロート下流円弧 ($0.382R_t$) を**半角 15°**の直線
   コーンに *接して* つなぐ (円弧終端角 = 15° → **1 階微分連続**, スロートを跨ぐ折れなし)。
   出口を $M_e=4$ ($\varepsilon=10.72,\ R_e=32.74$ mm) で固定し、C1 を満たす軸長
   $L_n=85.37$ mm を解く (= ちょうど $\varepsilon=10.72$ の 15° フルコーン長)。
2. **Rao ベルを同じ $R_t/L_n/R_e$ で作り直す**: 円弧を $\theta_n$ まで回し 2 次 Bézier で
   $(L_n, R_e)$ へ。ベル長は 15° コーンと同じなので **100% 長ベル**。出口角は
   $\theta_e=8°$ (100% 長相当; `--thetae` で可変) を採用。

→ 差は純粋に**発散部の形状**のみ。ベルは出口流が軸方向に近い ($\theta_e=8°$ vs コーン 15°)
ため**発散損失が小さく推力が高い**はず、という古典的結果を定量化する。

推力は出口面の運動量フラックス + 圧力項を軸対称積分:
$$F=\int_{A_e}\big(\rho u_x^2 + (P-P_\infty)\big)\,2\pi r\,dr$$
発散効率 $\lambda = \int\rho u_x^2\,dA / \int\rho|U|u_x\,dA$ も併せて評価する。

## 既定諸元 (`mesh/nozzle_params.txt`)

| 項目 | 値 |
| --- | --- |
| スロート半径 $R_t$ | 10.0 mm |
| 出口半径 $R_e$ | 32.74 mm (両ノズル共通) |
| 膨張比 $\varepsilon$ | 10.72 (設計 $M_e=4.0$, $\gamma=1.4$) |
| 軸長 $L_n$ | 85.37 mm (両ノズル共通) |
| コニカル半角 | **15°** (C1: 円弧に接続) |
| ベル $\theta_n / \theta_e$ | 34.2° / 8.0° (100% 長ベル) |
| 収縮半角 $\beta$ | 30° |
| チャンバー直管長 | 20 mm |

## メッシュ生成方針 (重要)

壁は **厳密プリミティブ** で構成する: 直線部 → `Line`、スロート円弧 → `Circle`、
ベル放物線 → `Bezier(N,Q,E)`。各セグメントを 1 ブロックとする**多ブロック構造格子**で、
壁の折れ (kink) とスロート $T=(0,R_t)$ を必ずブロック角点に置く。

> 初版で点列を Catmull-Rom `Spline` でつないだところ、接合部の重複点と曲率不連続で
> スロート付近が大きく **オーバーシュート** (喉径 10→7.85 mm に潰れ) し、コニカルでは
> 直線がスロートを跨いで **負体積セル** を生んだ。厳密プリミティブ + 多ブロックで解消
> (converter の `CW cells detected: 0`)。これが 23/27 番で出ていた圧縮波の主因でもある。

## 形状生成

```bash
cd mesh
python3 make_nozzle.py --thetae 8           # 既定 (Me=4, 15°コーン, 100%ベル θe=8°)
# 例: チャンバーを伸ばす / コーン半角や膨張比を変える
python3 make_nozzle.py --thetae 8 --chamber 60
python3 make_nozzle.py --cone-half 12 --eps 25 --thetan 32 --thetae 6
```

出力: `contour_{bell,conical}.csv` / `{bell,conical}.geo` / `nozzle_params.txt` / `nozzle_preview.png`。

## 実行手順

```bash
# 正準形状: theta_n=30 (100%長ベルに適), theta_e=8, cone 15deg
cd mesh && python3 make_nozzle.py --thetan 30 --thetae 8                 # inviscid 用 (一様格子)
python3 make_nozzle.py --thetan 30 --thetae 8 --prog-r 0.95 --ny 100      # 粘性/RANS 用 (壁密格子)
cd .. && ./run_case.sh <run_dir> <bell|conical>   # convert + 等エントロピー IC + forge
python3 postprocess_thrust.py        # 出口面推力・発散効率の比較
python3 make_plots.py                # residual / Mach 場 / 出口プロファイル / near-wall vis_turb 図
```

数値共通: SLAU・軸対称、入口 $P_t=4$ MPa/$T_t=1500$ K、出口背圧 20 kPa (等エントロピー出口
静圧 26.4 kPa より低く → 全域超音速で外挿)。IC は均一ではなく**準 1 次元等エントロピー場**を貼って
起動 (均一だと 40 倍圧力比の過渡衝撃で発散)。**2 次精度 (convMethod=1) / cfl 0.3** が安定
(clean mesh でも cfl 0.5 は step 81 で発散; 喉付近の過渡不安定で cfl≤0.3)。

| run | 内容 | 壁 BC | 粘性 |
| --- | --- | --- | --- |
| `run_0006_bell_inv_cl` / `run_0007_conical_inv_cl` | inviscid (発散損失=形状効果を分離) | slip | 0 |
| `run_0004_bell_visc` / `run_0005_conical_visc` | laminar (実摩擦) | no-slip | Sutherland |
| `run_0016_bell_imp_sst` / `run_0017_conical_imp_sst` | RANS SST (陰解法・clean inlet) | no-slip | Sutherland+SST |

> いずれも壁密格子 (`--prog-r 0.95 --ny 100`, 第1セル ~10µm, y+≈0.6)。
>
> **(1) 入口 $\omega$ の適正化 (freestream $\mu_t$ 過大の修正)**: 初期設定 $k=1,\ \omega=1000$ は
> 入口で $\mu_t/\mu_{lam}\approx176$、コアまで $\approx120$–180 と**過大** (非物理的にコアが乱流粘性で
> floodされる)。水力直径の発達管相関 ($\ell=0.07D_h$) は大径チャンバでは逆に $\mu_t/\mu_{lam}\sim80$–1200
> と更に過大で不適。代わりに**物理的 freestream 比 $\mu_t/\mu_{lam}=10$** を狙い、$k=1$ (強度 ~1.7%)
> を保って $\omega=\rho k/(10\mu_{lam})\approx1.8\times10^4$ を入力。結果、入口 $\mu_t/\mu_{lam}$ 8.9、
> 発散部コア 0.5–4 とクリーンになり、壁 BL の乱流 (近壁 $k\sim6500$, $\mu_t/\mu_{lam}$ ピーク ~26) は維持。
>
> **(2) 乱流着火 (turbulent ignition) は遅い**: $k$/$\omega$ 拡散は有効粘性
> $\mu_{lam}+\sigma\mu_t$ ([scalarTransport_d.cu](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu)) で
> 正しく実装されているが、層流的 IC からは $\mu_t\approx0$ で初期は**層流拡散のみ**のため $k$ の
> 立ち上がりが遅い (自己加速するまで潜伏。順圧勾配=喉の再層流化が助長)。陽解法では 12k step で
> 近壁 $k\approx0$ (未発達)、40k で部分発達、~60k で全域発達。
>
> **(3) 陰解法 (block DPLUR, `timeIntegration:11`) で cfl_pseudo≈4**: 発達済場から warm-start し
> 上記 clean inlet で収束。**ベルは cfl_pseudo=4 で安定**、**コーンは cfl_pseudo=2** (4 は step90 で
> 発散; コーンの方が陰解法 CFL 上限が低い)。推力は ±0.01% で収束。本表は陰解法・clean inlet の値。

## 結果 — 推力比較 (同じ喉 / 軸長 / 出口面積)

NaN/発散チェック: 全 run NaN なし、最終場は物理的 ($P>0,\ T>0,\ P_{max}\approx P_t$)。
**mdot はベル/コーンで一致 (≤0.2%)** = 同一チョーク喉の妥当性 (出口面積 $A_e=\pi R_e^2=3367$ mm²)。

| モデル | ベル $F$ [N] | コーン $F$ [N] | **ベル推力ゲイン** | $\lambda$ ベル/コーン | 出口流れ角 ベル/コーン |
| --- | --- | --- | --- | --- | --- |
| inviscid (slip) | 1976.2 | 1962.9 | **+0.68%** | 0.996 / 0.983 | 3.7° / 10.0° |
| laminar | 1978.9 | 1957.2 | **+1.11%** | 0.997 / 0.984 | 3.5° / 9.8° |
| RANS SST (陰解法・clean inlet) | 1971.3 | 1953.0 | **+0.94%** | 0.997 / 0.984 | 3.4° / 9.7° |

- **ベルはコーンより ~0.7–1.1% 高推力** (inviscid/laminar/turbulent で頑健)。出口流がより軸方向
  (3.4° vs ~10°) で発散損失が小さいため。
- **妥当性**: コーンの $\lambda\approx0.983$ は古典理論 $\tfrac12(1+\cos15°)=0.9830$ と一致。
  出口流れ角もコーン $\approx\tfrac23\times15°$ の source-flow 値と整合。
- **摩擦損失は小さい (<0.5%)**: 同一壁密格子の inviscid を基準に viscous-inviscid を取ると、
  laminar はベル +0.1%/コーン −0.3%、発達済 RANS はベル −0.2%/コーン −0.5%。ノズル BL は薄く
  強順圧勾配で部分再層流化するため (favorable pressure gradient relaminarization)。発達済
  乱流で摩擦損失は現実的に増えるが小さく、ベル優位 (+1%) は維持。
- $\theta_n=30°$ (100% 長ベル) で出口 Mach はおおむね均一 ($\theta_n=34°$ ではコア過膨張で
  $M\approx5.7$ の不均一が出ていた)。

## 陰解法 CFL / 前処理スタディ (`sweep_implicit/`)

inviscid ベルで `timeIntegration: 11` (block DPLUR, `blockDPLUR: 1`, `nStepInner: 20`) を
等エントロピー IC から起動し、`cfl_pseudo` と `lowMachPrecond` を掃引:

| | cfl_pseudo=5 | 10 | 20 | 50 | 100 |
| --- | --- | --- | --- | --- | --- |
| precond OFF (`lowMachPrecond: 0`) | 収束 | 遅い | NaN | NaN | NaN |
| precond ON  (`lowMachPrecond: 2`) | NaN | NaN | NaN | NaN | NaN |

- **陰解法の実効 CFL 上限は `cfl_pseudo`≈5** (10 は遅延収束、≥20 で発散)。`case/08.bump` の 50 に
  比べ低いのは、超音速 MUSCL コアが硬いため (case 27 の implicit+MUSCL 発散と整合)。
- **前処理 (`lowMachPrecond: 2`) は全 CFL で即発散**: 低マッハ用の手法で、超音速コアには不適
  (case 27 run_0005 と同じ)。
- **RANS SST でも陰解法は有効** (`run_0016`/`run_0017`): 発達済場から warm-start し、ベルは
  **cfl_pseudo=4 で安定収束**、コーンは **cfl_pseudo=2** (4 は発散、コーンの陰解法 CFL 上限が低い)。
  陽解法 (cfl 0.3) より大幅に大きい擬似時間刻みで収束。
- cf. [[implicit-blockdplur-config]] / [[nozzle-implicit-cfl-ceiling]]。

主要成果物: `mach_comparison.png`, `exit_profiles.png`, `rans_visturb_profile.png`,
各 `run_*/residual_history.png`, `mesh/nozzle_preview.png`。

## 軸中心 k 過大の調査 — **機構特定済み: 二層構造** (経緯の誤診は訂正済)

詳細・経緯は plan
[`architecture-axisym-axis-singularity.md`](../../plans/accepted/architecture-axisym-axis-singularity.md)。

### 層 (a): 実用格子 (軸セル ≥ 0.25mm) の軽度スパイク → **Kato–Launder で解決**

- 平均流は滑らか (Uy 振動なし)。軸第1セル $k$=17 vs 核 4.5 (`run_0017`)。
- 大きな軸方向ひずみ + フープ + 生産リミタ正帰還 + 軸無フラックス滞留の**生産駆動**。
- `katoLaunder: 1` で $k$ 17→1.9、壁 BL・推力不変 (下表)。3D の平坦 (第1セル/核 0.85,
  `run_3d_bf_lowbg`) と整合。

### 層 (b): 軸 ≤10µm の細分で**平均流の数値不安定** (未根治, 実用上は回避)

2 面クラスタ格子 (`run_gx_A`: 壁=軸 9.6µm / `run_gx_B`: 3.2µm, Bump) で:

- **$U_y$ が近軸チェッカーボード (偶奇) 振動**: x=0.51mm 列で +0.00, **−6.77, +7.93**, −2.32,
  +6.26 ... (物理値 ~0.5 m/s)。stored 勾配は場の FD と一致 = 勾配計算は正しく、**場が振動**。
- $(\partial_r U_y)^2\sim10^{11}$ → $S^2$ → $k$ 爆発 (A: $k$~1400–5500 で振動・収束せず;
  B: rms_roK **1.6e11** 爆発 → 近軸 $k$=0 崩壊)。**軸細分で悪化** = 数値不安定。
- inviscid / baseline (軸 0.25mm) は $U_y$ 滑らか → 細分で顕在化する潜在モード。
- 機構: 半径方向は近軸で深い低 Mach ($u_r\to0$) → 風上流束の圧力-速度結合が消え偶奇モードが
  生き残る (古典的チェッカーボード)。3D は軸が内部点でこの構造なし。
- 失敗した対処: `lowMachPrecond: 2` (超音速コアで発散, 既知)・`lowMachPrecond: 1` (細軸で
  step8 NaN)。根治候補 (near-axis 限定の圧力-速度結合回復) は plan §3。

**実用指針**: 軸対称 RANS ノズルは `katoLaunder: 1` + **軸の過細分を避ける** (軸セル ≥ ~0.1mm;
壁解像は `--prog-r` の壁側クラスタで行う)。

### 経緯メモ (誤診の訂正)

①「planar 面積 vs $r$ 重み体積の不整合」→ 誤り (flux 面積・体積は両方 $r$ 重み済,
[`variables.cpp`](../../solver_density_cuda/variables.cpp))。②「拡散勾配の planar 面積が bug」→
誤り (勾配は planar 作用素で正しい)。③ uniform 格子 (`run_gr_ny100/200`) のグリッド試験 →
壁粗で $k$ 非収束のため無効。最終特定は 2 面クラスタ格子での $U_y$ チェッカーボード観測。

### Kato–Launder の検証 (層 (a) への対処)

`katoLaunder:1` ($P_k=\mu_t S\Omega$) は実用格子の層 (a) を解決する (`run_0018_conical_kl`):

  | | 軸 $k$ | 軸 μt/μlam | 核 $k$ | 壁BL μt/μlam (x=50) | 推力 F | mdot/λ |
  | --- | --- | --- | --- | --- | --- | --- |
  | `katoLaunder:0` | 17.0 | 37.0 | 5.02 | 13.0 | 1953.0 | 1.2998/0.9842 |
  | `katoLaunder:1` | **1.93** | **10.7** | **0.47** | **12.9** | 1953.1 | 1.2998/0.9842 |

層 (b) (軸過細分のチェッカーボード) には効かないため、軸セルは ≥0.1mm を保つこと。

> メモ: 3D は convert が $k=\omega=0$ で IC を書くため起動時に $1/\omega$ 発散する。IC で $\omega$ を
> シード (例 18000) すれば安定。butterfly 格子は壁第1層が全周 8µm 一様 (squircle 版は対角で細・
> 対称面で粗だった)。

## 計算 run 一覧

> 本ケースは複数系列の run が同居する (推力比較 = `run_0001_bell`/`run_0002_conical` 系、k スパイク調査 =
> `run_0016〜0018`・`run_3d_*`、**近軸 混合精度の根治 = 下表 = `run_0001_mixed_lam_slau`〜`run_0015_prod_double`**)。
> 番号が衝突する別系列があるので、目的は slug と本表で判断すること。

### node 粘性壁の近壁振動 / no-slip クロージャ調査 (2026-06-20, plan [diffusion-node-wall-viscous-distance.md](../../plans/accepted/diffusion-node-wall-viscous-distance.md))

node-centered viscous の壁近傍振動の真因切り分け。全 run **explicit RK3・cfl=0.1・convMethod=0・bndFirstOrder=1・
laminar・axisym・conical**、参照と同一 IC/BC。**結論: 真因は粘性壁 flux 距離ではなく壁∩出口コーナーの矛盾ゴースト
(質量/エネルギー非流出→roe 先頭発散)。** dcc 退化は twall 出力ノイズの原因ではあるが発散の主因でない。
バルク解 (中心線 Mach) は元スキームで Euler/cell と一致し健全。既定挙動は不変 (`nodeWallViscGradFlux=0`)。

| run_* | 目的・設定差分 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_0020_visc_node_wallghost_ab` | 旧 mirror-ghost (基準, `nodeWallViscGradFlux=0`) | 安定。twall_x TV/range=14.8・符号反転 172・スパイク −1.77e4。中心線 Mach=Euler/cell。**基準** | active |
| `run_0019_visc_node_wallgradflux` | 純勾配 ∇u·S 法線項 | 発散 (~13k, 断熱熱漏れ) | 破棄予定 |
| `run_0021_visc_node_wallgradflux_adiabfix` | 純勾配 + 断熱熱流束=0 | 壁せん断 0.37Pa (cell の 1e4 倍過小=no-slip 喪失) | 破棄予定 |
| `run_0022_visc_node_wallsmoothdist` | 滑らか距離 2·vol/sss | 発散 (~9.8k, ドラッグ弱) | 破棄予定 |
| `run_0024_visc_node_normdrag_only` | 滑らか距離+転置項除去 | 発散 (~9.8k, 転置無関係) | 破棄予定 |
| `run_0025_visc_node_floor005` | dcc floor `max(dcc,0.05·vol/sss)` | 安定だが twall ノイズ残存 (TV/range=15.8) | 破棄予定 |
| `run_0023_visc_node_walldirichlet` | nodeWallDirichlet (残差射影) | 発散 (~22.6k, roe 先頭) | 破棄予定 |
| `run_0026_visc_node_cleandirichlet` | クリーン Dirichlet (毎ステージ KE 除去射影) | 全域崩壊 | 破棄予定 |
| `run_0027_visc_node_dirichlet_freeze` | クリーン Dirichlet+壁ノード全保存量凍結 | 発散 **3 倍遅延** (~68.5k)→**コーナー積算が主因と確認** | 破棄予定 (診断証拠) |

> 診断: `FORGE_VISC_WALL_DIAG=1` で壁半割面の `dn`/`dcc`/接線オフセットを集計 (退化 `|dn|∈[1e-8,3e-4]` 4〜5 桁変動を実証)。
> 正攻法は plan [discretization-node-boundary-ghostless.md](../../plans/active/discretization-node-boundary-ghostless.md) Phase 2 のコーナー面ごと BC。

#### SU2 流ゴーストレス再設計の検証 (2026-06-20, plan §7-8)

| run_* | 設定 | 結果 |
| --- | --- | --- |
| (旧 nozzle.h5 27472 + 現 HEAD) | single-owner, 元壁 | 安定 (cfl 0.50) — 現 HEAD は旧メッシュを正常に回す |
| (fresh conical.msh 34138, ny=80) | single-owner / multi-marker 双方 | **step0 NaN→崩壊** = fresh-fine メッシュ独立の不安定 (私の変更と無関係) |
| `run_0030_coarse_ghostless_su2` | 粗 conical (ny=40,13858) + multi-marker + nodeWallDirichlet=1 | **40k 安定**, P≈Pt, T 1827(残 overshoot)。ただし粗は single≡multi (コーナー未分割) で multi-marker 実質無効=「粗+Dirichlet」由来の安定 |

> **重要**: 当初 plan の「multi-marker emit が発散原因」は**誤り**。single-owner baseline が同一に発散し、
> 真因は fresh-fine conical.msh 自体の不安定 (別問題)。コーナー修正の純粋検証には「コーナーを持ち
> かつ node viscous で安定なメッシュ」が必要で未確保。詳細は plan §8.1.2.x。

#### メッシュ均一化 + B 確定 + 自走検証 (2026-06-20, plan §8.1.2.3.1, §9)

**メッシュ修正**: `make_nozzle.py` の `NX_ZONE['throat_dn']` 44→10 + `--prog-x` 既定 1.0→1.012 で
スロート下流の軸方向 dx 段差 (19倍) を解消 (vol min/mean 3.7e-4→4.8e-2)。`mesh/conical.msh`/`bell.msh`
再生成 (各 24624 ノード, `.msh` は gitignore 生成物で make_nozzle.py がソース)。

| run_* | 設定 | 結果 |
| --- | --- | --- |
| `run_0031_corner_single_dir` | single-owner + Dirichlet (均一メッシュ, 30k) | 解正常だが出口リップで **max cfl=2.92e25 (退化)** |
| `run_0032_corner_multi_dir` | **multi-marker + Dirichlet (均一, 100k)** | **cfl 0.1・収束 PASS**。中心線 Mach=SU2-lam ~2-5%, 推力 1889N, no-slip BL 正常 |
| `run_0033_wallclust_bump_node` | **Bump0.4 壁寄せ + node visc + Dirichlet (10k)** | 安定・**Tmax 1508 (Tt 1500 +0.6%)** = T overshoot ほぼ解消, 推力 1907N |
| `run_0034_cfl{0p2..2p0}_wc` | 壁寄せ CFL 掃引 (explicit RK3) | 0.7 安定 / 0.8 marginal / **1.0 発散** → explicit 上限 ≈0.7 |

> **結論**: B (multi-marker + nodeWallDirichlet) を node 壁の主軸に採用。SU2-lam と一致・no-slip BL 正常・
> 回帰なし (case08/24)。**T overshoot は近壁解像不足アーチファクトで、壁寄せ (Bump0.4) でほぼ消える**
> (均一 +13% → 壁寄せ +0.6%)。explicit CFL 上限 ~0.7。詳細は plan §9。

#### node RANS (SST) を働く壁レシピで実走 + GG/LSQ/Kato 比較 (2026-06-21)

これまで node は層流のみで **RANS (SST) は設定だけで未計算** (`run_fx_sst_node`/`run_divufix_sst_node_*` は residual 無し=未実走、しかも nodeWallDirichlet 未設定)。働く壁レシピ (run_0033: `nodeWallDirichlet=1`+壁寄せ+axisCentroidShift) に SST を載せ、converged cell SST (`run_0017`) から cross-mesh restart して実走した。

| run_* | 設定 (共通: node, SLAU, RK3, cfl0.1, convMethod0, nodeWallDirichlet=1, conical 壁寄せ) | VERDICT (20000 step) | 近軸 roK max / T max |
| --- | --- | --- | --- |
| `run_0038_node_sst_wc` | GG (gradLSQ=0) + **katoLaunder=1** | **NOT CONVERGED (still converging)** 全列 falling・NaN無し (rms_ro 3.6dec, roOmega 6.8dec, roUy 2.9dec) | 1.12e5 / 1500 |
| `run_0039_node_sst_nokato` | GG + **katoLaunder=0** (cell run_0017 と条件一致) | **同上** (rms_ro fin 2.19e-8、run_0038 とほぼ同一) | **1.64e5** / 1500 |
| `run_0040_node_sst_lsq` | **LSQ (gradLSQ=1)** + katoLaunder=0 | **同上** (NaN無し、roUy 3.0dec) | 1.08e5 / 1500 |

> **結論**: **node RANS SST は働く壁レシピで安定収束に向かう** (3本とも NaN無し・全残差 falling・T≤Tt overshoot無し)。
> いずれも 20000 step では未収束 (要追加 step) だが健全。
> - **katoLaunder**: 収束・バルク場は不変、近軸 k スパイクのみ ~30% 縮小 (kato0 1.64e5 → kato1 1.12e5)。安定性には不要。
> - **LSQ (gradLSQ=1)**: **nodeWallDirichlet 併用なら発散しない** (旧 `run_lsq_on` の step24716 発散は壁レシピ無しが原因)。近軸 k 最小・roUy 収束わずかに良。
> 比較図 `node_sst_gg_vs_lsq_kato.png`、各 `run_*/residual_history.png`、各 run の `CONVERGENCE_VERDICT.txt`。

### 近軸 混合精度・閉形式 FVS の調査と根治 (2026-06-13〜14, plan [precision-mixed-axisym.md](../../plans/archived/precision-mixed-axisym.md))

laminar conical 第一セル `Uy`(x=40mm) が物理値 −15 でなく float 陰解で −0.6 固着する問題の切り分け〜本番化。
**結論: 真因は block-DPLUR 線形 solve の精度。閉形式 FVS を既定化 (float ~10%速)、`implicitSolvePrecision=1` で根治。**

| run_* | 目的・設定 | 主要結果 | 状態 |
| --- | --- | --- | --- |
| `run_axis_lam_slau` | 参照: 純 float 陰解 (laminar conical) | Uy0=−0.645 (固着) | ref |
| `run_axis_lam_slau_double` | 参照: global double | Uy0=−15.1 (正) | ref |
| `run_axis_lam_slau_expl` | 参照: explicit float | Uy0=−14.9 (正) | ref |
| `run_su2cmp_forge_sst` | 参照: SST conical (float, k スパイク) | k_axis=9.16 | ref |
| `run_0001_mixed_lam_slau` | double 残差+**float** solve (itref 単純分割) | Uy0=−0.641 **固着→棄却** | active |
| `run_0002_dsolve_fres_lam_slau` | **float 残差+double solve** (対称実験) | Uy0=−14.88 ★真因=solve精度を確定★ | active |
| `run_0003_dsolve_sst` | SST + double solve | k_axis 9.16→7.4 (縮小・残留) | active |
| `run_0004_float_timing` | 純 float baseline 速度 | 91.5s (×1.0) | active |
| `run_0005_dsolve_bs64` | double solve, block size 64 | 502.9s (×5.5, bs 効果薄) | active |
| `run_0006_dsolve_restruct` | double solve, R/L 排除再構成 | 426.9s (×4.67) | active |
| `run_0007_dsolve_closedform` | **double solve, 閉形式 A±** | 240s (×2.62), Uy0=−15.03 | active |
| `run_0008_dsolve_cf_sst` | SST 閉形式 double solve | 249s, k_axis=7.45 (run_0003 と一致) | active |
| `run_0009_dsolve_shareddiag` | 閉形式 double + diag を shared | 239s (スピル除去も速度不変) | active |
| `run_0010_dsolve_df` / `run_0011_dsolve_df_fix` | compensated double-float (部分) | 固着 (精度不足) | active |
| `run_0012_dffull` | full double-float | 固着 (~48bit 不足→df 棄却) | active |
| `run_0013_float_cf` | **float 閉形式** (R/L 排除) 速度 | 82.8s (×0.90 = 10%速), 数値 legacy 等価 | active |
| `run_0014_prod_float` | **本番** `implicitSolvePrecision=0` (既定 float) | 82.6s, Uy0=−0.6449 | active |
| `run_0015_prod_double` | **本番** `implicitSolvePrecision=1` (double solve) | 234.5s, Uy0=**−14.98** (根治) | active |

### 近軸固着の切り分け: 粘性 vs convMethod (2026-06-14)

「近軸 Uy 固着の原因は粘性か 2 次精度 (convMethod) か」を分離する 2×2 マトリクス。**全 run cell・float
(`implicitSolvePrecision=0`)・conical・cfl=2・blockDPLUR=1・8000 step・同一等エントロピー IC**、変えるのは
物理 (`viscMethod`/`visc`) と `convMethod` のみ。第一セル Uy(x≈27.8mm) で固着を判定 (物理値 ≈ +13)。

| run_* | 物理 | convMethod | 収束 (check_convergence) | 第一セル Uy | 固着? |
| --- | --- | --- | --- | --- | --- |
| `run_disent_eul_cm0` | Euler | 0 (1次) | PASS (rms_roUy 4.6dec) | **+13.6** | なし |
| `run_disent_eul_cm1` | Euler | 1 (2次) | NOT CONV (rms 2.1dec 停滞) | **+17.9** | なし |
| `run_disent_lam_cm0` | laminar | 0 (1次) | NOT CONV (rms 1.8dec 停滞) | **+0.86** | **固着** |
| `run_disent_lam_cm1` | laminar | 1 (2次) | NOT CONV (rms 1.4dec 停滞) | **+1.41** | **固着** |

**結論**: 近軸 Uy 固着の原因は **粘性 (`viscMethod=1`)** であり convMethod ではない。laminar は cm0/cm1 とも
固着 (第一セル ~+1、第二セルで +110 へ暴騰=偽勾配)、Euler は cm0/cm1 とも健全 (軸から滑らかに立ち上がる)。
convMethod=1 (2次) は固着とは無関係の**良性の残差プラトー** (リミタのリミットサイクル) を生むだけで、
Euler cm1 は残差が止まっても場 (Uy=+17.9) は正常。plan [architecture-axisym-axis-singularity.md] の真因
(粘性が近軸 r→0 で剛性を増し float block-DPLUR の `D⁻¹` が悪条件化、陽解/Euler は `D⁻¹` 不使用で正常) と整合。

### 粘性対角の幾何是正で float のまま固着解消 (2026-06-14, 上記の真因の根治)

上の切り分けで「固着=粘性の LHS 由来」と判明したので、block-DPLUR の**粘性対角の形**を精査したところ、
`face_area·(2ν/delta)` が `delta=dcc·ss²/dcc_dot_s` (面積) を長さ²扱いして **≈2ν に潰れ**、軸対称近軸で
r 重みを失い・ゼロ面積(対称)面にもスプリアス項を載せていた (residual 側 [viscousFlux_d.cu] は
`ip<nNormalPlanes` で軸面を除外)。これを residual の粘性流束 Jacobian と整合する **`2ν·ss²/dcc_dot_s`
(= `2ν·delta/dcc`、∝r でゼロ面積面では消える)** に是正 (`timeIntegration_d.cu` の scalar/block/precond 3 箇所)。

検証 (基準は `run_disent_lam_cm1`、build-viscfix=是正 float / build-verify=未修正):

| run_* | ビルド・精度 | 収束 | 第一セル Uy | 判定 |
| --- | --- | --- | --- | --- |
| `run_disent_lam_cm1` | float 未修正 | NOT CONV (リミタ) | **+1.4** (固着) | baseline |
| `run_disent_lam_cm1_viscfix` | **float 是正** | NOT CONV (リミタのみ) | **+17.9** | **固着解消** |
| `run_disent_lam_cm1_double` | float 未修正 + `implicitSolvePrecision=1` | NOT CONV (リミタ) | +18.1 | 参照(double根治) |
| `run_disent_lam_cm0` | float 未修正 (1次) | NOT CONV (停滞) | +0.86 (固着) | baseline |
| `run_disent_lam_cm0_viscfix` | **float 是正** (1次) | **PASS (rms_roUy 4.5dec)** | **+12.1** | **固着解消＋収束回復** |

**結論 (旧結論の更新)**: 近軸固着の本質は **float 精度ではなく LHS 粘性対角の幾何不整合**だった。是正後は
**float のまま** double solve (`run_0015_prod_double`) と一致 (+17.9 vs +18.1、cm1) し、1次では NOT CONV→PASS に回復
(cm0)。`implicitSolvePrecision=1` は悪条件 LHS を倍精度で押し切る対症療法で、LHS を直せば**追加 FP64 不要・float 速度**で
根治できる。plan [precision-mixed-axisym.md] / [architecture-axisym-axis-singularity.md] の「double solve が必要」は
本是正で更新対象 (要 methods/plan 反映・他ケース回帰)。

#### SST (軸対称 RANS) でも軸中心 k スパイクが float のまま消える (本丸の確認)

本来の調査対象だった軸対称 SST の「軸中心 k スパイク」を、同一 IC (`run_su2cmp_forge_sst` の `nozzle.h5`)・
同一設定 (cfl_pseudo=2, 20000 step) で baseline(float 未修正) / viscfix(float) / double solve の3者比較:

| run_* | 第一セル Uy | k_axis | k_axis/k_core | スパイク |
| --- | --- | --- | --- | --- |
| `run_su2cmp_forge_sst` (baseline float) | +1.06 (固着) | 7.76 | 1.7 | あり |
| `run_disent_sst_viscfix` (**float 是正**) | **+19.30** | **4.65** | **1.0** | **消失** |
| `run_disent_sst_double` (`implicitSolvePrecision=1`) | +19.85 | 4.19 | 0.9 | 消失 |

viscfix(float) は Uy 固着解消・k スパイク消失 (k(r) 平坦) で **double solve と一致**。乱流まで含めて float のまま根治。
注: 3 run とも完全収束ではない (`rms_roUy` rising、`rms_roK/roOmega` は残差ノルムが NaN になる既知アーティファクト=
場は finite。baseline/double も同挙動)。よって「収束」とは報告せず、同一土俵での k スパイク/Uy 診断の比較として扱う。
plan の「SST k スパイクは double solve でも縮小・残留 (9.16→7.4)」は、本是正でも近軸では同等に解消することを示す。

### 粘性対角是正後の scalar DPLUR CFL スイープ (2026-06-14)

粘性流束 Jacobian 是正後、これまで全く収束しなかった **scalar DPLUR (`blockDPLUR: 0`)** の安定性が
変わったかを検証。基準 `run_disent_lam_cm1_viscfix` (block DPLUR, cfl 2.0, viscfix float build) の設定を
そのまま使い `blockDPLUR: 0` のみ変えて CFL をスイープ (全 run: SLAU, 軸対称, laminar, conical メッシュ +
等エントロピー IC, `timeIntegration: 11`, `nStepInner: 20`, 2 次 + limiter 2, 8000 step, build-viscfix)。

| run_* | CFL | NaN/発散 | 収束 (check_convergence) |
| --- | --- | --- | --- |
| `run_disent_scalar_cfl0p1` | 0.1 | NaN なし | NOT CONV (rms_ro/roUx **上昇**, roUy 停滞) |
| `run_disent_scalar_cfl0p2` | 0.2 | NaN なし | NOT CONV (**全列上昇** — 緩慢発散, roUy +2.9 roe +5460) |
| `run_disent_scalar_cfl0p5` | 0.5 | **DIVERGED** ~step 557 | NaN |
| `run_disent_scalar_cfl1`   | 1.0 | **DIVERGED** ~step 246 | NaN |
| `run_disent_scalar_cfl2`   | 2.0 | **DIVERGED** ~step 85  | NaN |
| `run_disent_scalar_cfl5`   | 5.0 | **DIVERGED** ~step 35  | NaN |
| `run_disent_scalar_cfl10`  | 10.0| **DIVERGED** ~step 32  | NaN |
| `run_disent_lam_cm1_viscfix` (基準) | 2.0 (**block**) | NaN なし | NOT CONV だが安定 (rms_ro 1.3dec/roUx 2.0dec 後プラトー) |

**結論**: 粘性対角の是正をもってしても **scalar DPLUR は依然として収束しない**。CFL≤0.2 で辛うじて NaN を
免れるが残差は下げ止まり〜上昇 (収束せず)、CFL≥0.5 で完全発散し、発散開始 step は CFL とともに単調に早まる
(古典的な CFL 安定限界)。一方 **block DPLUR は同 CFL 2.0 で安定**に回り 1〜2 dec 下げてプラトーする。
よって陰解法の実用パスは引き続き block DPLUR であり、粘性 Jacobian 是正は block の近軸固着を解消したが、
scalar DPLUR の対角近似 (運動量連成を落とす) の弱さを救うものではない。残差図は各 run の `residual_history.png`。

#### 切り分け: 陽解法は安定なのに scalar DPLUR が発散する真因 (2026-06-14)

「陽解法が回るのに陰解法 (scalar) が発散するのは式間違いでは」という疑義を、同一 conical メッシュ +
等エントロピー IC で対照実験して切り分けた。

| run_* | スキーム | CFL | 内反復 | 結果 |
| --- | --- | --- | --- | --- |
| `run_disent_expl_cfl0p3` | explicit RK3 | 0.3 | – | **安定** (rms_ro −1.2dec) |
| `run_disent_expl_cfl0p5` | explicit RK3 | 0.5 | – | **安定** (rms_ro −1.2dec) |
| `run_disent_scalar_cfl0p5_inner1` | scalar DPLUR | 0.5 | **1** | NaN なしだが残差**上昇** |
| `run_disent_scalar_cfl0p5` | scalar DPLUR | 0.5 | 20 | **NaN @ step 557** |
| `run_disent_lam_cm1_viscfix` | block DPLUR | 2.0 | 20 | 安定 (プラトー) |

確定事項:
- **陽解法 (RK3) は CFL 0.5 で安定**。よって発散は物理/メッシュでなく **scalar DPLUR の LHS** に起因。
- **内反復を増やすほど悪化** (1 sweep=上昇のみ, 20 sweeps=NaN)。収束した scalar-LHS の解そのものが不安定化方向
  = LHS が真の Jacobian の不十分/不整合な近似。

##### 仮説1 (軸対称ソース Jacobian 欠落) — 実装・検証したが**棄却**

当初「scalar カーネルが軸対称フープ源 (`res_roUy += (P−τθθ)·A_planar`) の Jacobian を持たない
(block は `diag_block[2][2]` に陰化済) のが真因」と考え、scalar カーネルに `isAxisymmetric`/`A_planar`/
per-cell `gamma` を渡し roUy 対角へ同形項 `A_pl·((γ−1)u_y + 2μ/(ρ r_eff))` を陰化した
([timeIntegration_d.cu] `implicit_defect_correction_d`, TP=per-cell γ 対応)。

結果は **before/after でほぼ不変** (`run_disent_scalarfix_*`): cfl0.5 NaN@557→557, cfl1 @246→246,
cfl2 @85→82。軌跡は摂動するが発散段は同じ。よって**軸対称ソース項は本ケース発散の主因ではない**。
(この Jacobian 追加自体は block との整合・他軸対称ケース向けに正しい改善なので残置。)

##### 真因 — 出口コーナーでの 2 次 MUSCL オーバーシュート (空間的に切り分け)

`run_disent_scalarfix_cfl2_diag` (dense 出力) で発散セルの座標を追うと、**軸ではなく出口面・外壁コーナー**
(x≈0.085=出口, y≈0.032≈R_e、`CELLS/centCoords` で確認) で先に壊れる: P が下限 1 Pa に貼り付き
(過膨張)→ ρ が負 → Ux が大きな負 (逆流) に発散。剛フープ源の軸近傍ではない。

決定打は **1 次精度 (`convMethod: 0`) で scalar DPLUR が安定収束**すること:

| run_* | 空間 | CFL | 結果 |
| --- | --- | --- | --- |
| `run_disent_scalarfix_cfl0p5` | 2次 MUSCL | 0.5 | NaN@557 |
| `run_disent_scalarfix_cfl2`   | 2次 MUSCL | 2.0 | NaN@82 |
| `run_disent_scalarfix_o1_cfl0p5` | **1次** | 0.5 | **安定 −2.4dec** |
| `run_disent_scalarfix_o1_cfl2`   | **1次** | 2.0 | **安定 −1.6dec** |

→ scalar DPLUR の**式は構造的に正しい** (1次なら cfl2 でも収束)。発散は **2次 MUSCL が出口コーナーの強膨張で
オーバーシュート**し、それを scalar の**弱い (式分離・1次凍結) LHS が減衰できない**ことに起因する
(defect-correction で LHS=1次・residual=2次の不整合)。block は真の 5×5 Jacobian の式間連成で同 2次残差を
cfl≈5 まで減衰でき、scalar のスペクトル半径対角はできない。陽解法 RK3 は広い安定域で 2次残差を直接捌けるため
安定。**式の符号ミスではなく、scalar LHS の前処理能力不足**が結論。

**乱流 (SST)**: 別経路 (`applySSTPointImplicit`, `D_φ=V/Δτ+V·src_jac+transport_diag`) で mean-flow の
scalar/block 選択に依存しない。RANS でも mean-flow を 2次 scalar DPLUR で回せば同じ出口コーナーで発散する
(SST 側に固有のバグはない)。

**実用指針**: 2次精度の陰解法は **block DPLUR** を使う。scalar DPLUR を使うなら **1次** (起動・ロバスト用)
に限る。残差図は各 run の `residual_history.png`。

### 軸対称 viscous 発散項の整合化 (divu fix) の before/after (2026-06-14, plan [diffusion-viscous-shear-flux.md](../../plans/accepted/diffusion-viscous-shear-flux.md))

planar 面の体積粘性 `-2/3 μ divu` がデカルト発散のみで**フープ項 $u_r/r$ を欠落**し、完全発散 `axisym_divU` を使う
$\tau_{\theta\theta}$ 源項と不整合だった件の修正 ([viscousFlux_d.cu])。**BEFORE=build-viscfix (修正前) / AFTER=build-lsq (修正後)**、
差分は divu のみ。非軸対称は `isAxisymmetric==1` ガード下でビット不変。

| run_* | mode | ビルド/精度 | 収束 (check_convergence) | divu fix 効果 (vs BEFORE) | 状態 |
| --- | --- | --- | --- | --- | --- |
| `run_lsq_gg` | node laminar expl RK3 40k | build-lsq 修正前 (gradLSQ=0) | NOT CONV (rms_ro 5e-8, roUy/roe 低下中) | — (= node BEFORE) | active |
| `run_divufix_node_after` | node laminar expl RK3 40k | **build-lsq 修正後** | 同上 (rms_ro 5e-8) | **P max 0.001% / 全量 <0.006% / 平均 <0.0003%** | active |
| `run_divufix_cell_before` | cell laminar impl cfl2 8k | build-viscfix 修正前 | NOT CONV (plateau, 1.3–2.1dec) | baseline | active |
| `run_divufix_cell_after`  | cell laminar impl cfl2 8k | **build-lsq 修正後** | NOT CONV (plateau) | P max 1.6% / Uy 3.5% (喉部局在) ※**未収束プラトー同士で過大** | active |

**結論**: divu 修正は**正しい整合改善**だが、**層流ノズルの収束場では効果は実質ゼロ** (node 収束比較で P 0.001%)。
理由は $\tfrac23\mu_{lam}(u_r/r)$ が圧力・対流項に比べ桁違いに小さいため。cell の 1.6% は**未収束プラトー同士の比較**で
過大評価 (両 run とも NOT CONV)。**乱流 (μ_t ≫ μ_lam) や強圧縮の収束場でこそ効きうる**ため、実害確認は SST 収束ケースの
before/after が本命 (未実施)。残差図は各 run の `residual_history.png`。

> 関連: 同セッションで LSQ 勾配 (`gradLSQ`, plan [discretization-lsq-gradient.md](../../plans/active/discretization-lsq-gradient.md)) を
> node viscous で試行 (`run_lsq_on`/`run_lsq_diag`) → **近壁 M 退化で発散・棄却** (詳細は plan §9)。`run_lsq_gg` は
> その GG ベースライン (gradLSQ=0) で、本 divu 比較の node BEFORE を兼ねる。

#### SST (乱流) での divu fix 効果

| run_* | mode | ビルド | divu fix 効果 (vs before) |
| --- | --- | --- | --- |
| `run_divufix_sst_cell_before` / `_after` | cell SST (impl 20k, **未収束** roe~42) | viscfix / **lsq** | res_20000 (μ_t max 1.05e-3≫μ_lam): **P 1.1% / Uy 4.4% / k 2.6% / vis_turb 1.5%** (喉部局在, 平均<0.02%) |
| `run_divufix_sst_node_before` | node SST (expl) | viscfix | **step~1300 で NaN 発散** → node SST は検証不可 |

**SST 所見**: cell SST では μ_t 発達後に divu fix が乱流量 (k/vis_turb) を ~1.5–2.6% 動かす (層流より大きいが依然局所・未収束)。
膨張ノズル核は加速流で乱流が弱く (核 μ_t 小)、フープ膨張粘性の効きは限定的。**node SST は node-viscous の既知脆弱性
(§7.2)+SST の ω 剛性で baseline build (修正前) から発散**するため、divu/fx いずれの検証も不可。

### node 面補間を中点 fx=0.5 に固定 (median-dual 標準化, 2026-06-14, plan [discretization-median-dual.md](../../plans/active/discretization-median-dual.md) §9)

node 内部双対面の面補間係数を `dualFaceCent` 射影比でなく **fx=0.5 (ノード間中点 $\phi_f=½(\phi_A+\phi_B)$)** に固定
([calcStructualVariables_d.cu], cell 不変)。安定な node 層流 (40k) で `run_divufix_node_after` (幾何 fx) vs `run_fx_node` (fx=0.5)、
差は fx のみ (両者 divu fix 込み, build-lsq vs build-lsqfx)。

| 指標 | 幾何 fx | fx=0.5 |
| --- | --- | --- |
| 近壁 dUxdy roughness (median/90/99pct) | 0.98 / 4.75 / **12.84** | 0.78 / 3.81 / **8.23** (−36%) |
| 場の差 (vs 幾何fx) | — | P 0.14% / ro 0.76% max・平均<0.26% |
| 局所最大 | — | Ux 出口近傍 near-wall で 1376→430 (敏感域, §7.2 exit-lip) |

**所見**: fx=0.5 は**近壁 checkerboard を全体に低減** (期待どおり標準 median-dual) する一方、**出口リップ近傍 near-wall の解を
実質的に変更**する (局所 Ux 1376→430)。

**SU2 クロスチェック (壁圧, `run_su2cmp_su2_lam`, axisym laminar 同条件)**: **fx ON/OFF は壁圧で区別不能** (<0.5% 差)、
両者とも SU2 に平均 **3.2%** 一致 (最大 +15% は未収束の超音速出口 x=85mm)。→ **fx=0.5 は安全** (SU2 一致を悪化させず
checkerboard を下げる)。局所 Ux 差は圧力場に伝播しない近傍速度の細部で、**SU2/forge とも未収束のため near-wall 速度の
厳密な是非は未確定**。**フラグ `nodeMidpointFx` (既定 0=OFF) で opt-in 実装** (cell 無関係)。

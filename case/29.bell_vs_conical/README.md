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
| `run_0008_bell_rans` / `run_0009_conical_rans` → 継続 `run_0014`/`run_0015` | RANS SST (乱流 BL 発達済) | no-slip | Sutherland+SST |

> いずれも壁密格子 (`--prog-r 0.95 --ny 100`, 第1セル ~10µm, y+≈0.6)。
>
> **RANS の乱流着火 (turbulent ignition) は遅い**ので注意。$k$/$\omega$ 拡散は有効粘性
> $\mu_{lam}+\sigma\,\mu_t$ ([scalarTransport_d.cu](../../solver_density_cuda/cuda_forge/scalarTransport_d.cu)) で
> 正しく実装されているが、層流的 IC からは $\mu_t\approx0$ のため初期は**層流拡散のみ**で $k$ の
> 立ち上がりが極端に遅い (自己加速するまで潜伏)。順圧勾配 (喉の再層流化) と低い freestream seed
> がこれを更に長引かせる:
> - 12k step: 近壁 $k\approx0$, $\mu_t/\mu_{lam}\sim1.3$ — **未発達** (層流とほぼ同じ)。
> - 40k step: x≳40mm のみ着火 ($k\sim10^4$)、x≤25mm はまだ $k\sim10^2$ — **部分発達**。
> - ~50–60k step (warm-start 継続 `run_0012→0014`): 全 x で着火し $\mu_t/\mu_{lam}$ ピーク
>   ~25 (BL 端) の正常な乱流 BL に収束 (`rans_visturb_profile.png`)。
>
> 発達後の RANS 推力は未発達 (40k) 比 −0.2% (ベル) と僅かに低い (摩擦が現実的に増える)。本表は
> **発達済** (`run_0014`/`run_0015`) の値。

## 結果 — 推力比較 (同じ喉 / 軸長 / 出口面積)

NaN/発散チェック: 全 run NaN なし、最終場は物理的 ($P>0,\ T>0,\ P_{max}\approx P_t$)。
**mdot はベル/コーンで一致 (≤0.2%)** = 同一チョーク喉の妥当性 (出口面積 $A_e=\pi R_e^2=3367$ mm²)。

| モデル | ベル $F$ [N] | コーン $F$ [N] | **ベル推力ゲイン** | $\lambda$ ベル/コーン | 出口流れ角 ベル/コーン |
| --- | --- | --- | --- | --- | --- |
| inviscid (slip) | 1976.2 | 1962.9 | **+0.68%** | 0.996 / 0.983 | 3.7° / 10.0° |
| laminar | 1978.9 | 1957.2 | **+1.11%** | 0.997 / 0.984 | 3.5° / 9.8° |
| RANS SST (発達済) | 1972.8 | 1953.4 | **+0.99%** | 0.997 / 0.984 | 3.4° / 9.7° |

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
- cf. [[implicit-blockdplur-config]] / [[nozzle-implicit-cfl-ceiling]]。

主要成果物: `mach_comparison.png`, `exit_profiles.png`, `rans_visturb_profile.png`,
各 `run_*/residual_history.png`, `mesh/nozzle_preview.png`。

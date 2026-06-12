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
[`architecture-axisym-axis-singularity.md`](../../.github/plans/architecture-axisym-axis-singularity.md)。

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

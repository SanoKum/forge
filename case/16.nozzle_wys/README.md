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

> **⚠️ 重要 (2026-06): `mesh/nozzle_H.geo` は Fig.3 の実験ノズルと別形状**だった
> (`13.nozzle_H` とバイト同一の別 Wyslouzil 論文ノズル。スロート半 2.25mm・A/A\*=1.69 で過発散)。
> JCP 113,7317 (2000) Fig.3 の正しいノズルは **スロート全高 5mm・収束 38mm・発散 95mm 直線壁・
> 壁間 1.8°・A/A\*≈1.58** で `mesh/make_nozzle_fig3.py` が生成する (`nozzle_fig3_2d/3d`)。
> **凝縮検証は run_0047 以降 (修正形状) を使うこと**。run_0001〜0046 (nozzle_H) は旧形状で、
> dry が実験より過膨張・凝縮も相殺で「合って見えていた」だけ。詳細は下記「## 修正形状 (Fig.3) run」。

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

## 非平衡凝縮 (H2O) run 一覧 — Wyslouzil Fig.3 検証

Wyslouzil et al. JCP 113, 7317 (2000) **Fig.3 (pv0=1.0 kPa 水)** 条件で、N2 キャリア中の
希薄水蒸気 (Y_H2O=0.01095) の非平衡凝縮を 4 モーメント法で計算し、凝縮あり/なしの静圧比を比較する。
凝縮潜熱で中心線静圧が dry 等エントロピー線より上振れする現象 ([condensation plan](../../.github/plans/condensation-nonequilibrium.md))。

- 入口: `inlet_Pressure` Pt=59070 Pa, Tt=286.65 K, Y(N2,H2O)=(0.98905, 0.01095)。
- **carrier=CPG**: 低温 (<200K) で NASA-9 TP が外挿不可・発散するため、N2 キャリアは CPG (γ=1.4) で解く
  (この温度域の N2 cp はほぼ一定で CPG が正確かつ頑健)。H2O は移流種 (vapor budget) として追跡し、
  凝縮潜熱・物性 (Murphy-Koop psat / Hertz-Knudsen 成長) は CPG エネルギー経路に組み込む。
- 検証図: `fig3_compare.png` (非粘性), `fig3_compare_visc.png` (非粘性+粘性SST 重ね合わせ)。

| run | 物理 | 凝縮 | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0007_h2o_cond` | 2D 非粘性 (CPG+species) | off | dry baseline。収束 (rms_ro↓3.5dec) | active |
| `run_0008_h2o_cond_on` | 2D 非粘性 | **on** | 潜熱で T 143→185K, p/p0 +19〜28%。cond/dry 比が exp と ~5% 一致 | active |
| `run_0009_h2o_sst_dry` | 2D 粘性+SST 乱流 | off | dry SST baseline (Wyslouzil 条件)。T≤Tt | active |
| `run_0010_h2o_sst_cond` | 2D 粘性+SST 乱流 | **on** | 凝縮効果 (ユーザ要望)。等温CNT + Hertz-Knudsen (baseline) | active |

**核生成/成長モデル感度 (SST, 凝縮 on, `condKantrowitz`×`condGrowthModel`)** — `compare_models.png`:

| run | Kantrowitz | 成長則 | onset 位置 | 備考 |
| --- | --- | --- | --- | --- |
| `run_0010_h2o_sst_cond` | off | Hertz-Knudsen | x≈1.5cm | baseline。等温CNT は onset 早すぎ・bump 過大 |
| `run_0011_h2o_sst_kw_hk` | **on** | Hertz-Knudsen | x≈2.3cm | Kantrowitz が J 抑制 → onset 遅延 (実験寄り)、bump 控えめ |
| `run_0012_h2o_sst_nokw_gyar` | off | **Gyarmathy** | x≈1.3cm | Gyarmathy は成長急 → bump が最も鋭く過大 |
| `run_0013_h2o_sst_kw_gyar` | **on** | **Gyarmathy** | x≈2.0cm | **onset 位置・peak が実験に最も一致** (Kantrowitz 遅延 + Gyarmathy 急成長) |

知見: **Kantrowitz 非等温補正は核生成 J を抑え onset を下流へ遅らせる**(等温CNT は onset が早すぎ実験の
bump 位置を外す)。**Gyarmathy(熱伝導律速)は Hertz-Knudsen より成長が急**で p/p0 の立ち上がりが鋭い。
両者を併せた `run_0013` (Kantrowitz+Gyarmathy) が x≈2〜3.5cm の onset/peak を実験と最もよく再現
(x=3.1cm: forge 0.353 vs exp 0.354)。下流 x≳4cm は全モデルほぼ収束し、実験より ~0.02 低い (膨張側オフセット)。

**Gyarmathy Knudsen 係数 C 感度 (`condGyarmathyC`, base=Kantrowitz+Gyarmathy)** — `compare_gyarC.png`:

| run | C | onset | x=3.1cm p/p0 |
| --- | --- | --- | --- |
| `run_0014_h2o_sst_gyarC0p0` | 0 (連続体, Kn 補正なし) | 最も早い (x≈1.3) | 0.342 |
| `run_0015_h2o_sst_gyarC1p59` | 1.59 | 早い | 0.348 |
| `run_0013_h2o_sst_kw_gyar` | **3.18 (標準)** | x≈2.0 | **0.353** (exp 0.354) |
| `run_0016_h2o_sst_gyarC6p36` | 6.36 | 遅い | 0.328 |
| `run_0017_h2o_sst_gyarC12p72` | 12.72 | 最も遅い (x≈3.5) | 0.284 |

C は成長率の Knudsen 抑制を制御し、**小さいほど成長が速く onset が早まり、大きいほど遅く弱くなる**
(g プロファイルの立ち上がり位置が C で単調にシフト)。**標準値 3.18 が onset/peak を実験に最も良く再現**し、
ずらすと悪化 (C<3.18 で早すぎ・C>3.18 で遅すぎ)。下流 x≳5cm は最終的に全量凝縮するため C 依存は小。
→ **この系では Gyarmathy 標準係数 3.18 がほぼ最適**であることを確認。

**液滴温度 T_d 考慮 (一温度 vs 二温度, `condTwoTemp`, Hertz-Knudsen)** — `compare_2temp.png`:

| run | モデル | 液滴温度 |
| --- | --- | --- |
| `run_0010_h2o_sst_cond` | isoCNT+HK | 一温度 (T_d=T_g) |
| `run_0018_h2o_sst_nokw_hk_2T` | isoCNT+HK | **二温度 (Hill T_d)** |
| `run_0011_h2o_sst_kw_hk` | Kantrowitz+HK | 一温度 |
| `run_0019_h2o_sst_kw_hk_2T` | Kantrowitz+HK | **二温度** |

知見: **この希薄水/N2 系では液滴温度の影響は小さい**(p/p0 差 ≤0.008、概ね ~0.002 = <1%)。液滴超過温度
T_d−T_g は**成長活発な onset 前線で局所的に最大 ~15K**(過冷却 ~25–30K の一部)に達するが、中央値は ~0.2K
(下流は凝縮完了で蒸気枯渇→成長停止→T_d→T_g)。**N2 キャリアが潜熱を効率よく奪う**ため巨視的影響は
Kantrowitz/Gyarmathy のモデル選択より一桁小さい。方向は予想どおり成長を僅かに遅らせ onset を微小に下流へ
(g プロファイルで二温度=破線が僅かに遅延)。→ **純蒸気凝縮では T_d が必須だが、キャリア希薄凝縮では一温度近似で十分**、という結論。実装詳細は [docs/condensation/theory.md §4](../../docs/condensation/theory.md)。

**複数 H2O 分圧スイープ (モデル汎化, Kantrowitz+Gyarmathy 固定)** — `compare_multicond.py` / `compare_multicond.png`:

入口 H2O 質量分率のみを変えて pv0 を 1.00 / 0.50 / 0.26 kPa とし、Wyslouzil Fig.3 multi-P (`wyslouzil_fig3_multiP.csv`, P0=59.07 kPa) と中心線 p/p0 を比較。**単一条件で較正したモデルを他分圧に汎化できるかの検証**。

| run | pv0 | 入口 Y_H2O | onset (exp) | x=3.1cm p/p0 (forge / exp) |
| --- | --- | --- | --- | --- |
| `run_0013_h2o_sst_kw_gyar` | 1.00 kPa | 0.01095 | x≈2cm | 0.353 / 0.354 |
| `run_0020_h2o_sst_kg_0p50kPa` | 0.50 kPa | 0.005461 | x≈3.5cm | 0.279 / 0.297 |
| `run_0021_h2o_sst_kg_0p26kPa` | 0.26 kPa | 0.002835 | x≈5cm | 0.274 / 0.293 |

知見: **分圧依存 (凝縮 onset の下流シフトと圧力上昇幅の序列) を 3 条件すべてで再現**。pv0 が下がるほど onset
が下流へ移り bump が弱まり dry に漸近する、という実験傾向を捉える。forge は絶対 p/p0 が一様に exp より
~5–8% 低い (既知の dry baseline オフセット = forge dry が非粘性等エントロピー、実験は粘性排除厚で上振れ) が、
**凝縮物理の相対挙動・分圧応答は単一条件較正なしで汎化**する。
(注: `run_0021` は step6000 で外部中断したが、p/p0・g とも step2000→6000 で 4 桁不変＝定常化済みのため発達場として採用。0.26 kPa は g max≈6e-4 と凝縮ごく僅か。)

**非平衡凝縮の検証結果**:

- **非粘性 (run_0008)**: 凝縮ありの中心線 p/p0 は dry 比で 1.19〜1.28 倍 (スロート下流 2.5〜8 cm)、
  実験の cond/isentrope 比 1.18〜1.23 と ~5% 以内で一致。**凝縮効果 (比) は定量一致**するが、絶対 p/p0 は
  非粘性 dry が実験等エントロピーよりやや低い (粘性排除厚さ未考慮)。
- **粘性+SST 乱流 (run_0010 vs run_0009, ユーザ要望)**: 粘性排除厚さで dry baseline が押し上がり、
  **絶対 p/p0 も実験と一致**する (下表)。凝縮 onset (x≈1.7cm の潜熱スパイク) も再現。x≳6cm の下振れは
  壁解像不足ではなく (y+ は良好、下記)、2D 近似/ノズル面積比など膨張側のオフセット (dry/cond 共通)。

  | x [cm] | SST dry | SST cond | exp isentrope | exp cond |
  | --- | --- | --- | --- | --- |
  | 2.5 | 0.300 | 0.381 | 0.305 | 0.375 |
  | 3.0 | 0.279 | 0.349 | 0.290 | 0.355 |
  | 5.0 | 0.223 | 0.267 | 0.240 | 0.290 |

  → **「粘性+乱流で凝縮あり/なしを計算すると Fig.3 の差が出る」というユーザ仮説を定量的に確認**。
  比較図は `fig3_compare_visc.png` (非粘性破線 + 粘性SST実線 + 実験点)。

## SU2 クロスチェック (dry, 3 物理モデル) — `compare_su2_forge.py` / `compare_su2_forge.png`

同一ノズルの中心線 p/p0 を forge (SLAU, セル中心) と SU2 v8.5 (ROE, 節点) で比較。

| 物理 | forge run | SU2 run | forge↔SU2 一致 |
| --- | --- | --- | --- |
| Euler (非粘性) | `run_0001_slau_2d_imp` | `run_0020_su2_euler` | 全域 ≤0.6% |
| Laminar (粘性) | `run_0003_slau_2d_visc` | `run_0023_su2_lam_conv` | 全域 ≤0.6% |
| Turbulent (SST) | `run_0004_slau_2d_sst` | `run_0024_su2_sst_conv` | 全域 ≤0.8% |

→ **forge は 3 モデルとも SU2 と全域 <1% 一致**。粘性/SST は下流で等エントロピー線の上 (境界層排除厚) に乗り、
forge が SU2 と同等の BL 変位を出していることを確認。

⚠️ **教訓**: SU2 の初回 run (`run_0021_su2_lam` / `run_0022_su2_sst`) は前セッション中断で 975/752 iter で
SIGTERM 終了 (`su2.log` に `Exit Success` は出るが `rms[RhoE]≈-0.2` の未収束) しており、その未収束場と比べると
下流で偽の 20-25% 差が出ていた。途中解から継続収束 (`run_0023`/`run_0024`, `rms[RhoE]≈-1.4`, 出口積分量
ドリフト <0.2%) させると上記のとおり一致。**SU2 は `Exit Success` でなく iter 数・`rms[RhoE]` で収束判定すること**
([.github/forge-su2-cross-check.md](../../.github/forge-su2-cross-check.md) の収束確認節)。

## TP gas (thermalMethod 2) + 凝縮の検証 — 真因 = sub-200K NASA-9 外挿

ユーザ要望で TP 気相 + 凝縮を wys で計算しようとしたが dry TP すら発散。系統的切り分け
(`run_0025`〜`run_0042`、大半は破棄可の診断用) で**入口/壁/出口 BC・凝縮・種数・粘性・乱流・
時間積分・convMethod・float 精度・初期値をすべて棄却**。**真因は wys が出口 T≈159K (元の凝縮ゾーン ~27K) で
NASA-9 有効下限 200K を割り、低温外挿が TP を不安定化**すること。全温 Tt スイープで確定:

| run | 全温 Tt | 出口 T | 結果 |
| --- | --- | --- | --- |
| `run_0041_tp_n2_hot500` | 500 K | 278 K (>200K) | **400 step 完走・残差安定** |
| `run_0042_tp_n2_Tt360` | 360 K | ~194 K | step24 で発散 |
| (Tt=286.65, 通常) | 286.65 K | 159 K | step12 で発散 |

→ TP 実装は壊れていない (T>200K では wys でも動く)。**wys の極低温凝縮で CPG を強制してきた方針は正しい**。
詳細・全棄却リストは [`.github/plans/condensation-nonequilibrium.md`](../../.github/plans/condensation-nonequilibrium.md) の 2026-06-16 ログ。

## 壁面 y+ (SST メッシュ `nozzle_2d_sst.h5`)

`wall_yplus.py` で run_0009/run_0010 の壁第一層 y+ を収束場から実測 (y+=√(u_t·y1/ν), ν=分子粘性):

- 第一セル距離 y1 = **1.6〜9.0 µm** (スロート半高 2.25mm の ~0.1%)。
- **y+ ≤ 1 が壁面の 98.2%**、mean 0.61、max **1.88**、全面 **y+ ≤ 2**。
- 1 を超えるのは入口リーディングエッジ (x≈−6.3cm) のごく一部のみ。スロート下流 (超音速・凝縮領域) は
  y+ ≈ 0.5 で完全に解像 → **低Re SST の壁まで積分が成立**。`wall_yplus.png` に分布。
- **乱流粘性も正しく作動**: μ_turb/μ_lam は壁第一セルで ≈0 (粘性サブレイヤ、k→0/ω→∞ ゆえ)、対数層で
  median 4.5・max 29、コアで ~3。場全体で 79% のセルが μ_turb>μ_lam。第一セルが粘性サブレイヤ内ゆえ
  τ_w を有効粘性 (μ_lam+μ_turb) で評価しても y+ は不変 (= 壁解像が本物である証拠)。

(旧記述「近壁メッシュは Euler デモ用で y+~1 ではない」は、非粘性デモ用 `nozzle_2d.h5` に関するもので、
SST 専用メッシュ `nozzle_2d_sst.h5` には当てはまらない。上記実測のとおり壁解像済み。)

## 既知の課題

- 本ケース検証中に **RANS エネルギー方程式の乱流熱伝導欠落バグ**を発見・修正
  (`.github/plans/diffusion-turbulent-thermal-conductivity.md`)。修正前は近壁静温が
  449 K (全温 293 K 超過)、修正後 293 K に収束。

## 後処理スクリプト

- `postproc_centerline.py <run> <label>` — 中心線 Mach / 静圧を面積比等エントロピーと比較。
- `plot_residuals.py <run> <label>` — 全 rms 残差の片対数プロット。
- `wall_yplus.py <run> [--wall-physid 3]` — SST 壁第一層 y+ を収束場から算出・可視化。
- `compare_fig3.py` / `compare_fig3_visc.py` — Wyslouzil Fig.3 (凝縮あり/なし) と中心線 p/p0 を比較。
- `build_restart.py <src_res> <src_mesh> <dst_mesh> [--rok K --roomega W]` — 引き継ぎ初期場生成。

## 修正形状 (Fig.3) run 一覧 — run_0047 以降

正しい Fig.3 ノズル (`mesh/make_nozzle_fig3.py` → `nozzle_fig3_2d.h5`, 2D pseudo, front/back=slip) で再構築。
SLAU 陰解法 (timeIntegration:11, blockDPLUR, **nStepInner:5, cfl_pseudo:4**)。dry は実験 isentrope、
凝縮は実験 1kPa と比較。比較図: `dry_vs_exp.png` (laminar/SST vs isentrope)、`cond_models_compare.png` (4モデル vs 1kPa)。

| run | 物理 | 凝縮モデル | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0047_fig3_2d_lam_dry` | 2D laminar | off | dry。中心線 p/p0 が実験 isentrope と −2〜5% | active |
| `run_0048_fig3_2d_sst_dry` | 2D SST | off | dry。**実験 isentrope と ±1.5% (laminar より良)** → 凝縮は SST 上で実施 | active |
| `run_0049_fig3_2d_sst_hk` | 2D SST | Hertz–Knudsen | 全凝縮 (g~90%)。onset 過早・overshoot (0.44 vs exp 0.36 @2cm) | active |
| `run_0050_fig3_2d_sst_kwhk` | 2D SST | **Kantrowitz+HK** | **onset 遅延で実験最良** (Fluent UDF 知見と一致) | active |
| `run_0051_fig3_2d_sst_gyar` | 2D SST | Gyarmathy | 全凝縮。onset 過早・overshoot | active |
| `run_0052_fig3_2d_sst_kwgyar` | 2D SST | Kantrowitz+Gyarmathy | onset 遅延で実験良好 | active |

知見: **凝縮成長停止バグ修正後** (核生成を $r_{\rm nuc}=1.01r_*$ で生成し液滴を $r_*$ から離脱させる;
[docs/condensation/implementation.md](../../docs/condensation/implementation.md) §検証 case/16) に全モデルが全凝縮を再現。
下流 (x≳4cm) は全モデル実験 ±数%。**onset 域は Kantrowitz 有り (Kw+HK/Kw+Gyar) が overshoot を抑え実験に最良**。
残課題: x≈3cm のピークが実験よりやや高い (onset レート微調整)、3D 壁解像 (`make_nozzle_3d_wallres.py`) は未実施。

## 多成分 TP 発散の再検証 run 一覧 — run_0069 以降 (2026-06-18)

「多成分 thermally perfect gas (thermalMethod:2) が収束しない」件の原因切り分け。**結論 (0acca05 の
「sub-200K のみ・種数無関係」を訂正): 真因は 2 つ**。
1. **IC 生成器の運動エネルギーバグ**: `gen_ic_2sp_from_n2.py:69` が `ek=0.5*(roUx²+..)/ro` (solver 規約は
   `/ro²`)。高速域 (出口 |u|~500, ro~0.1) で roe が ~10⁵ J/kg 過小 → 読み戻しで多数セル T<200K に床張り付き。
2. **nSpecies≥2 コードパス固有の不安定**: ek を修正した整合 IC・全>200K・species 完全保存でも 2 成分は
   step6-7 で局所 ro<0 → NaN。真の単成分は同条件で安定。Y_H2O≈0 でも発散 (step15)、実 H2O で加速 (step7)。
3. sub-200K NASA-9 外挿は別の第 3 要因 (cold で全部落ちる。Tt スイープ `run_0041_hot500`[単成分]安定/
   `run_0042_Tt360`/`run_0043_Tt286` は単成分の話)。

棄却済 (再調査不要): float32 精度 (double 同一)・CFL/剛性 (explicit cfl0.05 でも発散)・implicit 分離更新。
全 run native ビルド (`.build-native/`)。**全て診断用・破棄予定**。

| run | 変更点 | 結果 (初 NaN) | 判定 |
| --- | --- | --- | --- |
| `run_0069_ver_min2sp_float` | 2sp, inviscid+explicit+laminar (最小化) | step4。ro<0・sonic→0 が最冷出口で発生 | 発散機構を特定 |
| `run_0070_ver_min2sp_double` | =0069 を **double** ビルド | step4、残差軌跡が float と6桁一致 | **①float精度 棄却** |
| `run_0071_ver_min_n2o2` | =0069 で H2O→O2 (IC不整合, 参考) | step3 | 参考 (IC不整合) |
| `run_0072_ver_min2sp_implicit` | =0069 を implicit | step4 | **④分離更新 棄却** |
| `run_0073_ver_2sp_zeroH2O` | 整合IC・**Y_H2O≈0** (実質純N2を2sp機構で) | step5 | **②H2O/種数 棄却** |
| `run_0074_ver_2sp_realH2O_consistIC` | 整合IC・Y_H2O=0.01095 | step4 | 2sp は発散 |
| `run_0076_ver_1sp_N2_stablecfg` | **単成分N2**・run_0056 と同一 config | step4 | cold で単成分も発散(交絡: ek+sub200K) |
| `run_0077_ver_2sp_stablecfg` | =0076 を 2sp に | step4 | cold は交絡で切り分け不可 |

注: cold (Tt=286.65) では sub-200K と ek バグの両方で 1sp も落ちるため切り分け不能。**hot (Tt=500K, 全>200K) が
クリーンな対照** (下表)。hot で 1sp は安定・2sp は発散し、種数効果が分離できる。

| run (hot Tt=500K) | 変更点 | 結果 | 判定 |
| --- | --- | --- | --- |
| `run_0081_ver_1sp_hot500_control` | **単成分N2**, res_400 から restart | **exit0・400step 安定・残差フラット** | 単成分 hot は安定 |
| `run_0080_ver_2sp_hot500` | N2+H2O, ek バグ IC | step7。res_0 で 9466 セル T<200 (IC破壊) | 要因1顕在 |
| `run_0083_ver_2sp_n2o2_hot500` | N2+O2 (生成エンタルピー≈0), ek バグ IC | step18 (序盤は残差低下) | O2 でも発散・H2Oより遅い |
| `run_0084_ver_2sp_hot500_fixedIC` | N2+H2O, **ek 修正 IC** | res_0 整合(全>200K, 残差=1sp並)・なお step7 | **要因2: 修正後も発散** |
| `run_0085_ver_2sp_hot500_zeroH2O_fixedIC` | **Y_H2O≈0**, ek 修正 IC | step15 (局所 ro<0・species完全保存) | **要因2はコードパス側** |

### 要因2 の特定 (run_0089〜0094, hot Tt=500K, ek修正IC, condensation無)

O2 と H2O の対照で **H2O 固有**と判明、cfl_pseudo 依存で **implicit 緩和ミスマッチ**と確定:

| run | 2nd種/Y | cfl_pseudo | 結果 |
| --- | --- | --- | --- |
| `run_0089_cmp_1sp` | なし | 2 | 40step 安定 |
| `run_0089_cmp_o2_1e6` | O2 / 1e-6 | 2 | **40step 安定** (h_O2−h_N2≈0) |
| `run_0089_cmp_h2o_1e6` | H2O / 1e-6 | 2 | step13 発散。H2O が入ったセルの T が偽上昇 |
| `run_0089_cmp_h2o_real` | H2O / 0.01095 | 2 | step5 発散 |
| `run_0091_fix_*` | (energy補正カーネル配線) | 2 | **悪化** (O2もstep12)=**二重計上**。デッドカーネルは修正でない |
| `run_0092b_h2o_cflp10` | H2O / 1e-6 | 10 | step2 |
| `run_0092b_h2o_cflp1` | H2O / 1e-6 | **1** | **60step 安定** |
| `run_0092b_h2o_cflp0p1` | H2O / 1e-6 | 0.1 | 安定 |
| `run_0093_h2o_real_cflp1_long` | H2O / 0.01095 | **1** | **400step 完走・T[288,504]・Y1保持** (SST残差はプラトー=単成分と同) |
| `run_0094_cold_real_fixedIC_cflp1` | H2O / 0.01095 **cold Tt=286** | 1 | step13 発散=**要因3(sub-200K)** が cold で残存 |

**結論 (3 要因)**:
1. **IC ek バグ** (`gen_ic_2sp_from_n2.py:72` を `/ro²` に修正済)。
2. **多成分 implicit 結合不安定**: roe は block-DPLUR、roY は別 point-implicit で更新→擬似時間緩和がミスマッチ→roY 変化に roe が対応せず Newton で T ジャンプ。H2O の |h|≈1.3e7 が増幅 (O2 は無害)。**回避策: 多成分 TP は `cfl_pseudo ≤ 1`**。恒久修正は species を block 同梱 or 緩和整合 (未実装)。対流流束自体は整合 (energy 補正配線で二重計上→O2 発散が証拠)。
3. **sub-200K NASA-9 外挿** (cold 固有): 出口<200K の極低温では CPG を使うか thermo 低温拡張。

`run_0070` 用 double バイナリは `.build-native/double/` (要 `FORGE_CUDA_BLOCKSIZE=128`)。デッド energy 補正カーネルは
[`.github/plans/thermophysics-multicomponent-tpgas.md`](../../.github/plans/thermophysics-multicomponent-tpgas.md) 参照。

### cfl_pseudo パラスタ + 等エントロピー IC + ramp (run_0095〜0100, hot Tt=500K, N2+H2O 1%)

「cfl を上げても収束するか」「初期値を等エントロピーに」の検証。**結論: cfl_pseudo≤1 は硬い上限**。

- **20step 刻み発達** (`run_0095`, cfl_pseudo=1, 600step): 場は定常・物理的 (T[288,504]・sub-200K=0・Y1=0.01095 保持・ro<0=0)。ただし残差はプラトー (rms_ro~1e-7, rms_roe~8e-3、機械ゼロには落ちない=単成分 run_0056 と同じ既存挙動)。cfl=1 vs 0.5 (`run_0096`) で最終場 max|ΔP|≈3.4kPa(~5%)・max|ΔT|≈13K=厳密収束前。
- **cfl_pseudo パラスタ** (warm-start IC, `run_0099_ws_cflp*`):

  | cfl_pseudo | 1 | 2 | 4 | 5 | 8 | 10 |
  |---|---|---|---|---|---|---|
  | 結果 | **安定(399)** | step6 | step3 | step2 | step2 | step2 |

- **等エントロピー IC** (`gen_ic_isentropic.py`: 混合 R/Cp/s0 から (Pt,Tt) 等エントロピー再構成。T[283,509]・sub-200K=0): 熱力学的には妥当だが **連続の式 (ρuA=const) 非整合** (軸速度を sqrt(2Δh) から出したため)。no-slip 粘性壁では壁せん断ショック、inviscid+slip でも質量不整合で **cfl=1 でも発散** (`run_0097`/`run_0098`)。warm-start (mass-consistent) の方が安定。
- **cfl ramp** (`run_0100_ramp_cflp*`: cfl=1 で落ち着いた `run_0099_ws_cflp1/res_400` から高 cfl 再起動): **cfl≥2 はやはり発散** (cfl=2→step7)。過渡限定でなく**構造的上限**。

**まとめ**: 多成分 TP は (1) CEA 有効域 (>200K) の高温条件 + (2) `cfl_pseudo≤1` で**安定・物理的な定常場**が得られる (打ち手は有効)。ただし残差はプラトー品質 (厳密収束・高 cfl には要因2の恒久修正=species を block 同梱/緩和整合が必要)。等エントロピー IC を使うなら ρuA 連続を満たす quasi-1D 構成が要る。

### エンタルピー基準オフセット (sensible-enthalpy datum) による安定化 (run_0101〜0109, 2026-06-19)

**要因2 (多成分 implicit 結合不安定) への打ち手**として、各化学種のエンタルピー基準を
`h_s(Tref=298.15K)=0` へ平行移動する **sensible-enthalpy datum** を実装 (config `physProp.thermoHrefTemp`、
solver は NASA-9 係数 `a7` に焼き込み・全経路整合、Python IC 生成器は `gen_ic_2sp_from_n2.py` 第6引数で同一基準)。
**狙い**: H2O の生成エンタルピー `h_H2O(298K)≈−13.4 MJ/kg` (N2/O2 は ≈0、`thermo_href_compare.png` 参照) が
roe(block-DPLUR)/roY(point-implicit) の擬似時間緩和ミスマッチを増幅して Newton 温度ジャンプを起こすため、
**桁違いの生成エンタルピーを除いて増幅を抑える**。非反応流では物理不変 (e_mix と Σh_s J_s* の基準移動が連続式で相殺)。

同一 hot N2+H2O(1%) 場 (`run_0081/res_400` 由来、全>200K)・同一メッシュ・同一 BC で、絶対基準 (abs) と
オフセット基準 (href) を `cfl_pseudo` スイープで比較 (全 native ビルド、要 `solverConfig.hpp` 変更後の full rebuild)。

| run | 基準 | cfl_pseudo | 結果 (初 NaN step) |
| --- | --- | --- | --- |
| `run_0101_abs_cflp1`   | 絶対       | 1  | **400step 完走** (control = run_0099 再現) |
| `run_0102_abs_cflp2`   | 絶対       | 2  | step6 発散 |
| `run_0103_abs_cflp4`   | 絶対       | 4  | step3 発散 |
| `run_0104_href_cflp1`  | オフセット | 1  | **400step 完走** |
| `run_0105_href_cflp2`  | オフセット | 2  | **400step 完走** ← abs は step6 で落ちる cfl で安定化 |
| `run_0106_href_cflp4`  | オフセット | 4  | step7 発散 (abs step3 より遅延) |
| `run_0107_href_cflp5`  | オフセット | 5  | step4 発散 |
| `run_0108_href_cflp8`  | オフセット | 8  | step3 発散 |
| `run_0109_href_cflp10` | オフセット | 10 | step3 発散 |

**結論 (打ち手は有効だが部分的)**:
1. **安定 `cfl_pseudo` 上限が 1→2 へ向上** (abs は cfl2 で step6 発散、href は cfl2 で 400step 完走)。高 cfl (4,5,8,10) でも発散ステップが一貫して後退 (cfl4: step3→7 等)。**生成エンタルピーの増幅が安定性に効いている直接証拠**。
2. **物理不変を確認**: 同一入力場の step0 再構成 (T,P,ρ,Ux) が abs と href で**機械精度一致** (T rel 1.2e-7, P 6.1e-8, ρ/Ux=0)。roe は基準分だけ平行移動 (`roe_abs−roe_href ≈ ρ·Y_H2O·h_ref(H2O)`, 整合 5%)。`roe/ρ` は abs の **−86〜−19 kJ/kg (負・桁落ち)** から href の **+61〜+127 kJ/kg (sensible・良条件)** へ。
3. **限界**: 上限は 2× 止まりで `cfl_pseudo≥4` は依然発散、残差はプラトー品質 (`check_convergence.py`=NOT CONVERGED、単成分 run_0056 と同等)。**要因2 の構造的ミスマッチ自体は未解消** — オフセットは増幅振幅を下げるだけ。厳密収束・高 cfl には species を 5×5 block 同梱 (full coupling) が引き続き必要。

図: `thermo_href_compare.png` (cp/h vs T・オフセット効果), `href_cfl_ceiling.png` (cfl 上限比較), 各 run の `residual_history.png`。
**run_0101〜0109 は診断用** (恒久保持は不要、結論は本表)。

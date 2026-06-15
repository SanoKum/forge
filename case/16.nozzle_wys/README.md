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

# 09. Taylor-Green Vortex (3D)

3 次元 Taylor-Green 渦 (TGV)。全方向周期境界の立方体 $[0,2\pi]^3$、構造化 hex $32^3$ (節点 $33^3$)。
低散逸スキームの **運動エネルギー (KE)・エントロピー総和の保存性** を見るための検証ケース。

- メッシュ: [`mesh/Taylor-Green.geo`](mesh/Taylor-Green.geo) (`nx=ny=33`, z 方向 32 層 extrude)。
- 初期条件: `initial: "Taylor-Green"` (無次元, $M_0=0.4$, $\rho_0=1$, $P_0=1/\gamma\approx0.714$)。
  $u=M_0\sin x\cos y\cos z$, $v=-M_0\cos x\sin y\cos z$, $w=0$。
- 数値: pure KEEP (`solver: KEEP`, `keepDissipation: 0` = Roe 散逸無し・非散逸中心流束のみ), 陽解法 RK4 (`timeIntegration: 4`), 乱流モデル無し。

## 重要な設定上の注意

- **無次元・低圧ケースなので EOS フロアを下げること**: 既定の圧力フロア `pMin=1.0` Pa は本ケースの圧力場
  ($P\in[0.65,0.78]$) を全域クランプし、初期場を破壊する (step0 で $P\equiv1.0$)。`physProp.pMin/roMin/tMin`
  を小さく (例 `1e-6`) 指定する。floor の config 化は [`plans/accepted/thermophysics-eos-positivity-floor-config.md`](../../plans/accepted/thermophysics-eos-positivity-floor-config.md) 参照。
- **node モードのメッシュは node 設定で変換すること**: `discretization: "node"` のとき `convertGmshToForge` は
  median-dual を構築して primal を置換する。cell 設定で変換した h5 を node run に流用すると primal hex 上で
  解いてしまい不正 (双対 nCV=35937 ≠ primal hex 32768)。各 node run dir で自身の node config を使って変換する。

## 計算 run 一覧

全 run 共通: pure KEEP は `keepDissipation: 0`、陽解法 RK4、乱流なし、`pMin: 1e-6`、出力間隔 50 step。

全 run で全運動量 $\sum\rho\mathbf{u}\,V$ (理論=0) と全質量を保存量チェックに用いる。**node の res は周期マージ双対
セルを重複格納する** (VALUE 長 35937 = ノード数、一意 CV は 32³=32768) ため、後処理は各境界 CV を多重度
$m=2^{(境界方向数)}$ で割って一意 1 個分に補正する (`plot_ke_entropy_history.py` の `cv_weight`)。cell は補正不要。

以下は **cell 周期 partnerCellID device 転送バグ修正後** ([`plans/accepted/boundary-cell-periodic-conservation.md`](../../plans/accepted/boundary-cell-periodic-conservation.md)) の結果。
保存量チェックは全運動量 $\sum\rho\mathbf{u}\,V$ (理論=0) と全質量。**node の res は周期マージ双対セルを重複格納する**
(VALUE 長 35937 = ノード数、一意 CV は 32³=32768) ため、後処理は各境界 CV を多重度 $m=2^{(境界方向数)}$ で割る
(`plot_ke_entropy_history.py` の `cv_weight`)。cell は補正不要。

| run_* | 離散化 | 設定差分 | 主要結果 (多重度補正後) | 状態 |
| --- | --- | --- | --- | --- |
| `run_0009_node_keep_pure_eul` | node (median-dual) | pure KEEP, **非粘性** | **質量厳密・運動量 ~1e-7・KE 0.3%・エントロピー ~1e-5 保存**。教科書通りの KEEP 保存性 | active ✅ |
| `run_0010_cell_pure_eul` | cell (primal hex 32768) | pure KEEP, **非粘性** | **修正後: 質量厳密・運動量 ~1e-7・KE 0.4%・エントロピー ~1e-5 保存**。node と一致 (修正前は step358 発散) | active ✅ |
| `run_0007_node_keep_pure_visc` | node (median-dual) | pure KEEP, 粘性 (Re≈160) | 質量厳密・運動量 ~1e-7 保存。**KE 物理減衰 K/K0→0.644、エントロピー +1.3e-2 (第二法則)** | active ✅ |
| `run_0008_cell_pure_visc` | cell (primal hex) | pure KEEP, 粘性 | **修正後: 運動量 ~1e-7 保存・KE 物理減衰 K/K0→0.661** (node とほぼ一致; 残差は primal/dual メッシュ差)。修正前は KE×8 スプリアス増殖 | active ✅ |
| `run_keep`, `run_keep_M0.1` | cell | 旧参照入力 (旧スキーマ config) | — | ref |

### KEEP 陰解法 (block-DPLUR dual-time) 物理CFL掃引 (2026-06-29)

KEEP を陰解法化するときの flux Jacobian 方針 ([`plans/active/convection-keep-revive-node.md`](../../plans/active/convection-keep-revive-node.md) §7)
の検証。**新規コードゼロ** (`solver: KEEP` + `dualTime: 1` + `blockDPLUR: 1`, `timeIntegration: 11`) で、LHS は
既存の Roe 分割ヤコビアン (`accumulate_split_jacobian_cf`) を流用。全 cell・非粘性・共通 T≈3.5。explicit は
物理CFL≈0.05 (dt=0.007)、implicit は dt を上げて物理CFLを掃引。掃引解析は
[`analyze_implicit_sweep.py`](analyze_implicit_sweep.py) → [`implicit_sweep_ke_entropy.png`](implicit_sweep_ke_entropy.png)。

| run_* | 物理CFL | 主要設定差分 | K/K0 (T末, 理想=1) | dS/\|S0\| | 状態 |
| --- | --- | --- | --- | --- | --- |
| `run_0011_cell_keep_expl_ref` | 0.05 | explicit RK4 (基準) | 1.0033 | -1.1e-6 | active ✅ |
| `run_0012_cell_keep_impl_cfl005` | 0.05 | implicit, dt=explicitと一致 | **1.0033 (explicitと一致)** | -5.7e-7 | active ✅ |
| `run_0013_cell_keep_impl_cfl02` | 0.2 | implicit dt=0.028 | 1.0032 | -8.0e-7 | active ✅ |
| `run_0014_cell_keep_impl_cfl05` | 0.5 | implicit dt=0.07 | 1.0032 | 1.8e-6 | active ✅ |
| `run_0015_cell_keep_impl_cfl1` | 1.0 | implicit dt=0.14 | 1.0031 | 1.1e-5 | active ✅ |
| `run_0016_cell_keep_impl_cfl2` | 2.0 | implicit dt=0.28 | 1.0020 | 4.7e-5 | active ✅ |
| `run_0017_cell_keep_impl_cfl4` | 4.0 | implicit dt=0.56 | 0.9983 | 1.8e-4 | active ✅ |
| `run_0018_cell_keep_impl_cfl8` | 8.0 | implicit dt=1.12 | 0.9865 | 6.1e-4 | active ✅ |
| `run_0019_cell_keep_impl_cfl16` | 16 | implicit dt=2.24 | 0.9542 | 1.8e-3 | active ✅ |
| `run_0020_cell_keep_impl_cfl2_sub50` | 2.0 | CFL2 で nSubIterDualTime 20→50 | 1.0020 (sub20と完全一致) | 4.7e-5 | active ✅ |
| `run_0021_cell_keep_impl_cfl16_sub80` | 16 | CFL16 で nSubIterDualTime 20→80 | 0.9536 (sub20とほぼ同一) | 1.8e-3 | active ✅ |

### keepDissType (ES 散逸レイヤ) L2 較正 run ([plan](../../plans/active/convection-keep-es-dissipation.md))

| `run_*` | 目的・設定差分 | KE drop (500step, M0.4, RK4) | 状態 |
| --- | --- | --- | --- |
| `run_0022_cell_keep_dissoff_l2` | keepDissType=0 (新バイナリで無散逸経路の不変確認) | **−0.33%** (run_0011 参照と同挙動, 差は atomicAdd ノイズ級) | ref (不変確認) |
| `run_0023_cell_keep_diss015_l2` | keepDissType=1, σ=0.15, lowMachPrecond=1 | **8.37%** (市松6桁減衰の対価) | ref (σ 較正上限) |
| `run_0024_cell_keep_diss005_l2` | 同 σ=0.05 | **2.71%** (市松~4桁減衰で十分 → **既定採用**) | ref (σ 較正・既定根拠) |

**所見**: (1) KEEP+block-DPLUR は **CFL≈16 (explicit の320倍) まで NaN なしで完走** — Roe-Jacobian LHS の流用で
陰解法が安定動作。(2) 物理CFL≤2 では explicit と KE/エントロピー保存が一致 (matched CFL.05 は K/K0=1.0033 で
explicit と4桁一致) → **LHS の Roe 散逸は収束解を汚染しない**。(3) CFL を上げると KE が単調減衰 (CFL16 で
-4.6%)・エントロピー増。**ただしこの劣化は内部反復不足 (LHS散逸漏れ) ではなく BDF2 dual-time の時間離散誤差**が
支配: CFL2 で内部残差は 6.8e-5→3e-9 と完全収束し sub20≡sub50、CFL16 でも sub20≈sub80 で KE が回復しない
(サブ反復を増やしても改善しない=温度離散誤差は irreducible)。→ **LHS の Jacobian 選択は解の質にほぼ無関係**で、
大 dt での精度限界は外側の時間積分 (BDF2) が決める。LES では物理 dt を時間精度が許す範囲に保つこと。

成果物: KE・エントロピー時間履歴 [`ke_entropy_history.png`](ke_entropy_history.png)、ポストスクリプト
[`plot_ke_entropy_history.py`](plot_ke_entropy_history.py) (旧 `plot_taylorGreen.py` の現行版=res ベース・多重度補正付き)。
各 run の VERDICT は `CONVERGENCE_VERDICT.txt` (非定常 TGV のため定常収束ツールは `NOT CONVERGED`=正常)。

### 結論

- **pure KEEP は cell・node とも KE/エントロピー/運動量を保存する** (修正後)。非粘性で KE 0.3–0.4%・エントロピー
  ~1e-5・運動量 ~1e-7、粘性で KE 物理減衰・エントロピー物理増加。cell と node が一致 (差は primal hex vs median-dual)。
- **根本原因 (修正済)**: cell 全周期で運動量が線形注入され seam が非周期化していたのは、**`bint_d["partnerCellID"]`
  (device) が host から未転送**だったため (`setPeriodicPartner` は host のみ充填、device コピーは yaml uniform 経路だけ)。
  periodic_d が未初期化 device partnerCellID を読み ghost が誤値に。`setPeriodicPartner` 直後に H2D コピーで根治。
  node は host の partnerCellID を直接使う (`buildPeriodicNodeGroups`) ため無傷だった。float/double・ghost 更新
  タイミングは無関係 (切り分け済)。
- 切り分けの決め手: ghost セル値 vs 解析 TGV(ghost 重心) が修正前 6112/6144 不一致 (誤差 0.43) → 修正後 0 (1.5e-7)。
- 補足: 過去 cell で「KEEP」と思っていた流束は legacy で `if(false)` 無効化されており実体は **MUSCL+Roe**。本物の
  KEEP 中心流束は `Revive KEEP` (79b4e67) で初めて有効化された (この periodic バグが顕在化したのもそのため)。
- 旧 `run_keep` が現バイナリで回らないのは solver 退行ではなく①config スキーマ進化②圧力フロア (`pMin` 既定 1.0 Pa
  が無次元低圧場をクランプ) であり、modern config + `pMin` 引き下げで解消する。

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

| run_* | 離散化 | 設定差分 | 主要結果 (多重度補正後) | 状態 |
| --- | --- | --- | --- | --- |
| `run_0009_node_keep_pure_eul` | node (median-dual) | pure KEEP, **非粘性** | 完走。**質量厳密保存・運動量 ~1e-7 (機械ゼロ)・KE 0.3%・エントロピー ~1e-5 で保存**。教科書通りの KEEP 保存性 | active ✅ |
| `run_0007_node_keep_pure_visc` | node (median-dual) | pure KEEP, 粘性 (Re≈160) | 完走。質量厳密・運動量 ~1e-7 保存。**KE 物理減衰 K/K0→0.64、エントロピー +1.3e-2 (第二法則どおり増加)** | active ✅ |
| `run_0008_cell_keep_diss_visc` | cell (primal hex 32768) | KEEP+Roe散逸, 粘性 | NaN は出ないが **全運動量が 0→(+105,−94,−15) と注入**され KE がスプリアス増殖 (K/K0→8.2)。**周期境界の非保存が原因** | active (不正例) |
| `run_0010_cell_keep_pure_visc` | cell (primal hex) | pure KEEP, 粘性 | step382 で **DIVERGED**。運動量注入 + 無散逸で即発散 | active (発散例) |
| `run_0012_cell_roe_visc` | cell (primal hex) | **ROE** (MUSCL+Roe=旧「KEEP」の実体), 粘性 | KE は一旦物理減衰するが**運動量はやはり注入** (t=7 で Px+11/Py−22)→step1039 NaN。**周期非保存はスキーム非依存** | active (切り分け) |
| `run_keep`, `run_keep_M0.1` | cell | 旧参照入力 (旧スキーマ config) | — | ref |

成果物: KE・エントロピー時間履歴 [`ke_entropy_history.png`](ke_entropy_history.png)、ポストスクリプト
[`plot_ke_entropy_history.py`](plot_ke_entropy_history.py) (旧 `plot_taylorGreen.py` の現行版=res ベース・多重度補正付き)。
各 run の VERDICT は `CONVERGENCE_VERDICT.txt` (非定常 TGV のため定常収束ツールは `NOT CONVERGED`=正常)。

### 結論

- **node (median-dual) の KEEP は厳密に保存的で正しい**: 非粘性で質量厳密・運動量 ~1e-7・KE 0.3%・エントロピー
  ~1e-5 保存、粘性で KE 物理減衰・エントロピー物理増加。pure KEEP の設計どおりの KE/エントロピー保存性を実証。
- **cell (collocated primal) の triply-periodic 境界が保存していない (実バグ)**: 周期 ghost 面の対流流束が
  反対称になっておらず、初期 (周期) 場ですら $\sum \mathrm{res}_{\rho u}\neq0$ → 全運動量が線形注入される。
  KEEP・ROE 両方で出るため**スキーム非依存・周期境界そのものの問題**。KEEP は低散逸ゆえ顕在化が早いだけ。
  → cell periodic flux の幾何/対応付けの修正が必要 (別タスク)。
- 補足: 過去 cell で「KEEP」と思っていた流束は legacy で `if(false)` 無効化されており、実体は **MUSCL+Roe**
  (`else` 分岐) だった。本物の KEEP 中心流束は `Revive KEEP` (79b4e67) で初めて有効化された。
- 旧 `run_keep` が現バイナリで回らないのは solver の退行ではなく①config スキーマ進化②圧力フロア (`pMin` 既定
  1.0 Pa が無次元低圧場をクランプ) であり、modern config + `pMin` 引き下げで解消する。

# 09.Taylor-Green (KEEP 保存性検証)

## 位置づけ

`09.Taylor-Green` は、**低散逸スキーム (純粋 KEEP) の運動エネルギー (KE)・エントロピー保存性**と、
**周期境界の保存性 (全運動量保存)** を確認するための検証ケース。対流スキーム (特に KEEP) や周期境界・
node/cell の共有コードを変更したときに使う。標準回帰 (`08`/`13`/`20`) が「定常場が基準とずれないか」を見るのに対し、
本ケースは**非定常の保存量 (KE・エントロピー・運動量) が物理的に振る舞うか**を見る。

3 次元 Taylor-Green 渦。全方向周期の立方体 $[0,2\pi]^3$、構造化 hex $32^3$。詳細・最新結果は
[`case/09.Taylor-Green/README.md`](../../case/09.Taylor-Green/README.md)。

## 検証対象 (何を見るか)

無次元設定 ($M_0=0.4$, $\rho_0=1$, $P_0=1/\gamma\approx0.714$) の TGV を pure KEEP・陽解法 RK4 で時間発展させ、
**全領域積分の保存量**を時間追跡する:

- **全運動量** $\displaystyle\sum_i \rho\mathbf{u}\,V_i$ — 周期箱では理論上 **0 を厳密保存**。崩れたら周期境界の非保存 (バグ)。
- **全運動エネルギー** $\displaystyle K=\sum_i \tfrac12\rho|\mathbf{u}|^2 V_i$。
- **全エントロピー** $\displaystyle S=\sum_i \rho\,c_v\ln(P/\rho^\gamma)\,V_i$。
- **全質量** $\displaystyle M=\sum_i \rho V_i$ (周期で厳密保存)。

期待される振る舞い:

| 構成 | KE | エントロピー | 運動量・質量 |
| --- | --- | --- | --- |
| **非粘性** (visc=0, thermCond=0) | ほぼ保存 ($K/K_0$ が数 % 以内・有界) | ほぼ保存 ($\sim10^{-5}$) | ~$10^{-7}$ (機械ゼロ) |
| **粘性** (visc=0.0025, Re≈160) | 物理減衰 ($K/K_0$ 単調減少) | 物理増加 (第二法則, $>0$) | ~$10^{-7}$ (機械ゼロ) |

pure KEEP が設計どおり KE/エントロピーを保存し (非粘性)、スプリアスな KE 生成や運動量注入を起こさないことが合否の核心。

## 計算方法 (準備・実行)

メッシュは `mesh/Taylor-Green.geo` (gmsh, `nx=ny=33`, z 32 層 extrude)。**cell と node でメッシュ変換が異なる**点に注意。

1. gmsh で `.msh` を生成: `gmsh -3 mesh/Taylor-Green.geo -o mesh/Taylor-Green.msh -format msh4`。
2. 各 run ディレクトリで `.msh` をコピーし、**その run の `solverConfig.yaml` を使って** `convertGmshToForge` で h5 化する。
   - **cell**: `discretization: "cell"` で変換 (primal hex, nCells=32768)。
   - **node**: `discretization: "node"` で変換すると median-dual を構築・置換する (双対 nCV=35937 = primal ノード数)。
     **cell 変換した h5 を node run に流用しない** ([[node-mesh-must-convert-with-node-config]])。
3. `forge` を実行 (run ごとに複製ディレクトリで)。

### 必須の設定 (solverConfig.yaml)

- `solver: "KEEP"` (純粋 KEEP・散逸なし)、`timeIntegration: 4` (RK4 陽解法)、`turbulence: {LESorRANS: 0}` (乱流なし)。
- `initial: "Taylor-Green"` ($M_0=0.4$)。`time.deltaT: {control: 0, dt: 0.007}`、`unsteady: 1` (uniform global dt)。
- **`physProp: {pMin: 1.0e-6, roMin: 1.0e-6, tMin: 1.0e-6}` を必須**。既定の圧力フロア `pMin=1.0` Pa は本ケースの
  低圧場 ($P\in[0.65,0.78]$) を全域クランプして初期場を破壊する ([[eos-pressure-floor-config]])。
- `bcondConfig.yaml` は全 6 面 periodic で、各面に並進ベクトル `floats: {dx,dy,dz}` ($=\pm2\pi$) を与える
  (旧スキーマの空 `floats:` では変換が通らない)。`probe.yaml` も必要 (空でよい)。

代表 run (修正後・pure KEEP):

- `case/09.Taylor-Green/run_0010_cell_pure_eul` (cell, 非粘性)
- `case/09.Taylor-Green/run_0009_node_keep_pure_eul` (node, 非粘性)
- `case/09.Taylor-Green/run_0008_cell_pure_visc`, `run_0007_node_keep_pure_visc` (粘性)

## 後処理・合否判定

[`case/09.Taylor-Green/plot_ke_entropy_history.py`](../../case/09.Taylor-Green/plot_ke_entropy_history.py) が
`res_*.h5` から $K(t)$・$S(t)$ を積分して履歴を描く (旧 `plot_taylorGreen.py` の現行版)。

- **node の res は周期マージ双対セルを重複格納する** (VALUE 長 35937、一意 CV は 32³=32768)。後処理は各境界 CV を
  多重度 $m=2^{(境界方向数)}$ で割って一意 1 個分に補正する (スクリプトの `cv_weight`)。cell は補正不要。
- **合否の目安** (非粘性): $|K/K_0-1|\lesssim1\%$、$|(S-S_0)/S_0|\lesssim10^{-4}$、全運動量 $\lesssim10^{-6}$、
  NaN なし完走。cell と node がほぼ一致すること (差は primal hex vs median-dual の離散化差)。
- **NaN/発散チェック**: 初手から発散していないか `residual_history.csv` を確認 (AGENTS 共通ルール)。

### 注意: 非定常ケースなので「定常収束」はしない

TGV は本質的に非定常で KE が物理的に減衰するため、`check_convergence.py` は `NOT CONVERGED` を返すのが**正常**。
保存性の合否は上記の積分量 (KE・エントロピー・運動量) で判断し、残差の下げ止まりを「未収束=異常」と解釈しない。

## node / cell 両方で確認する

対流流束・周期境界の共有コード変更時は cell・node の双方で回し、両者の保存量が一致することを確認する
([[verify-node-and-cell-both]])。周期境界の保存は **cell の partnerCellID device 転送修正**が前提
([boundary-cell-periodic-conservation](../../plans/accepted/boundary-cell-periodic-conservation.md))。

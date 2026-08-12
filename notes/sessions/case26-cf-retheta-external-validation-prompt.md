# [Codex 依頼] case/26 の検証基準見直し — 外部相関 $C_f(Re_\theta)$ の導入 (調査のみ)

> このファイルは Codex へコピペする作業プロンプト。実施後は結果を
> `plans/active/turbulence-node-wf-omega-source.md` と `case/26.flat_plate_sst/README.md` に反映する。

---

## 0. 依頼の要旨と制約

case/26 (乱流平板 SST) の検証基準を見直してほしい。node 壁関数問題の判定を、forge 自身の
cell 壁解像解との内部比較だけでなく、**外部相関 $C_f(Re_\theta)$ (Kármán–Schoenherr) を主基準**に
切り替える。

**制約 (厳守)**:

- **既存の保存済み結果だけで調査すること。新規計算を投入しない。**
- **ソルバ本体 (`solver_density_cuda/cuda_forge/`, `*.cpp`) を変更しない。**
- **例外: `solver_density_cuda/tools/check_quasisteady.py` への `cf_retheta` 量の追加だけは許可する** (§2.7)。
  case/26 系の相関ツール (`case/26.flat_plate_sst/tools/*.py`) の修正は今回やらない (提案に留める)。
- 成果物は「調査報告」+「`plans/active/turbulence-node-wf-omega-source.md` の検証・採否基準の更新」
  +「`case/26.flat_plate_sst/README.md` への検証基準・実測値の反映」。
- **新たな「真因確定」をしない。** 分離実験と収支が揃うまでは仮説として扱う
  (本セッションで既に 2 回、根拠不足の「真因」を出して撤回している。§7 参照)。

---

## 1. 主基準とする外部相関

NASA TMR の zero-pressure-gradient turbulent flat plate validation と同じ Kármán–Schoenherr 関係:

$$Re_\theta=\frac{\rho_\infty U_\infty\theta}{\mu_\infty},\qquad
\theta=\int_0^\infty\frac{\rho}{\rho_\infty}\frac{u}{U_\infty}\left(1-\frac{u}{U_\infty}\right)dy$$

$$C_{f,\mathrm{KS}}=\left[17.08(\log_{10}Re_\theta)^2+25.11\log_{10}Re_\theta+6.012\right]^{-1}$$

出典: <https://tmbwg.github.io/turbmodels/flatplate_val.html>

**呼称の統一 (今後の文書すべてで)**:

| 役割 | 対象 | 呼び方 |
| --- | --- | --- |
| 主たる絶対物理基準 | $C_f(Re_\theta)$ Kármán–Schoenherr | **外部相関** / **理論基準** (「厳密解」「真値」と呼ばない) |
| 補助的な絶対基準 | Schlichting $C_f(Re_x)$、普遍速度分布 | 補助基準 |
| 収支確認 | 壁面 $C_f$ と $2\,d\theta/dx$ | 収支 |
| 内部回帰基準 | forge cell 壁解像解 (`run_0007`) | **内部基準** (今後「真値」と単独で呼ばない) |
| SST 実装 verification | 条件を一致させた NASA TMR の SST/FUN3D/CFL3D 解 | verification 基準 |

---

## 2. 依頼事項

### 2.0 【前提ゲート】ZPG と「十分発達した乱流」であることを先に確認する

Kármán–Schoenherr は **zero-pressure-gradient の十分発達した乱流境界層**が前提。
case/26 は入口全圧・出口静圧で流れを作っているので、**名前だけで ZPG と仮定してはいけない**。
K–S を当てる前に次を調べ、**ゲートとして扱う**こと。

**(a) ZPG 判定**: 各 run で BL 外縁の $U_e(x)$、$p_e(x)$、および平板区間全体での
$\Delta U_e/U_e$ を出す。加速/減速があれば K–S 比較の前提が崩れる。
実質 ZPG とみなせるか、みなせないならどの x 範囲でならみなせるかを明記する。

**(b) 発達乱流判定**: $k(x)$、$\mu_t/\mu(x)$、$C_f(x)$ の streamwise 分布から
**「十分発達した乱流になった $x$」を特定**し、**それより下流だけを K–S 比較に使う**。
明示的な transition model は無いが、モデル自身の turbulence activation 位置が生じるため、
前縁直後は乱流として成立していない可能性がある。

**(c) 結論の書き方**: 有効域が短ければ、**K–S を唯一の合否基準にはできない**。
「条件を満たす station における主外部基準」と位置づけ、補助基準・収支・内部基準と併用する。

### 2.1 既存 run の $C_f(Re_\theta)$ を再計算する

最低限これらを比較 (パスは `case/26.flat_plate_sst/` 配下):

| run | 内容 | 使用 snapshot (本体) |
| --- | --- | --- |
| `run_0007_slau_muscl_innersweep` | cell 壁解像 (y+≈0.34) = **内部基準** | 6 個, 最終 `res_120000.h5` |
| `run_0017_ewt_yp30_wf` | cell y+30 壁関数 | 2 個, 最終 `res_10000.h5` |
| `run_0042_node_yp30_planar_2nd` | node y+30 基準 (第一層 1.768e-4) | 4 個, 最終 `res_20000.h5` |
| `run_0046_node_yp30_omgwfdir1` | node + `nodeOmegaWfDirichlet:1` A/B | 4 個, 最終 `res_20000.h5` |
| `run_0040_node_yp30_planar_2nd_long` | node y+ 依存 (第一層 3.400e-4) | 8 個, 最終 `res_40000.h5` |
| `run_0043_node_yp112_planar_2nd` | node y+ 依存 (第一層 6.588e-4) | 4 個, 最終 `res_20000.h5` |
| `run_0044_cell_yp30_wf_regr` | cell y+30 を**現行バイナリで**再走 (`wf_pk`/`Pk_diag` 出力あり) | 2 個, 最終 `res_5000.h5` |
| `run_0045_node_yp30_kwfdir0` | node + `nodeKwfDirichlet:0` A/B | 4 個, 最終 `res_20000.h5` |

各 run × 複数 x station で次を表にする:

$x$ / $Re_x$ / $\theta$ / $Re_\theta$ / 壁面 $C_f$ / $C_{f,\mathrm{KS}}$ / $C_f/C_{f,\mathrm{KS}}$ /
Schlichting 比 / $2\,d\theta/dx$ と壁面 $C_f$ の比

- **$d\theta/dx$ は 3 点差分でなく、保存されている streamwise 列を使った局所 fit か十分な点数の差分で
  評価すること。** 平板上には **199 station** ある (cell: x∈[0.00020, 0.99000], node: x∈[0.00040, 1.00000])。
- **前縁・出口近傍は除外**し、除外範囲を明記すること。

### 2.2 $C_f$ の定義を分離する

次の 4 つを混ぜないこと。表では**どの定義か列で明示**する。

1. **wall-resolved**: 分子粘性による壁面接線 traction ($\tau_w=\mu\,\partial u/\partial y$)
2. **wall-function**: ソルバが運動量式に課す modeled traction $\tau_w=\rho u_\tau^2$
3. **出力用 `twall`**: **独立検証量ではない可能性が高い** — §3.4 参照
4. **運動量積分**: $C_f=2\,d\theta/dx$

**node では、W–I 双対面に実際に加えられた接線方向の力から求めた $C_f$ が保存済み出力から取得できるか
確認すること。取得できない場合は「未確認」と明記し、`twall` を独立な証拠として代用しない。**

### 2.3 Kármán–Schoenherr の適用範囲を守る

NASA TMR の主比較域はおおむね $4000<Re_\theta<13000$。**この範囲に入る station を主判定**とし、
範囲外は参考値と明記する。

- **case/26 の平板長で有効 station が足りるかを調べること。**
  参考 (要検証の粗い見積り): README の x=0.6 で $\theta\approx1.03\times10^{-3}$ (cell) から外挿すると
  $Re_\theta$ は x=0.3 で ~2.6e3、x=0.6 で ~4.6e3、x=1.0 で ~6.7e3 程度。
  **上限 13000 には全く届かず、下限 4000 を超えるのは平板後半だけ**の可能性がある。
  足りない場合は「現ケースでは K–S 主比較域の下端しかカバーできない」と明記し、
  必要なら平板長/Re を変えた将来ケースの提案に留める (新規計算はしない)。

### 2.4 現状値をまず再確認する

README 記載値からの概算では x=0.6 付近で:

- cell y+30: $\theta\approx1.033\times10^{-3}$, $Re_\theta\approx4600$, $C_f/C_{f,\mathrm{KS}}\approx0.95$
- node: $\theta\approx8.75\times10^{-4}$, $Re_\theta\approx3900$, $C_f/C_{f,\mathrm{KS}}\approx0.75$

**これは概算なので実データから再計算すること。** 違っていれば、
使用した**定義・積分範囲・使用 snapshot** を示して訂正する。

#### $\theta$ 積分の規約 (これで統一すること)

$\delta_{99}$ は**境界層厚さの定義であって $\theta$ 積分の打切り位置ではない**。
$\theta$ は積分核 $\frac{\rho}{\rho_\infty}\frac{u}{U_\infty}(1-\frac{u}{U_\infty})$ が
**十分ゼロになる外部流まで**積分する。

- **原則として領域上端 ($y=0.2$) まで積分**し、**上端を下げても値が変わらないことを感度確認**する。
- 感度: 上端を $y=0.2$ / $0.1$ / $2\delta_{99}$ 程度で振り、$\theta$ の変化を表に出す。
- **node**: 壁点 ($y=0$, $u=0$) を**そのまま含める**。
- **cell**: 最下点がセル重心なので $(y,u)=(0,0)$ を**補って**積分する。
- **両者とも同じ $U_\infty$, $\rho_\infty$ の定義**を使う (BL 外縁の局所値でなく、
  §2.5 の自由流条件で固定するか、局所 $U_e$ を使うかを決めて明記する)。

### 2.5 現ケースと NASA ケースの条件差を確認する

最低限これらを比較表にする: $M$ / 単位長さ Reynolds 数 / SST variant / 流入 $k,\omega,\mu_t/\mu$ /
前縁の扱い / transition model の有無 / 上端・入口・出口 BC / plate length と得られる $Re_\theta$ 範囲。

**現ケースの実測値 (config から確定済み、再確認はしてよい)**:

- $M=0.2$, $P_t=100000$ Pa, $T_t=288.15$ K, 背圧 $P_s=97250$ Pa
- $U_\infty=67.782$ m/s, $\rho_\infty=1.1854$ kg/m³, $\mu=1.8\times10^{-5}$ Pa·s (定数, `viscMethod: 0`)
- $Re/m = 4.4636\times10^6$ → 平板長 1.0 m で $Re_{x,\max}=4.464\times10^6$
- 流入乱流 `k=0.3, omega=300` → **freestream $\mu_t/\mu = 65.9$**
- 領域: $x\in[-0.1,1.0]$, $y\in[0,0.2]$。底面は $x<0$ が slip (sym)、$x\ge0$ が no-slip wall。
  上面 slip、入口 `inlet_Pressure`、出口 `outlet_statPress`
- **明示的な transition model はない。乱流輸送方程式は全域で解くが、境界層が前縁直後から
  十分発達した乱流になっているとは限らない** (transition model 無しでもモデル自身の
  turbulence activation 位置が生じる。NASA の SST 平板も同様)

**★ NASA SST-Vm 標準ケースの freestream $\mu_t/\mu$ は約 0.009** で、現ケースの 65.9 と**4 桁違う**。
このため **NASA の SST 数値解との直接比較は apples-to-apples ではない**。
**K–S 相関との物理比較**と、**NASA SST 数値解による実装 verification** を分けて扱うこと。
後者をやるなら流入乱流量を合わせた別ケースが要る (今回は提案まで)。

### 2.6 相関式の不整合を調査する

リポジトリ内で 2 つの係数が混在している:

| ファイル | 係数 |
| --- | --- |
| `tools/postprocess_wall_law.py`, `tools/wall_law_modeled.py`, `tools/plot_uplus_truth_compare.py` | $0.0592\,Re_x^{-1/5}$ |
| `tools/compare_cf_bl_cell_node.py` | $0.0576\,Re_x^{-1/5}$ (コメントは「1/7 乗則」) |

出典と式の意味 (1/7 乗則の係数の流儀、$\delta/x=0.37Re_x^{-1/5}$ とセットの導出か等) を確認し、
**同じ「Schlichting」として扱ってよいかを整理**する。数 % の判定をするので係数差 (2.8%) を放置しない。
**今回は調査と計画更新まで。ツール修正は別途提案する。**

### 2.7 時系列の妥当性を確認する

使用する各 run について明示すること:

- `solver_density_cuda/tools/check_convergence.py` の VERDICT (`CONVERGENCE_VERDICT.txt` が
  `run_0007` 以外にはある。`run_0007` は無いので必要なら実行する — これは計算投入ではない)
- 報告する $C_f$, $\theta$, $Re_\theta$ の**準定常性**
- 使用した snapshot と run path

**`check_quasisteady.py` は `shock,asym,machmax,pmax` にしか対応しておらず $\theta$ / $C_f(Re_\theta)$ を
判定できない。** 一方 AGENTS.md は派生量の報告に同ツールの VERDICT を必須としている。
「ツール修正はしない」と「保存 snapshot から同等判定」は**この点で衝突する**ので、次で解消する:

- **`solver_density_cuda/tools/check_quasisteady.py` に `cf_retheta` 量を追加することだけ許可する**
  (再現可能性のため。他のツール修正は今回やらない)。追加した量で各 run の VERDICT を出すこと。
- **保存済み `CONVERGENCE_VERDICT.txt` を信用するだけでなく、`check_convergence.py` を各 run に
  再実行**して現行ツールの VERDICT を取り直すこと (`run_0007` は VERDICT ファイルが無い)。
- 「最後の 2 snapshot が近い」だけでは不十分。判定に使った窓幅を明示する。

> 注: 本 case の残差は block-DPLUR の近似ヤコビアン由来で**構造的にプラトー**し、
> 受理済み本番 `run_0007` を含め VERDICT は `NOT CONVERGED` になる。これは既知で、
> 派生量の準定常で判定する運用 (README 「収束判定」節)。ただし**その判定を今回きちんとやり直す**。

### 2.8 計画の採否基準を更新する

`plans/active/turbulence-node-wf-omega-source.md` の §7 (検証と判定) / §8 (完了条件) に少なくとも:

- **cell への接近だけを合格条件にしない** (現行 plan は既に「壁解像基準 1.000」に変えてあるが、
  外部相関ベースへさらに置き換える)
- **有効な $Re_\theta$ 範囲で $C_f/C_{f,\mathrm{KS}}$ を評価する**
- **壁面 $C_f$ と $2\,d\theta/dx$ の収支が閉じる**
- **case/26 の改善が case/40 の過大化を伴わない** (`nodeOmegaWfDirichlet` は case/40 で
  $\tau_w$/y+1 基準 1.237 倍に過大化して棄却された)
- **E1/E2/E3 は §7 の交絡を分離して判定する**
- **数値許容差は現状調査のばらつきを見て事前に明文化する** (後出しにしない)

---

## 3. 引き継ぎ: 既知の落とし穴 (ここで時間を溶かさないでほしい)

本セッションで実測により判明したもの。**いずれも「そう見えるが違う」系**。

### 3.1 `wf_pk >= 0` は壁ノードと第一内層ノードの**両方**に立つ

`ransWallFunction_d.cu:209` (`wf_pk[ic]`) と `:213` (`wf_pk[irep]`)。
「壁関数が効いている DOF」を選ぶマスクとして使うと壁ノードも巻き込む。

### 3.2 `nodeKwfDirichlet` の既定は **1** — 第一内層の $k$ は Dirichlet 固定

`solverConfig.cpp:460`。`roK_wf` に固定され `res_roK=0` にされる
(`ransSource_d.cu:236`, `update_d.cu:293`)。実測で `roK == roK_wf` が厳密一致。
→ **$\Sigma$`wf_pk`$\cdot V$ を「注入された生産量」と読んではいけない**。
`run_0042` では合計の **67.7%** がピン DOF 上で k 収支から消える。
(この誤読で「node は cell の 78% しか注入していない」という誤った結論を出し、撤回した。)

### 3.3 `nodeOmegaWfDirichlet` は $\omega$ ピンと limiter 迂回を**束ねている**

`turbulent_viscosity_d.cu:222`:
```cpp
vis_turb[ic] = wfPin ? rho * k_c / w_c
                     : rho * a1 * k_c / max(a1 * w_c, S_mag * F2);
```
→ このフラグの A/B は $\omega$ 単独の効果を測っていない。

### 3.4 `twall` / `res_wall_4_*.h5` の `utau` は**モデル目標値へ再スケール済み**

`viscousFlux_d.cu` の AddTauWall 再スケールは W–I 面の traction を $\tau_w=\rho u_\tau^2$ に
合わせ込む。したがって `twall` から出した $C_f$ は「課した値」であって独立検証にならない。
なお同再スケールは `tau_x *= scale` と **traction ベクトル全体**に掛けており、
**法線成分も同じ倍率 (本ケースで ~2.1) で増幅される** (コメントの「接線成分を再スケール」と不一致)。

### 3.5 node と cell は同じ primal メッシュでも $y^+$ が 2 倍違う

node の第一 DOF は**節点そのもの** ($y=h$)、cell の第一 DOF は**セル重心** ($y=h/2$)。
`run_0042` (node, $h$=1.768e-4) と `run_0044` (cell, $h$=3.400e-4) は
**第一 DOF の物理高さを揃えてある** (1.768e-4 vs 1.700e-4) 組。$y^+$ で揃えたい時はこれを使う。

### 3.6 $\theta$ 積分の下限・上限の扱いで値が変わる

node の列には**壁ノード ($y=0$, $u=0$)** が入る。cell の列は最下点がセル重心 ($y=h/2$) で
壁点を含まない。$\theta$ 積分の下限処理 (壁点を足すか、$u=0$ を補うか) と、
**上限 (BL 端の決め方: $u=\max$ の index か $0.99U_e$ か)** を揃え、明記すること。
node/cell 比較の 15% 差はこの処理に敏感な可能性がある。

### 3.7 座標の取り方が node と cell で違う

`res_*.h5` の `MESH/COORD` は node モードでは DOF と 1:1 (要素数一致)。
cell モードでは節点座標なので `MESH/CONNE` から重心を作る必要がある。
判定は `COORD.shape[0] == VALUE/ro.shape[0]`。
既存実装例: `case/26.flat_plate_sst/tools/plot_uplus_truth_compare.py` の `load()`。

### 3.8 `run_0017` は旧バイナリの出力で `wf_pk`/`Pk_diag` を持たない

現行バイナリの cell 対照が要るときは **`run_0044_cell_yp30_wf_regr`** を使う
($u_\tau$ は run_0017 と一致することを確認済み)。

### 3.9 平板 node メッシュは平面 2D (z 1 層) が正

押し出し (spanwise 2 ノード) メッシュは 2 次精度で spanwise 市松が無減衰になり壊れる
(`run_0032`–`run_0037` の発散の真因)。`run_0038` だけが押し出し版の 1 次。
`run_0040/0042/0043/0045/0046` は平面 2D。**混ぜて比較しないこと。**

### 3.10 `check_mesh_quality.py` は median-dual (node) の h5 に非対応

primal (cell) 変換した h5 に対して実行する。

---

## 4. 報告フォーマット (この順で)

1. **結論を最初に**: 「現時点で node/cell のどちらが外部相関にどの程度合うか」を 3–5 行で。
2. **定義**: 使った $C_f$ の定義、$\theta$ の積分範囲・下限上限処理、$d\theta/dx$ の評価法と除外範囲。
3. **表**: §2.1 の全列。run × station。$Re_\theta$ が K–S 主比較域内かを列で明示。
4. **適用範囲**: K–S 有効域に入る station がいくつあるか。足りない場合の扱い。
5. **条件差**: §2.5 の比較表。NASA SST 数値解と直接比較できない理由を明記。
6. **相関式の不整合**: §2.6 の整理と、ツール修正の提案 (今回は実施しない)。
7. **時系列の妥当性**: §2.7。run ごとの VERDICT・準定常判定法・窓幅・使用 snapshot・run path。
8. **計画変更点**: `plans/active/turbulence-node-wf-omega-source.md` の §7/§8 をどう書き換えたか (diff 要約)。

**新たな「真因確定」はしないこと。** 現時点の正確な言い方は:

> 第一内層ノードでは **$k$ が壁関数由来の Dirichlet 値に固定される一方、$\omega$ は free で
> 解像 strain による $P_\omega$ を受け、標準 SST limiter が適用される**。この組合せの周辺に
> 問題が局在するが、**状態量・limiter・源項の寄与は未分離**である。

---

## 5. 参考: 現在の plan と経緯

- 計画: [`plans/active/turbulence-node-wf-omega-source.md`](../../plans/active/turbulence-node-wf-omega-source.md)
  (作業順序: **調査 → plan 方針更新 → 実装前に `methods/` の現仕様・変更仕様を更新 → plan 最終化
  → 診断実装 → 2×2 factorial → E3**。AGENTS.md の開発フローに従い、実装前に methods を更新する)
- 経緯と実測値: [`case/26.flat_plate_sst/README.md`](../../case/26.flat_plate_sst/README.md)
  の「node × 2 次精度の発散」「反証と、第一内層ノード $\omega$ への局在化」節
- SST-2003 の $P_\omega$ 誤記 (NASA TMR で一次確認済): [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §3.1.1

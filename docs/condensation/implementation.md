# 非平衡凝縮 — forge 実装解説

理論は [theory.md](theory.md) を参照。本書は forge コードへの対応 (追加保存変数・移流・ソース・
二相 EOS・境界・段階実装) を記述する。実装計画は
[`.github/plans/condensation-nonequilibrium.md`](../../.github/plans/condensation-nonequilibrium.md)。

既存の RANS 2 方程式 (`roK`/`roOmega`) が「移流される追加保存スカラー + stiff ソース」の手本であり、
凝縮モーメントは**同じ骨格**で追加する。

---

## 1. 段階実装 (Phase 分け)

| Phase | 内容 | 状態 |
| --- | --- | --- |
| **Phase 1** | 4 モーメントを**受動スカラー** (ソース=0) として輸送する骨格。case/34 dry 回帰一致 | 本セッション |
| **Phase 2** | N2 凝縮物理 (CNT+Iland 核生成 / Goodheart 成長 / Hill $T_d$ / 二相 EOS / point-implicit ソース / $\mu_n$ 無次元化) | 次セッション |
| **Phase 3** | H2O (CNT+Kantrowitz / Hertz–Knudsen) + case/16 Wyslouzil 検証 | 後続 |

Phase 1 では `solverConfig` の `condensation`/`nCondSpecies` で機構を ON/OFF し、モーメントは
気相速度で運ばれる受動スカラー (拡散なし、ソース=0) として移流するだけ。ソース=0 ゆえ
$g\equiv0$ で気相に影響せず、dry 計算と場が一致することを回帰確認する。

---

## 2. 追加保存変数 (多成分一般化)

RANS の固定 2 本 (`roK`,`roOmega`) と違い、モーメントは**凝縮種ごと 4 本 × 可変種数**。化学種輸送の
動的登録 (`variables::registerSpecies`, [variables.cpp](../../solver_density_cuda/variables.cpp)) を手本に
`registerCondensation(nCondSpecies)` を新設し、種 $s$ ごとに次を `cellValNames`/`c`/`c_d` へ追加する
([variables.hpp](../../solver_density_cuda/variables.hpp) の roK/roOmega 構成を 1:1 踏襲):

- 保存量: `rog_<s>`, `roQ2_<s>`, `roQ1_<s>`, `roQ0_<s>` (Phase 2 で無次元 `roμn_<s>` に読み替え)
- RK ステージ: `*N_<s>`, `*M_<s>`
- 残差: `res_rog_<s>`, `res_rog_<s>_m` ほか
- point-implicit 対角: `src_jac_g_<s>`, `transport_diag_g_<s>` ほか
- primitive: `g_<s>`, `Q0_<s>` ほか (`dependentVariables` で $\rho\phi/\rho$ から復元)

`allocVariables` の前に 1 度だけ呼ぶ (`registerSpecies` の隣)。`output_cellValNames` に
`rog_<s>`/`g_<s>` 等を追加して可視化対象にする。`nCondSpecies<=0` のときは何も登録せず既存経路を保つ。

---

## 3. 移流 — 既存 ScalarTransportDesc の流用

汎用スカラ輸送コア
[`scalarTransport_d.{cu,cuh}`](../../solver_density_cuda/cuda_forge/scalarTransport_d.cuh) の
`ScalarTransportDesc` (1 次風上移流 + point-implicit `transport_diag` 対角 + 任意拡散) は物理非依存。
RANS の `ransTransport_d.cu` (`buildScalarDescs`/`ransTransport_d_wrapper`,
[ransTransport_d.cu:83-115](../../solver_density_cuda/cuda_forge/ransTransport_d.cu)) と同型の
**`condensationTransport_d.{cu,cuh}`** を新設する:

- 種 × 4 モーメントぶんの `ScalarTransportDesc` を構築し `scalarTransportResidual_d()` を回す。
- モーメントは分子拡散しない受動スカラーなので **`diffusion=0`** (SST の $\sigma_k/\sigma_\omega$ 拡散は不要)、
  `floor=0`。Phase 1 は `src_jac=0`。
- ゲート: `condensationEnabled(cfg) = (cfg.condensation==1 && cfg.nCondSpecies>=1)`。
- `assembleResidual` ([main.cpp](../../solver_density_cuda/main.cpp)) の RANS transport 呼出の直後に
  `condensationTransport_d_wrapper` を、explicit RK / segregated point-implicit に
  `condensationTimeIntegration_d_wrapper` を、roK/roOmega と同じ位置で差し込む。

### 対流フラックスと二相 EOS (Phase 2)

Phase 2 で圧力が $p=(\rho-\rho g)RT$ となり気相へ弱く逆結合する。**SLAU を初手**とし、圧力流束は
$p$ を二相 EOS で評価するだけ、$\rho g$ は質量流束に乗る受動スカラーとして扱う ($\xi_g$ を陽に
Jacobian へ入れない)。保存量での圧力微分:

$$
\frac{\partial p}{\partial\rho}=\chi+\kappa K,\quad
\frac{\partial p}{\partial(\rho u_k)}=-\kappa u_k,\quad
\frac{\partial p}{\partial(\rho E)}=\kappa,\quad
\frac{\partial p}{\partial(\rho g)}=\xi_g
$$

($K=\tfrac12|\mathbf u|^2$、$\kappa,\chi,\xi_g$ は [theory.md](theory.md) 5 節)。flux Jacobian では
運動量行に $n_k\,\partial p/\partial U_j$、エネルギー行に $u_n\,\partial p/\partial U_j$ が入り、
**既存の $(\gamma-1)$ が全て $\kappa$ に置換、$(\rho g)$ 列に $\xi_g$ が新規追加**される。Roe は固有構造が
$\kappa,\chi,\xi_g$ で変わるため一般 EOS Roe (Vinokur–Montagné 流) が必要で後段。

---

## 4. ソース項 (Phase 2)

`ransSource_d.cu` の point-implicit パターン (消散ヤコビアン `src_jac` を対角に組む,
[ransSource_d.cu](../../solver_density_cuda/cuda_forge/ransSource_d.cu)) を手本に
`condensationSource_d.cu` を新設。device 側で:

- 物性 $p_{sat}(T),\ \rho_l(T),\ L(T),\ \sigma(T)$ を**種ごと構造体** (`condProperties_<sp>`) から評価。
- 核生成 $J$ (CNT × 種ごと補正 enum: `CNT`/`CNT_Iland`/`CNT_Kantrowitz`) を $T$ で評価、$r_*$ 算出。
- 成長 $dr/dt$ (enum: `gyarmathy`/`goodheart`/`hertz_knudsen`) を Hill $T_d$ 反復で評価。
- ソース $S_{Q_n}=J r_*^n + n\rho Q_{n-1}\,dr/dt$、$S_g$ を各モーメント残差へ加算。
- 源項ヤコビアン $\partial S/\partial U$ を `src_jac` に積み (point-implicit)、stiff を陰化。

モデル選択は **enum + switch (device)**、係数は種ごと構造体。新モデル追加は enum と関数を足すだけ。

---

## 5. 一温度 二相 EOS の温度逆算 (Phase 2 — 実装済)

[dependentVariables_d.cu](../../solver_density_cuda/cuda_forge/dependentVariables_d.cu) の
**`thermalMethod==2` (thermally-perfect) 分岐**に二相 EOS を追加した
([condensationEOS_d.cuh](../../solver_density_cuda/cuda_forge/condensationEOS_d.cuh) の
`cond_T_from_e_onetemp`)。**凝縮ケースは `thermalMethod 2` が必須** (定比熱近似は使わない。気相
$e_v(T)$ は NASA-9 thermo を再利用。NASA-9/CEA は ~200K 未満は外挿だが N2 は $c_p$ 平坦で許容)。

### 一温度近似 (T_v=T_d=T)

初期実装は気相温度 $T_v$ と液滴温度 $T_d$ を分けず $T_v=T_d=T$ とする。**$T_d$ は輸送変数に
追加しない** (Hill $T_d$ 式・$(T_v,T_d)$ 2 変数局所 Newton は後続拡張)。`cond_T_from_e_onetemp` は
T 1 変数 Newton で、拡張時に 2×2 Newton へ置換できる設計。

### 混合内部エネルギーと温度反転

論文形 $e=(1-g)e_v(T)+g\,e_l(T)$、液相 $e_l(T)=e_v(T)-L(T)$ より $e=e_v(T)-g\,L(T)$。保存量から
$e_{in}=\rho e/\rho-\tfrac12|\mathbf u|^2$ を作り、

$$
G(T)=e_v(T)-g\,L(T)-e_{in}=0,\quad G'(T)=c_v(T)-g\,L'(T)
$$

をセルごとに Newton で解く。$g$ は総液相質量分率 $\sum_s \rho g_s/\rho$ ($g_s$ は device `rog` 配列、
`condensationInit_d` が構築)。$g=0$ で $G=e_v(T)-e_{in}$ となり `thermo_T_from_e` と一致 → **単相 TP に
厳密縮約** (実装は $g<10^{-12}$ で従来 `thermo_T_from_e` を呼び bit 同一を保証)。

### 圧力 (実在気体補正なし)

液滴は圧力を持たないとして $p=(1-g)\rho R_v T$ (van der Waals 補正は将来オプション)。
$\rho e$・$Ht$ も $e_{mix}=e_v-gL$ で再構成。frozen 音速は当面気相 $\sqrt{\gamma_v R_v T}$ (loose coupling 近似)。

### 検証 (host unit test 済)

- **(b) g=0 厳密縮約**: `cond_T_from_e_onetemp(g=0)` が `thermo_T_from_e` と T 完全一致 (diff 0、50–300K)。
- **(c) g>0 安定・物理**: $e_{in}=e_v(50\text{K})$ 固定で $g$ を 0→0.08 に上げると潜熱で $T$ 単調上昇
  (50→72K、論文 40→56K と同オーダー)、Newton 残差 ~$10^{-9}$、$p=(1-g)\rho RT$ 正。
- 既存 `thermalMethod 0/2` 単相経路は未改変 (Phase 1 の run_0004=thermalMethod 0 は不変)。

### 将来 (二温度拡張)

$e=(1-g)e_v(T_v)+g\,e_l(T_d)$ と Hill $T_d$ 式を組み、$(T_v,T_d)$ の 2 変数局所 Newton へ。
気相は当面 1 成分だが、多成分凝縮では $e=e_v-\sum_s g_s L_s(T)$、圧力の $R_v$ も気相混合で一般化する。

---

## 6. 連成・陰解法 — 分離 (segregated)

モーメント輸送は **NS (気相 5×5 ブロック) と分離して解く**。roK/roOmega・化学種と同じく
block-DPLUR 後段の point-implicit 更新に乗せる (`applySSTPointImplicit` 同型の
`applyCondensationPointImplicit`, [update_d.cu](../../solver_density_cuda/cuda_forge/update_d.cu))。
これは論文の fractional-step (均質対流部 → 凝縮ソース部) とも整合する。

point-implicit 対角 (時間 + 移流、Phase 1 は source=0):

$$
D = \underbrace{\frac{V}{\Delta\tau}}_{\text{時間項}} + \underbrace{\sum_f\frac{\max(\dot m_f,0)}{\rho}}_{\text{移流項}\,=\,\texttt{transport\_diag}},\qquad
\delta(\rho\phi)=\frac{\text{Res}}{D}
$$

Phase 2 の二相 EOS による気相逆結合 ($p$ が $g$ 依存) は密結合せず**毎反復の従属変数 $T/P$ 再計算で
渡す loose coupling**。圧力微分 $\xi_g$ の陰解法対角への取込みと気相+モーメント密結合 (block) は
収束不良時のフォールバックとして検討する。

---

## 7. 境界・初期・残差

- **境界**: 入口は dry (液相モーメント=0)。`ransBoundary_d.cu` の Dirichlet(=0)/Neumann(zero-grad)/
  wall(zero-grad) パターン ([ransBoundary_d.cu](../../solver_density_cuda/cuda_forge/ransBoundary_d.cu)) を
  手本に `condensationBoundary_d.cu` を新設。inlet=0 固定、outlet/slip/axis=zero-gradient、壁は当面
  droplet 付着無視で zero-gradient。
- **初期**: `setInitial.hpp` の H2D コピー対象にモーメントを追加 (既定 0)。restart 読込
  (`readValueHDF5`) は roK/roOmega 同様、無ければ 0 フォールバック。
- **残差**: `main.cpp` の `residualEquationNames()`/`gatherResidualSnapshot()` に `rms_rog_<s>` 等を
  条件付き追加 → CSV 列に出力。

---

## 8. 精度・無次元化 (Phase 2)

Phase 1 は source=0 ゆえモーメントは常時 0 で精度無関係。Phase 2 で
$\mu_n=Q_n/(N_{ref}r_{ref}^n)$ ($N_{ref}=10^{18}$/kg, $r_{ref}=1$ nm) を導入し、保存量を $\rho\mu_n$ と
する。**まず全 float**で進め、無次元化で桁を O(1) に抑える。$r=(\mu_3/\mu_0)^{1/3}r_{ref}$ 等の
モーメント間比も無次元量どうしで安定。float で桁落ちが顕在化したら、その時点でモーメントの
double 化 (混合精度 cond storage、`flow_float` とは別の `cond_float=double`) をフォールバックとして
導入する。

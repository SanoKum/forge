# SST 断熱壁の SU2 式熱的結合 (Taw primitive overlay) — experimental mode 2

## メタ

- **area**: `turbulence / boundary`
- **status**: `done`  <!-- 2026-08-11 ユーザ決定: mode 3 (defect-flux 閉包) を正式採用。
  断熱壁 SST 壁関数 (node) の生産閉包 = sstThermalWallFunction: 3。mode 1 (output-only) は
  比較 fallback として維持。ω irep ピン (nodeOmegaWfDirichlet) の正式化と壁 μt 監査は
  別途 follow-up (§10c) -->
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5(f)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) 同節
- **related_plans**:
  - [`turbulence-sst-adiabatic-taw-fluxmodel.md`](turbulence-sst-adiabatic-taw-fluxmodel.md)
    (前身: output-only fallback (§0) と W–I 注入初回試行の失敗記録 (§9)。本計画はその後継)
  - [`architecture-node-boundary-gradient-dof-only.md`](architecture-node-boundary-gradient-dof-only.md)
    (前提: GG/LSQ owner-state-only — 本計画でも**変更しない**)
- **created**: `2026-08-11`
- **owner**: `sano`

## 1. 目的

SST automatic wall treatment の断熱壁について、**SU2 v8.5.0 と同じ熱的結合** (壁 primitive 温度の
Taw 上書き → 内部粘性熱流束が Taw を端点温度として参照) を experimental mode
(`sstThermalWallFunction: 2`) として実装し、壁半 CV のエネルギー収支検証を通過したら正式採用する。
output-only 版 (mode 1, commit 23f03169) は安定な比較基準として維持し、直ちに置き換えない。

## 2. 前提の訂正 (何を root cause と仮定しないか)

初回試行 (前身 plan §9) の発散について、以下を**前提にしない**:

- 「Taw 端点置換そのものが原理的に誤り」— SU2 も勾配計算後に壁 primitive 温度を Taw へ置換し
  内部粘性熱流束に使うが、同じ崩壊は起きていない。
- 「T[W] への依存が消えることだけが root cause」
- 「forge over-relaxed 補正と SU2 corrected-gradient の差だけが root cause」— 発散点周辺の
  補正係数差は約 1.08〜1.13 倍で、単独で 15 K への崩壊を説明するとは考えない。

確定しているのは**旧実装では壁ノードが先に異常冷却した**こと (step1000 で T[W] min
1100.7→27.1/15.5 K) のみ。調査対象は Taw 熱伝導・壁せん断仕事・壁側 μt・SST 壁変数・
陰解法 Jacobian を含む**壁半 CV 全体**。

## 3. SU2 v8.5.0 の処理順 (コード確認済み, 2026-08-11)

`CNSSolver::Preprocessing` ([CNSSolver.cpp](../../.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp) L71–127):

1. `CommonPreprocessing` — 保存量から primitive `T_state` を再構築
2. `SetPrimitive_Gradient_GG/LS` — **`T_state` で勾配計算**
3. limiter
4. `SetTau_Wall_WF` (L798–, Nichols & Nelson 2004) — **勾配計算の後**に:
   - Crocco–Busemann (eq.11): `T_Wall = T_Normal / (1 + β·U⁺ − Γ·U⁺²)` を Newton ループ内で更新し、
     **`nodes->SetTemperature(iPoint, T_Wall)` で壁ノード primitive 温度を直接上書き**
     (保存量 rhoE は据え置き = primitive overlay)。断熱では q_w=0 → β=0。
   - White–Christoph y⁺ + Spalding 合成で U_Tau を Newton 収束 (`relax` 緩和,
     `MinYPlus` 未満は自動 off, 不収束時は Y⁺=30/EddyVisc=1/U_Tau=1 の safe fallback)
   - `Eddy_Visc_Wall = μ_lam·(1 + dY⁺_White/dy⁺ − κe^{−κB}(1+κU⁺+κ²U⁺²/2) − μ_lam,N/μ_lam,W)`
     を `EddyViscWall` に格納
   - `Tau_Wall = (1/ρ_w)(Y⁺·μ_lam/d)²` を `SetTau_Wall`
5. 内部粘性流束 (`CAvgGrad_Base` corrected-gradient) が端点温度として overlay 済み primitive
   (壁ノードでは Taw) を参照
6. 物理壁面: HEAT_FLUX 指定値 / 断熱 q_w=0
7. 壁は**運動量行のみ**強 Dirichlet (`DeleteValsRowi`)、エネルギー行は生かす

さらに `CTurbSSTSolver::SetTurbVars_WF` ([CTurbSSTSolver.cpp](../../.external/su2-src/SU2_CFD/src/solvers/CTurbSSTSolver.cpp) L540–):
Normal_Neighbor の k/ω を `ω_new=√(ω_vis²+ω_log²)`, **`k_new = ω_new·EddyViscWall/ρ`** で緩和つき
Dirichlet → **SU2 の壁/第一内点 μt は壁法則整合の `EddyViscWall` に収束する** (§10 監査の核心)。

forge mode 2 もこの順序に従う: T_state → GG/LSQ (owner-state) → Taw overlay → 粘性流束が overlay を参照。

## 4. 設計

### 4.1 設定モード (baseline を失わない)

- `sstThermalWallFunction: 0` — Taw 計算・出力とも無効
- `sstThermalWallFunction: 1` — **現行 output-only (生産 baseline, commit 23f03169)**。
  `Tsb=Taw_diag`、場には非介入
- `sstThermalWallFunction: 2` — **experimental SU2 coupled**。`Tsb=Taw_diag` + 勾配計算後の
  primitive overlay を内部粘性流束へ使用

検証完了後に mode 2 を正式化するか判断する。

### 4.2 GG/LSQ は owner-state-only を維持 (変更禁止)

- GG: 全 primitive の非周期境界面値は owner node の状態 `phi[W]`。`bvar`/`Tsb`/`Taw` を参照しない。
- LSQ: 実在する内部隣接 node のみ。境界疑似点を M にも RHS にも追加しない。
- 温度勾配は必ず `T_state[W]` から計算する。

### 4.3 Taw は勾配計算後の working primitive overlay

- 状態配列 `T[W]`・保存量 `roe[W]` は**上書きしない** (SU2 の SetTemperature と等価な効果を、
  状態配列を汚さず per-node marker で実現する — forge は SU2 と違い primitive を毎ステージ
  再構築しないため、状態配列上書きは次ステージへ漏れるリスクがある)。
- per-node 配列 `Taw_Prim_Overlay` (非対象 −1 / 対象 = `Taw_diag[W]`) を
  `applySstThermalWallFunction` (mode 2 のみ) が書く。概念的には
  `T_flux[i] = (overlay[i]>0) ? overlay[i] : T_state[i]`。
- `Tsb` を内部流束カーネルから直接参照しない (bvar の境界面配列と node DOF の対応を
  流束層へ持ち込まない)。

### 4.4 対象辺は SU2 と同じく全入射内部辺

wall–interior / wall–wall を含め、overlay 対象壁ノードを端点に持つ**全内部辺**が自然に
`T_flux=Taw` を参照する (SU2 の primitive endpoint 参照と同じ)。`Taw_Rep_Id` による代表辺限定は
復活させない。

### 4.5 温度 corrected-gradient は SU2 式

mode 2 の overlay 対象辺では、旧 forge の compact 項へ Taw を単純代入せず SU2 式を実装する:

```
d      = x_j − x_i                       (ノード間ベクトル)
g      = 0.5·(gradT_state[i] + gradT_state[j])   (算術平均, Taw 注入前の DOF 状態勾配)
ΔT     = T_flux[j] − T_flux[i]           (overlay 済み端点温度)
g_corr = g + (ΔT − g·d)·d/(d·d)
heatflux = k_eff_face · (g_corr·S)       (S は ic0→ic1 向き内部双対面ベクトル, 残差は従来 ±F)
```

ただし corrected-gradient 差だけを root cause と仮定しない (§2)。まず SU2 と式を一致させ、
残る差を収支として測定するための変更である。

### 4.6 面平均の SU2 差

SU2 は端点 primitive・物性を算術平均する。mode 2 の overlay 対象辺では
`gradT_face`・`k_eff_face` (=0.5(λ0+λ1)+0.5(cp0+cp1)·0.5(μt0+μt1)/Pr_t) を**算術平均**で評価する。
非対象辺は forge 従来 (f 補間 + over-relaxed) のまま — この scoping (SU2 は全辺同式) は
「SU2 との差」として明記し、必要になれば全辺 SU2 式の A/B を追加する。
viscous work の face velocity・cp・μt の平均方式差も同様に §11 診断の測定対象。

### 4.7 物理壁面と壁エネルギー行 (mode 2 でも不変)

- 物理壁面→W の境界半割面熱流束 = 0 (断熱, 既存 `adiabaticWall` 分岐)
- 壁 `res_roe[W]` は生かす / エネルギー行 decouple しない / `roe[W]`・`T_state[W]` を Taw へ
  ピンしない (SU2 と同じ構成)

## 5. §10 SST 壁関数整合の監査 (最重要容疑)

旧実装が SU2 と異なる最大の候補は corrected-gradient よりも**壁側の乱流熱伝導と壁せん断仕事**。
既存データ (崩壊点付近 node 6564):

| 量 | SU2 (run_0042) | forge Taw OFF | forge Taw ON (旧注入) |
| --- | --- | --- | --- |
| wall μt | 1.86e-4 | 9.53e-4 | 1.42e-3 |
| 第一内点 μt | 1.12e-3 | 8.86e-4 | — |

差は壁ノード側に集中。Taw 熱流束は `k_eff = k_lam + cp·μt/Pr_t` を使うため壁側 μt 差は必ず監査。
**構造的原因 (2026-08-11 特定)**: SU2 は `EddyViscWall` (White–Christoph 微分による壁法則整合 μt)
を `SetTurbVars_WF` で壁/第一内点の k/ω に焼き込む → 壁 μt = EddyViscWall。forge は壁ノード ω を
Menter ブレンドでピンするが **k[W] は解かれたまま** → 壁 μt = ρk[W]/ω[W] は壁法則と無関係な値に
なる。監査チェックリスト (SU2 と比較):

- [ ] wall-function active 判定と MINYPLUS (SU2 は y⁺<MinYPlus で自動 off / forge は無条件)
- [ ] Normal_Neighbor 選択 / U_t / u_tau / Tau_Wall
- [ ] **EddyViscWall 相当の壁 μt** (forge に無い — 導入要否を §11 の収支で判断)
- [ ] 壁・Normal_Neighbor の k/ω (SetTurbVars_WF: SU2 は両点 Dirichlet + 緩和 / forge は
  nodeKwfDirichlet で第一内点のみ, μt_wall は Reichardt g 由来)
- [ ] wall-function relaxation (SU2 relax / forge なし)
- [ ] Newton 不収束時 fallback (SU2 Y⁺=30/EddyVisc=1/U_Tau=1 / forge 층流則)
- [ ] Nichols–Nelson (White–Christoph+Spalding) と forge Reichardt の u⁺(y⁺) 差
- [ ] Crocco–Busemann 陰形式 `T_N/(1−Γ U⁺²)` と forge 線形形式 `T_rep + r U_t²/2cp` の差

Taw だけを移植し、せん断・渦粘性・k/ω を別モデルのまま混ぜたことが収支不整合の原因でないかを
§11 で確認する。

## 6. §11 凍結場収支診断 (動的 run より先に必ず実施)

run_0055 (または run_0058) の収束場を凍結し、状態を更新せず残差を 1 回 assemble
(実装は Python による流束再構成で可 — 本セッションで SLAU/粘性の再構成手法は確立済み)。
壁ノードごとに出力: node id / x,y,z / volume / T_state[W] / T_rep / Taw / ρ[W],ρ[I] /
μt[W],μt[I] / k_eff_face / u_tau / Tau_Wall / U_t。各 W–I 辺: 伝導熱流束 / 粘性仕事 τ·u /
総粘性エネルギー流束 / 辺 id・相手ノード・向き。壁ノード集計: Σq_cond / Σ(τ·u) / q_wall·A(=0) /
total res_roe[W]。

同じ凍結場で 4 方式比較:

- **A**: 通常の T_state (output-only baseline)
- **B**: Taw + 旧 forge over-relaxed 式 (初回試行の再現)
- **C**: Taw + SU2 corrected-gradient 式
- **D**: C + 乱流熱伝導を laminar/turbulent 成分に分離

→ 壁冷却の主成分 (corrected-gradient / turbulent conductivity / viscous work 不足 / wall μt /
別項) を特定する。**旧発散点付近 node 6564, 6759, 11829 は必ず個別表示**。

## 7. §12 陰解法 Jacobian の確認

SU2 はエネルギー熱伝導 Jacobian を組み FGMRES/ILU で解く。forge block-DPLUR は粘性対角を
ν_eff ベースの共通近似としており Taw overlay を線形化していない。確認項目:

- [ ] Taw coupled residual に対するエネルギー対角の符号と大きさ
- [ ] dT/d(ρE) を含む熱伝導 Jacobian
- [ ] overlay を lagged 固定値とした場合の wall/interior 微分
- [ ] Taw の T_rep/U_t/cp 依存を lag するか線形化するか
- [ ] 現 DPLUR 対角が熱拡散 stiffness を覆っているか

最初の動的試験は explicit または十分小さい CFL・warm restart・1 次から・数十 step で早期
NaN/壁温チェック。低 CFL 安定 × 通常 CFL 発散なら Jacobian/時間積分不整合として切り分ける。
**クランプで Taw や T[W] を止める方法は root cause を隠すため採用しない**。

## 8. §13 動的 A/B 検証

凍結場収支が理解できてから新規 run (例: `run_0059_node_yp30_taw_su2coupled_lowcfl`,
`run_0060_node_yp30_taw_su2coupled_prod`)。投入前に divergence-and-startup.md とメッシュ品質
VERDICT を確認。確認項目: 初期数十 step NaN/Inf なし / T_state[W] が低温床へ向かわない /
壁 ρ・μt・k/ω が暴走しない / 壁半 CV で Σ(q_cond+τ·u)+q_wall·A→0 (定常時) / 内部辺 ±F 保存 /
物理壁面 q_w=0 / Tsb=Taw / GG/LSQ が Taw 非参照 / check_convergence.py / check_quasisteady.py /
run_0058 (output-only) 比較 / SU2 run_0042 と primitive 壁温・保存エネルギー復元壁温・第一内点温度・
wall/first-node μt・Tau_Wall・k/ω・massflux・η を分けて比較。未収束 run 同士を「一致」と表現しない。

## 9. §14 受入条件 (mode 2 → 正式 mode 1 昇格の条件)

- [ ] 発散・低温床・異常密度なし
- [ ] 収束または少なくとも報告量が準定常
- [ ] 壁半 CV エネルギー収支が説明可能 (q_cond と τ·u の釣り合いを確認)
- [ ] GG/LSQ owner-state-only を維持
- [ ] 物理壁面 q_w=0 を維持
- [ ] SU2 との差が、壁温出力だけでなく第一内点・μt・TauWall まで説明可能
- [ ] OFF / output-only / coupled の違いが文書化されている
- [ ] methods・active plan・case README が同期している

満たすまで commit 23f03169 の output-only 版を生産 baseline として維持する。

## 10. §11 凍結場収支診断の結果 (2026-08-11 実施済み)

run_0055 (Taw OFF 収束場, res_12000) を凍結し、壁 221 ノードの半割 CV エネルギー収支を Python で
再構成 (GG 勾配は保存 `dPdx` と 99%ile 相対誤差 1.9e-3 で照合済み)。スクリプト:
`wall_budget_frozen.py` / `su2_wall_budget.py` (`case/40.nozzle_design_tool/` に保存)。

**方式別の壁ノード収支 (平均値, W, ノードへ入る向き正):**

| 方式 | q_cond | τ·u work | net | 備考 |
| --- | --- | --- | --- | --- |
| A: 現行 forge (T_state, f補間) | −4.00 | +3.61 | **−0.37** | baseline, ほぼ釣合 |
| B: Taw+旧 over-relaxed 式 | −12.9 | +3.61 | **−9.3** | run_0053 の再現。無補償ドレイン |
| C: Taw+SU2 corrected-gradient (forge壁μt) | −11.9 | +3.61 | **−8.3** | **B とほぼ同じ = 式差では救えない** |
| E: C + SU2壁μt + 算術平均速度 work | −7.32 | +5.44 | **−1.88** | SU2 場再構成と同値まで縮小 |
| (参考) SU2 自身の場での再構成 | −8.61 | +6.72 | **−1.88** | 過渡/再構成誤差レベル |

- 式差 (C0−A) と平均方式差 (A2−A) は各 ~0.08 W (1–2%) で微小。**旧発散の主因は式ではなく、
  (i) forge 壁側 μt 過大による k_eff 過大 (Δ≈−4.6 W)、(ii) f 補間 (壁側重み 0.668) による
  面速度過小→仕事過小 (Δ≈−1.8 W)、(iii) 壁ノードが Taw 近傍まで昇温する正常な過渡 (~−2 W)**。
- D 分解: 乱流伝導が支配 (turb −10.8 W vs lam −1.11 W 平均) → 壁側 μt が直接効く。
- 対象ノード (netC → netE): 6564: −10.8→−2.35 / 6759: −10.5→−2.15 / 11829: −14.5→−3.67 W。

**§10 監査の確定事項 (SU2 run_0042 flow_fix2 実測):**

| 量 (node 6564 相当) | SU2 | forge (OFF) |
| --- | --- | --- |
| 壁 primitive T (=Taw 上書き) | 1428.2 K | (Tsb=1409.7 K) |
| 壁 **保存量復元** T | **1385.7 K** | 1125.5 K |
| 第一内点 T | 984.9 K | 995.2 K |
| 壁 μt (流束実効値) | **1.86e-4** | 9.53e-4 |
| 壁 ρk/ω (リミッタ前) | 4.08e-3 | — |
| 第一内点 μt | 1.12e-3 | 8.86e-4 |

- SU2 の壁 μt 1.86e-4 の正体 = `SetTurbVars_WF` が壁ノードに書く k_wf (EddyViscWall 由来,
  ρk/ω=4.1e-3) を **SST ひずみリミッタが ~22 分の 1 に切った値**。forge は k[W] を解いたまま
  (ピンなし) のため壁 μt が 5 倍過大。forge 既存の Reichardt ν_t,wall=ν(1/g−1) で壁 k を
  Dirichlet すれば同じ連鎖で ~1.7e-4 を再現できる (見積り、SU2 比 ~10%)。
- SU2 の収束場でも壁エネルギー行は厳密には閉じない (rms[RhoE] が 10^1.66 で prtau — SU2 自身の
  残差履歴でプラトーを確認)。壁保存量 T_cons は Taw の ~40 K 下で静定しており、発散はしない。
- **結論**: mode 2 は「overlay + corrected-gradient + 算術平均」に加えて **壁ノード k_wf
  Dirichlet (SU2 SetTurbVars_WF の両点書込) が必須**。片方だけでは −8 W 級ドレインが残り
  旧発散を再現する。

## 10b. §12 短尺動的試験の結果と未解決の設計論点 (2026-08-11)

mode 2 実装後の warm restart (run_0055 res_12000 起点)・低 CFL (cfl_pseudo 0.5)・implicit 試験:

| run | 条件 | 結果 |
| --- | --- | --- |
| `run_0059_..._lowcfl` | mode 2, 200 step | NaN なし完走。ただし T[W] が −0.7 K/step で単調冷却継続 |
| `run_0062_..._lowcfl_long` | mode 2, 5000 step | T[W] min 1125→113.6 K へ単調降下 (旧 run_0053 と同軌道、EOS 床向き)。T_I はむしろ微降 (対流が持ち去り、BL 加熱で釣り合う経路は存在しない) |
| `run_0063_..._mutw02` | + 実験 env `FORGE_TAW_WALL_MUT_SCALE=0.2` (壁側 μt を SU2 実効値相当へ) | 減衰速度は半減するが依然単調降下 (T6564 幾何減衰 比~0.9, 漸近 ~500 K 見込み=非物理) |

- 壁ノード k_wf ピン自体は機能 (k[6564]=9142≈k_wf) だが、forge の Reichardt k_wf は
  log 層平衡 k=u_τ²/√β* と同スケールで、SU2 の収束壁 k (3681) より 2.5 倍大きい。
  SU2 の壁 μt 1.86e-4 は **μt=ρa₁k/S (Bradshaw ひずみリミッタ, 壁で F2=1) で正確に再現**
  (ρa₁k/S = 1.98e-4)。ただし SU2 の壁 k/ω 収束値自体が SetTurbVars_WF の公称式から乖離した
  反復力学の産物であり、クリーンな式では再現できない。
- **SU2 の「安定」の実像**: SU2 自身も壁エネルギー行は収束していない (rms[RhoE] が 10^1.66 で
  プラトー)。primitive 上書きにより壁 T_cons (保存量復元 1386 K) は**流束にも出力にも使われない
  zombie DOF** となり、微小 local dt でゆっくり漂うだけで場を汚染しない。forge mode 2 は
  T_state[W] が EOS・対流・次ステップに実効なので、同じ漂流が場を直接冷却する。
- **未解決の設計論点 (ユーザ判断待ち)**: 指示 §9 (roe[W]/T_state[W] を Taw へピンしない・
  エネルギー行を生かす) は SU2 の構造の正確な写しだが、SU2 の実挙動は「エネルギー行は生きて
  いるが収束せず、primitive 上書きが漂流を無害化している」。forge で同等の無害化を得る選択肢:
  (a) SU2 完全準拠 = EOS/対流/出力も overlay 温度を参照する (実質、壁 T_state の primitive 層
  置換 — §9 の禁止に抵触するか要判断)、(b) 壁エネルギー行の扱いを SU2 と変えて真に釣り合う
  閉包を設計する (output-only mode 1 はその一形態)、(c) 壁 μt モデルと仕事項をさらに詰めて
  net→0 の完全平衡を目指す (本診断では forge 側の残差 −2 W 級が式一致後も残った)。

## 10c. mode 3 — 保存的壁層エネルギー流束閉包 (defect-flux, 本命, 2026-08-11 起案)

§10b の設計論点に対するユーザ決定: SU2 完全準拠 (a) でも aggregation でもなく、
**「実温度 T_W と回復温度 Taw のずれに比例する保存的 W–I 全エネルギー流束」** を実装する
(codex レビューと合意済み)。

### 設計

定応力 Couette 層の恒等式 $q(y)+\tau u(y)=q_w=0$ (断熱) より、層が Crocco 平衡
($T_W=T_{aw}$) のとき W–I 双対面の全 (伝導+仕事) エネルギー流束は厳密ゼロ。よって対象辺の
「線形伝導 + $\tau\!\cdot\!u$ 仕事」の**合計**を defect 流束

$$F^{wl}_e = H_T\,(T_{aw,rep}-T_W)\,(\vec S_e\cdot\hat m),\qquad
H_T=\frac{\rho_{rep} c_p u_\tau}{T^+(y^+_{rep})}\;(\text{Kader 原式})$$

に置換し `res_roe[W] += F` / `res_roe[I] -= F` (±厳密保存)。運動量行 (AddTauWall)・物理壁面
$q_w=0$・GG/LSQ owner-state-only・壁エネルギー行 (生かす)・T/roe (ピンなし) は全て不変。

- **負帰還**: $\partial F/\partial T_W = -H_T A<0$。mode 2 の死因 (T_W 非依存=復元力ゼロ) を解消
- **±保存**: ピン方式の死因 (片側捨て=湧き出し) を解消。Taw が T_rep と動いても W–I 間ゼロサム
- **面積分配**: 双対 CV 面閉性 $\sum_e \vec S_e\cdot\hat m = A_{wall}$ が厳密 (実測 0.994、
  残り 0.6% は W–W' 辺射影と曲率)。多重辺・非構造へ自動一般化、二重計上なし
- **$T^+$ 誤差不感**: 平衡点 $T_W=T_{aw}$ は $T^+$ に依存しない (緩和速度のみ)。圧縮性ベルでの
  Kader +87% 課題が答えを汚さない
- **壁 μt 非依存**: k_eff が熱経路から消え、§10 の壁 μt 5 倍問題を熱閉包として回避
- 淀み ($u_\tau\to0$ / $y^+<0.1$): $H_T=\lambda_{rep}/y$ (純伝導) へ退避 (既存 energyWf と同ガード)
- y⁺ ブレンド (低層 resolved / buffer smoothstep): Kader all-y⁺ 形が $H_T\to\lambda/y$ に自動退化
  し case/40 壁は全域 $y^+\approx80\text{–}210$ のため**初期実装では省略** (contingency として設計のみ)
- 陰解法: $\partial F/\partial T_W$ をエネルギー対角へ point-implicit 追加可 (完全 Newton 不要)。
  $H_T A$ は置換前の線形伝導係数より弱いため剛性は増えない — まず対角追加なしで試験

### 凍結場 scheme F (run_0055, 実施済み 2026-08-11)

| 量 (221 壁ノード) | 中央値 | 範囲 |
| --- | --- | --- |
| $H_T A$ [W/K] | 0.030 | 0.011–0.046 |
| $F$ (T_W=1125 K 時, 壁加熱向き) [W] | +5.9 | +0.01–+9.4 |
| $A_{eff}/A_{wall}$ | 0.994 | 0.993–0.995 |
| 平衡オフセット $\varepsilon=R_{rest}/(H_T A)$ [K] | **0.06** | −5.1–+41 |

$|\varepsilon|<5$ K が 212/221、$\geq20$ K は 4 ノード (コーナー由来と推定、動的試験で確認)。
対象ノードの予測平衡壁温: 6564→1410.5 K / 6759→1412.5 K / 11829→1396.2 K (≈Taw ✓)。

### 検証結果 (2026-08-11 実施, 実装済み)

- **mode 0 回帰** (run_0064): run_0055 と node 再現性帯 (~2e-6 rel) で一致 — mode 0/1/2 経路不変
- **短尺** (run_0065, warm+cfl_pseudo 0.5, 200 step): T[W] 単調上昇 1125→1246 K (幾何減衰,
  外挿 ≈1400 K)、内点 T_I=986 K 不動 (ドレイン・過熱なし)、NaN なし — 負帰還設計どおり
- **生産条件** (run_0066, cfl_pseudo 4.0, 12000 step, warm): **T[W] が step1000 で 1418.2 K に
  プラトーし 12000 まで不変** (node 6564; Taw_diag 1416.5 K, オフセット +1.7 K ≈ 凍結場予測
  +0.9 K)。壁平均 1451.3 K ≈ Taw 平均 1450.9 K。**SU2 primitive 壁温 1428 K と 10 K 一致 —
  今回は出力上書きでなく解いた実状態同士の一致**。非壁場は run_0055 と rel mean ~1e-4
  (roK のみ 2.7e-3) で保存。剛性問題なし (対角追加不要のまま cfl_pseudo 4 安定)
- **残差床 (コーナー) の決着 (2026-08-11 追試)**: baseline 比 6–12 倍の残差床
  (rms_roOmega 0.14 vs 0.012) を診断した結果:
  - **凍結場 ε 外れ値は動的に自己解消**: コーナー壁ノードの実測平衡オフセットは +7.5/+7.6/−3.7 K 等
    (凍結場予測 +41 K → 実解は状態が自己調整)。最大は出口リップ node 14364 の +22 K のみ。
    しかも当該ノードは y⁺=220–285 で淀み分岐外 — 「淀み fallback 適用範囲」問題は実在しなかった
  - **真の内訳 = 出口リップ第一内層列の ω 残差リミットサイクル**: baseline に既存 (±0.5) の振動が
    ±8 に増幅されたもの。**16/14365 ノード、x∈[62.5,69.8]mm・壁近傍に完全局在**。状態は
    6000→12000 で rel 1.5e-4 のみ変化 (実質静止)
  - **壁 k ピン仮説は A/B で棄却** (run_0067: 壁ノード roK_wf Dirichlet を mode 3 へ拡張 →
    rms_roOmega 0.141→0.165 で改善せず、壁温・場は不変)。mode 3 は最小構成 (壁 k ピンなし) を維持
  - **生産量への影響なし**: `check_quasisteady.py` VERDICT **ALL STEADY** (machmax/pmax drift
    0.0%)。出口 massflux 0.204731 kg/s (baseline +1e-5 rel)、運動量+圧力積分 310.208 N
    (baseline +2.4e-5 rel)、時間 drift 1.4e-6
  - 結論: 残差床上昇は出口リップ特異点の既存微小サイクルの増幅で、収束・定常性・生産量の
    いずれの障害でもない。既知の局所アーチファクトとして記録し、mode 3 の採用判断は
    これをブロッカーとしない
- **リップサイクルの打ち手 (2026-08-11 追試, run_0068/0069)**: 振動ノードは全て第一内層 (irep)
  行 = 「k のみピン・ω 自由」の半拘束行だった。既存フラグ `nodeOmegaWfDirichlet: 1`
  (SU2 SetTurbVars_WF は k・ω 両方を第一内点に書く) で **サイクル完全消滅**:
  - run_0068 (mode 3 + ω ピン): rms_roOmega **0.141→9.4e-4** (baseline 0.012 の 12 倍良)、
    rms_roK 1.9e-6・rms_roUx 2.5e-6 と**全成分 baseline 以下の残差床**。壁温プラトー 1421 K 維持
  - **代償 = 近壁乱流閉包の変更**: BL の k 場が大幅に変わり (near-wall 24%, BL 外縁 ~2 倍)、
    推力積分 310.20→309.63 N (**−0.18%**, η 許容 ±0.3% 内)
  - **帰属分離 (run_0069, mode 0 + ω ピン)**: 推力 309.623 N ≈ mode3+ω ピン → −0.18% は
    **ω ピン閉包そのものの効果で mode 3 とは独立** (mode 3 単独の推力影響は +2e-5 のまま)。
    ただし mode 0 (冷壁) との組合せは残差悪化 (rms_roOmega 0.23) — ω ピンは mode 3 と
    組み合わせたときだけ全成分が締まる
  - **採用判断の構造**: (i) mode 3 単独 (リップサイクルは無害と定量済み) と (ii) ω ピン併用
    (残差床最良・SU2 の SetTurbVars_WF に最も忠実・ただし BL 閉包変更で要再検証
    [case/26 平板 Cf・case/40 η vs SU2/y+1]) は独立に選べる。ω ピンの正式化は別途
    検証系列を経ること
- **ω irep ピンの検証結果 → 正式化見送り (2026-08-11 追試)**: y+1 壁解像 (run_0052, 解像
  |twall|) を真値アンカーに τ_w=ρ_rep·u_τ² (コードが AddTauWall で課す momentum 実効値) を比較:
  **SU2 wf 0.985 / forge 生産 (mode 3, ピンなし) 0.945 / forge+ω ピン 1.237**。ω ピンは第一内点
  ω を SU2 側へ寄せる (3.6 倍ズレ→1.2 倍) が、**μt が SU2 の 3.3 倍に過大化** (SU2 は第一内点でも
  Bradshaw ひずみリミッタが μt を 2.7 倍カットするのに対し forge はリミッタ不発 — ひずみ評価差) →
  BL 過混合 → U_t(rep) 上昇 → **τ_w +24% 過大**。推力 −0.18% はその帰結 (摩擦過大) で悪い方向。
  リップ ω サイクルは無害と定量済みのため、**残差床と引き換えに壁摩擦の真値一致を壊す取引は不成立。
  `nodeOmegaWfDirichlet` は既定 0 のまま** (剛性対策オプションとして残置、使用時は τ_w +24% の
  代償を明記)。注意: 壁ノード密度で τ_w を計算すると逆の結論に見える (ON が SU2 に接近) —
  必ず ρ_rep 基準で比較すること。付随知見: (a) forge 生産の τ_w は真値 −5.5% で健全、
  (b) node × 亜音速平板 y+30 (case/26 run_0030) は数値レシピ未確立 (convMethod 1 発散・MUSCL は
  Cf 10% ドリフトのプラトー) — ω ピン検証には nozzle 真値比較で足りたが、平板系は別課題

### 恒久 follow-up (mode 3 で解決しないもの)

- **壁 k/μt の SU2 5 倍差** (§10): mode 3 では熱経路から外れるが、SST 場の近壁 k・運動量の
  法線応力成分・将来の壁 μt 依存モデルには残る未解決差。SU2 の壁 μt=ρa₁k/S (ひずみリミッタ)
  再現と壁 k の扱い (SetTurbVars_WF 両点書込) の是非は別途監査を継続すること。
- 等温壁 `sstEnergyWallFunction` は伝導のみ q_w 置換で仕事項を残す (厳密には
  $q_{face}=q_w-\tau u_f$)。圧縮性ベル +87% 過大の一部要因の可能性 — mode 3 と同じ総流束置換で
  直せるか要検討 (スコープ外 follow-up)。

## 11. 変更ログ

- `2026-08-11` — 起票 (ユーザ指示による方針転換: output-only を最終形とせず SU2 熱結合を
  experimental mode 2 として再実装する)。SU2 v8.5.0 の処理順・SetTau_Wall_WF・SetTurbVars_WF を
  コードで確認し §3 に記録。壁側 μt 差の構造的原因 (EddyViscWall の有無) を §5 に特定。
- `2026-08-11 (2)` — §11 凍結場収支診断・§10 監査を実施 (§10 節に結果)。mode 2 実装:
  `sstThermalWallFunction` 0..2 拡張 (`solverConfig`)、`Taw_Prim_Overlay` 配列 (`variables.hpp`)、
  `compute_wall_friction_sst_d` に tawCoupled (overlay 書込 + 壁ノード roK_wf 拡張)、
  `viscousFlux_d` に overlay 辺の SU2 corrected-gradient + 算術平均 k_eff/面速度。
  mode 0/1 は全ゲートで nullptr/不変。
- `2026-08-11 (3)` — mode 0 回帰 (run_0061: run_0055 と node 再現性帯 ~2e-6 rel で一致) と
  §12 短尺動的試験 (run_0059/0062/0063, §10b 節)。**mode 2 は現状 experimental のまま非採用**:
  壁 μt を SU2 実効値へ寄せても壁ノードの単調冷却が止まらず、SU2 の安定機構が
  「primitive 上書きによる zombie 保存量の無害化」であることを特定。§10b の設計論点 (a)/(b)/(c)
  のユーザ判断待ち。生産 baseline は引き続き mode 1 (output-only, 23f03169)。
- `2026-08-11 (4)` — ユーザ決定により §10c mode 3 (defect-flux 閉包) を起案→凍結場 scheme F 検証
  (ε 中央値 0.06 K)→実装 (`Taw_HTnx/y/z`, `compute_wall_friction_sst_d` の H_T 格納,
  `viscousFlux_d` の W-I エネルギー行置換, config 0..3)→動的検証 (run_0064/0065/0066)。
  **T[W]=1418 K プラトー・SU2 と実状態 10 K 一致・内部場保存を達成**。残課題 = コーナー局所の
  残差床上昇 (§10c 検証結果)。正式採用 (mode 1 置換) はコーナー処理決着後に判断。
- `2026-08-11 (5)` — コーナー処理の決着 (§10c「残差床の決着」節): ε 外れ値は動的自己解消・
  淀み fallback 問題は不存在・残差床は出口リップ 16 ノードの既存 ω サイクル増幅と特定。
  壁 k ピン仮説は run_0067 A/B で棄却 (mode 3 は最小構成維持)。`check_quasisteady` ALL STEADY・
  massflux/推力 baseline 比 ~1e-5 で不変。**mode 3 の技術的ブロッカーは解消 — 正式採用
  (生産レシピの mode 1→3 切替) はユーザ判断待ち**。

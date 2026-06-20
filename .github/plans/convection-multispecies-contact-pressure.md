# 多成分 TP 接触面の圧力振動 (face-state 整合化 → PEP 切り分け)

## メタ

- **area**: `convection`
- **status**: `draft`
- **related_docs**:
  - `docs/convection/theory.md`
  - `docs/convection/implementation.md`
- **related_plans**:
  - `.github/plans/species-in-dplur-session-prompt.md` (本課題の発端: contact limit-cycle を別 issue に分離)
- **created**: `2026-06-20`
- **owner**: `CFD Dev`

## 1. 目的

多成分 thermally-perfect (TP) 気体の高 CFL 計算で、He/空気 contact 混合層に残る
**擬似時間 limit-cycle (残差床 ~4e-6 @ cfl4)** の根本原因を確定し、最小の正しい修正を入れる。
現状の SLAU は face state が **mixed-order**: `ρ_f, p_f, u_f` は 2 次再構成だが `Y_{s,f}` は
セル値の 1 次で、そこから `R_mix, T, h, c` を再計算している。本計画は **(a) 1D 材料接触面ベンチで
原因を厳密に切り分け**、**(b) face composition を flow と整合 (2 次再構成・全用途で同一)** にして
標準保存形 SLAU を face-state の面で正す。**double-flux / 完全 PEP 流束は本計画では実装しない**
(必要性が示されたら次 issue)。

## 2. スコープ

- **やる**:
  - 1D He/空気 material-contact ベンチ (stationary/moving × sharp/smooth, CFL sweep) を追加。
  - species を flow と同じ face 位置へ 2 次再構成する診断モード (`speciesFaceReconstruction`)。
  - face composition を SLAU flow 流束 (thermo: R_mix,T,h,c,H) と species 流束で **同一**に使う。
  - face 組成の positivity (共通 θ_Y スケーリング) / simplex (Σ∇Y=0 接空間) / 非実現可能 fallback。
  - limiter 結合モード (`faceLimiterCoupling`: 変数別 / species 個別 / ρ,p,u,Y 共通) の比較診断。
  - 既定 off でビット不変・回帰 (step-0 完全一致)。
- **やらない (本計画外・将来 issue)**:
  - double-flux / frozen-γ / 左右セル別 energy flux / Fujiwara-Terashima PEP / KEEP-PEP 中心流束。
  - full `5+N_s` implicit block / chemistry 結合変更 / axisymmetric source 変更。
  - `axisTimestepBeta` (near-axis 安定化) は既存設定を維持し本計画で変更しない。
  - CFL=8 達成のためだけの過剰安定化 (評価軸は「所定残差・工学量誤差までの wall-clock」)。

## 3. 関連 docs と前提 (これまでに確定済みの診断)

`species-in-dplur-session-prompt.md` の結論: contact limit-cycle は **near-axis (軸対称) とは別機構**で、
**数値 limit-cycle** (Test4: 振幅が pseudo-CFL に ~100×依存, cfl1 で消滅)。本計画の前提となる追加診断:

- **化学種移流は既に 1 次風上** (`scalar_advection_first_order_d`)。SLAU thermo は face で `Y` をセル値 (1 次)
  使用 (`convectiveFlux_d.cu` の `thermalMethod==2` 分岐, `Tl=P_L/(ρ_L·Rmix_cell[ic0])`)。**= mixed-order**。
- **CPG (thermalMethod!=2) は単一 γ で `h=γ/(γ-1)·P_L/ρ_L`、組成を一切使わない** → CPG は組成-熱力学の
  clean control (mixed-order を持たない)。
- **reconstruction を contact で 1 次化すると振動は増える** (`FORGE_CONTACT_1ST` hard も
  `FORGE_CONTACT_BLEND` 連続も A_YHe 2.1e-2→3.7e-2)。chatter フリーでも悪化 → **mixed-order/MUSCL/
  limiter chatter は (この方向では) 反証**。古典の「多成分 contact 圧力振動」(Abgrall) と整合
  (1 次保存形は接触面で最悪)。
- **relax を下げると振幅低下** (D3, 0.7→0.4 で A_YHe 2.1e-2→1.2e-2) → 擬似時間/defect-correction も関与。

> 重要な未確認点: 上記は **flow を 1 次へ下げる**方向のみ。**species を 2 次へ上げて整合させる**方向
> (本計画 Case B/C) は未検証。両者は同じ mixed-order の異なる解消法であり、別挙動の可能性がある。
> これを 1D ベンチで決着させるのが本計画の核心。

## 4. 設計方針

### 4.1 順序の差分 (助言からの唯一の変更): 1D ゲートを先に

助言は「§1-3 (species 2 次再構成 + positivity/simplex/fallback, 実装が重い) を実装 → §4-6 で検証」。
本計画は **cheap-diagnostic-before-expensive-implementation** の原則で順序を入れ替える:

1. **先に 1D ベンチを既存コードで回す**: Baseline A (flow2次/species1次=現行) と **Case D (全量 1 次,
   `convMethod:0`)** を CFL sweep で比較。これは**実装ほぼ不要**で、原因の最上位フォークを決める:
   - **Case D (全量 1 次) で ε_p がほぼ消える** → 振動は reconstruction 由来 → §4.2 (species 2 次整合,
     Case B/C) を実装する価値あり (Gate 1 候補)。
   - **Case D でも ε_p が残る** → 標準保存 SLAU の material-contact PEP 不足 (Gate 2) → §4.2 は効かない
     見込み。重い実装を回避し、double-flux を次 issue へ。
   - 既存 case/28 証拠 (全量 1 次で悪化) は Gate 2 を示唆。1D はそれを confound 無しで厳密確認する。
2. 1D ゲートが「reconstruction 由来 (Gate 1)」を示した場合のみ §4.2 を本実装する。

### 4.2 species を flow と同じ face 位置へ 2 次再構成 (整合 face composition)

理論・式は `docs/convection/theory.md` の追記節へ (本 plan には差分のみ)。要点:

- 各内部面で全 species を左右セルから face へ再構成: `Y_{s,L}=Y_{s,i}+ψ_{i}∇Y_{s,i}·(x_f-x_i)` (右も同様)。
  既存 gradient/Venkat limiter を流用。config `speciesFaceReconstruction` (0=現行1次, 1=2次)、既定 0。
- **positivity (共通 θ_Y)**: `ΔY_{s,f}=Y_{s,f}^raw-Y_{s,i}`、`ΔY<0` 成分で `θ_s=(Y_{s,i}-Y_min)/(-ΔY_{s,f})`、
  `θ_Y=min(1,min_s θ_s)` を全 species に同一適用 (成分間の相対再構成を保つ)。
- **simplex (Σ Y=1)**: セルで `Σ∇Y_s` が丸め以上にずれないか確認 → 必要なら接空間射影
  `∇Y_s ← ∇Y_s - Y_{s,i}Σ_k∇Y_k` (**式は無検証採用せず、現行の独立 species 数・保存変数に合わせて導出**)。
  face 残差のみ `Y_{s,f}←Y_{s,f}/Σ_k Y_{k,f}` を許容し、renorm 前 `ΣY-1`・補正量・発生 face 数を記録。
- **fallback**: 実現可能性を満たせない face はその側だけ `Y_{s,f}=Y_{s,cell}` (局所 1 次, SU2-NEMO 流)。
  fallback 数を集計。
- **limiter 結合** (`faceLimiterCoupling`: 0=変数別, 1=species 個別 2 次, 2=ρ,p,u,Y 共通
  `ψ_common=min(ψ_ρ,ψ_p,ψ_u,ψ_Y)`)。**1 と 2 を比較し、共通 limiter を即 production 既定にしない**。

### 4.3 同一 face composition を全用途へ

再構成・制約後の `Y_{s,L/R}^f` を **thermo (R_mix,T,h,c,H 再計算)・SLAU flow 流束 (mass/momentum/energy)・
species 流束**で共通使用。**(T,h,c,H) を独立 MUSCL しない** (基本再構成は ρ,p,u,Y、secondary は EOS 再計算)。
species 流束は同一 face 組成の upwind: `F_{ρY_s}=ṁ_f·Y_{s,upwind}^f`。**`Σ_s F_{ρY_s}=ṁ_f=F_ρ` を
丸め内で診断/assert** (個別 clamp で壊さない)。

## 5. 実装ステップ

1. **コード survey 文書化** (`docs/convection/implementation.md` に現行 face-state 生成フロー: gradient →
   ρ,p,u MUSCL → Venkat → face 組成 (1次) → R_mix → face T,h,c,H → SLAU mass/mom/energy → species flux)。
   触る主ファイル: `cuda_forge/convectiveFlux_d.cu`, `gradient_d.cu`, `limiter_d.cu`, `scalarTransport_d.cu`。
2. **1D He/空気 material-contact ベンチ追加** (`case/` に新規 or `tests/`): stationary/moving × sharp/smooth,
   `p_L=p_R=p_0`, `u_L=u_R=u_0`, 組成・ρ・T のみ変える (`p_0=ρR_mix T`)。Baseline A / Case D を既存コードで実行。
3. **1D ゲート判定** (§4.1)。Gate 2 なら §4-6 を中止し double-flux を次 issue 化、結果報告。
4. (Gate 1 のとき) **species face 2 次再構成 + positivity/simplex/fallback** 実装 (`convectiveFlux_d.cu`,
   config `speciesFaceReconstruction`/`faceLimiterCoupling`)。`docs/convection/{theory,implementation}.md` 先行更新。
5. **1D で Case B/C/D を比較** (§5.2)、CFL sweep (§5.3)、指標 (§5.4) を表化。
6. **case/28 cfl4 で baseline / 整合版 / 共通 limiter を比較** (§6 指標)。`axisTimestepBeta` は既存維持。
7. **回帰** (§9)、docs 同期、commit 分離。

## 6. 検証

- **単体 / ビルド**: native build。1D ベンチは新規 case/test。
- **1D 指標** (初期 `p_0,u_0`): `ε_{p,∞}=max|p_i-p_0|/p_0`, `ε_{p,2}` (体積加重 L2), `ε_u=max|u_i-u_0|`,
  `ε_{ΣY}=max|ΣY-1|`, `Y_min/Y_max`, **全エネルギー保存 `ε_E=|Σ(ρE)V-Σ(ρE)^0 V|/|Σ(ρE)^0 V|`**。
- **case/28**: residual floor, contact `A_{Y_He},A_p,A_{u_r},A_T`, massflux, 軸上 p/T, ΣY=1,
  face fallback 数, renorm 数, contact 厚さ, **wall-clock/step**, **所定残差までの総 wall-clock** (主要評価)。
- **判定基準**: 下記ゲート (§7)。最大安定 CFL は評価軸にしない。

## 7. 判定ゲート

- **Gate 1 (mixed-order 主因)**: face 整合版で `ε_p, A_p, A_{u_r}` が baseline より **1 桁以上低下** →
  mixed-order が主因。double-flux/完全 PEP は不要の可能性大。
- **Gate 2 (PEP 不足)**: 2 次 species+整合 thermo **でも**、全量 1 次 Case D **でも** ε_p が残る →
  標準完全保存 SLAU の material-contact PEP 不足。**この場合のみ** Fujiwara/Terashima 型多成分 PEP 適合
  条件を次 issue で調査。
- **Gate 3 (limiter 非同期)**: 共通 limiter だけが大きく効く → limiter 非同期切替/chatter。hard sensor を
  即 production 化せず連続 limiter coupling を検討。
- **Gate 4 (defect-correction)**: relax / inner iteration だけが効く → 高次 RHS と近似 LHS の
  defect-correction limit-cycle。full PEP より **LHS 改善を優先**。

> 現状証拠は Gate 2 + Gate 4 を示唆 (1 次化で悪化=PEP, relax で低減=defect-correction)。1D ベンチで確定する。

## 8. 影響範囲

- 触る: `cuda_forge/convectiveFlux_d.cu` (SLAU 経路), config (`speciesFaceReconstruction`,
  `faceLimiterCoupling`), 新規 1D ベンチ, `docs/convection/{theory,implementation}.md`, `docs/index.md`。
- 既定 off で既存ケース不変 (単成分 CPG/TP, 非軸対称, 既存多成分, case/23, case/28 baseline)。
- near-axis (`axisTimestepBeta`) と独立。

## 9. 完了条件 / 回帰

- [ ] `docs/convection/theory.md` 更新 (mixed-order, 整合 face composition, positivity/simplex)。
- [ ] `docs/convection/implementation.md` 更新 (現行フロー survey + 変更)。
- [ ] 1D ベンチ + CFL sweep + 指標表、case/28 比較、fallback/renorm 統計。
- [ ] **回帰**: `speciesFaceReconstruction=0` で step-0 完全一致 (face state/dt_local/residual/flux)。
      単成分 CPG/TP・非軸対称・既存多成分・case/23・case/28 baseline がビット不変
      (atomicAdd 非決定性のため長時間完全一致は不要、step-0 一致で no-op 証明)。
- [ ] 判定ゲートの結論 + 次の推奨を報告。
- [ ] commit 分離: `test: add multispecies material-contact benchmark` /
      `feat: reconstruct consistent multispecies face composition` / `test/docs: validate contact pressure`。
- [ ] `.github/plans/README.md` 同期、本 plan `status` 更新、§10 に変更ログ。
- [ ] **コード変更後は自動で PEP 実装へ進まず、結果と解釈を先に報告する。**

## 10. 変更ログ

- `2026-06-20` — 初稿。診断フラグ (`FORGE_CONTACT_1ST/BLEND/LOG`, commit `5dcdc80`) と既存切り分け結果
  (mixed-order をこの方向では反証, Gate 2/4 示唆) を反映。助言からの差分は §4.1 (1D ゲートを §4.2 実装の前に置く)。

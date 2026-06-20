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

## 11. 実装ログ (S2/S3 face 整合再構成)

- **periodic_d 修正 (commit 9ff8471)**: ghost を partner セル全状態 (roY/TP thermo 含む) でコピー。
  多成分 periodic の seam 破綻 (12%) を解消 → 壁反射なしのクリーン moving-contact ベンチ。
- **Phase A+B (commit 7ceaade)**: ∇Y_s (Green-Gauss) + config `speciesFaceReconstruction` (既定0=ビット不変)。
  flux TP 分岐で Y_L/Y_R を ρ と同じ limiter (limiter_ro) で face 再構成 (clamp+正規化) し thermo (R_mix,T,h) に使用。
  - **結果**: thermo-only S2 は contact 圧力振動を**減らさない** (periodic He/N2, ε_p 比 ~1.1-1.33 で僅悪化、
    中心補間 S2 と同様)。ΣY=1 維持。**これは thermo だけで、energy 流束は再構成 Y・species 質量流束は 1 次 Y の
    不整合を生む**ため。完全な判定には Phase C (S3) が必要。
- **Phase C (S3) の設計課題 (未実装)**: species 流束も同一 face 組成にするには ① res_roY の zero/集計が
  `speciesTransport` (convectiveFlux の後) にあるため、flux 内集計は上書きされる → 残差組立順の再構成が要、
  または ② scalar advection 側で同一再構成を複製 (interp_dispatch が convectiveFlux ローカル) が要。
  どちらも中規模。thermo-S2 が効かない現状を踏まえ、S3 実装の是非を判断ポイントとする。

## 12. case/28 検証 + S3 判断 (2026-06-20)

- **1D 再評価 (sin でなく tanh smooth contact, 5-10 セル, 周期, periodic_d 修正後)**: S2 (thermo 整合) は
  contact 振動を**減らさない**(ε_p 比 ~1.1-1.4, 一貫)。
- **case/28 cfl=4 で S0 vs S2 直接検証**:
  - S0: rms_ro=4.63e-6, A_YHe=2.8e-2, A_ur=7.6, A_P=2808
  - S2: rms_ro=**3.23e-6 (↓30%)** だが A_YHe=7.6e-2, A_ur=19.7, A_P=6094 (**contact 振幅は↑**)
  - → S2-thermo は global 残差を下げるが contact 振動を上げる = **thermo/advection 不整合**(energy 流束は再構成 Y・
    species 質量流束は 1 次 Y)を生む。1D と整合。**S3 (移流も同一 face 組成) が必須。**
- **S3 の restructure コスト**: species 移流と implicit `transport_diag` が `scalarTransportResidual_d` で
  一体。flux 内で advection を組む (exact consistency) には res_roY zero 順 + diag 温存の再構成が要
  (case/28 は speciesImplicitCoupling:1 で diag 必要)。中規模。
- **判断**: thermo-only が一貫して悪化する一方 global 残差は下がる両義性。S3 を実装して完全整合で判定するか、
  シグナルを受け 別路線 (Gate 2 double-flux / Gate 4 高 CFL implicit) へ倒すか、要判断。

## 13. 自走セッション結果 (2026-06-20, S2/S3 + ψ_Y 実装・検証)

ユーザ就寝中の自走で実装・検証。**結論: S2 (thermo 整合) は case/28 で効く。S3 (移流整合) は flux は正しいが
高 CFL で LHS 不整合により発散 → species LHS の 2 次化が次の鍵。**

### 大発見: S2 は長時間 settled では効く (前の「悪化」は restart 過渡だった)
case/28 cfl=4, settled (step 5000-10000), `speciesFaceReconstruction:1` (S2 thermo, ρ-limiter):
- S0: rms_ro=4.36e-6, A_ur=2.9, A_P=1252
- **S2: rms_ro=2.64e-6 (-40%), A_ur=2.2, A_P=638 (-50%)** ← 圧力振動が半減・残差も低下
- bounded: max|ΣY-1|=1.2e-7, Ymin=0, Ymax=0.767。**= mixed-order 30K 補正が実問題で効く。**
- 1D explicit synthetic は依然 S2 僅悪化 (高 CFL implicit 多次元の case/28 を代表しない)。

### ψ_Y リミタ: min(ψ_ρ,ψ_Y) は発散 → ρ-limiter (同一リミタ) が正しい
- ユーザ指摘どおり ψ_Y (`limiter_Y`, Venkat on ∇Y) を実装 (`limiter_d.cu` で limiter_r1_d 流用)。
- だが Y 再構成に **min(ψ_ρ,ψ_Y)** を使うと case/28 で step~1141 発散。理由: ρ は ψ_ρ・Y は min(<ψ_ρ) で
  **ρ-Y のリミタが食い違い再び mixed-order**。→ **Y にも ψ_ρ を使う (ρと同次数=consistent) が正解**で安定+有効。
  case/28 では ρ-limiter でも Y は bounded のまま (overshoot 顕在化せず)。robust な代替は全変数共通
  min(ψ_ρ,ψ_p,ψ_u,ψ_Y) (faceLimiterCoupling=2, flow 再構成も変える大改修) で将来課題。ψ_Y は診断用に残置。

### S3 (移流も同一 face 組成): flux は正・高 CFL で LHS 発散 (AI 予言的中)
- 実装: convectiveFlux が upwind 再構成 face 組成を `Yface[ip*nSpecies+s]` へ書き出し (ΣY_f=1)、
  `species_advection_faceY_d` がそれで移流。res_roY zero/diag 順は不変 (clean)。config
  `speciesFaceReconstruction:2` = S2+S3。
- **保存 OK** (1D: ΣY=1 を 6e-8 で維持), **低 CFL 安定** (1D + case/28 cfl1,2 NaN なし)。
- **だが case/28 cfl=4 で発散** (step~1289)。原因 = **2 次 species 流束 (RHS) と 1 次 upwind の
  transport_diag (LHS) の defect-correction 不整合**。AI の警告どおり。→ **species LHS を 2 次/整合化**すれば
  高 CFL も通る見込み。現状は **S3 は cfl≤2 専用**、S2 は cfl=4 で production 可。

### production 推奨
- **`speciesFaceReconstruction:1` (S2, ρ-limiter)** = cfl=4 で安定・圧力振動半減。**実用候補。** 既定 0 でビット不変。
- `:2` (S2+S3) は cfl≤2 か、species LHS 2 次化後に高 CFL。
- commit: 7ceaade (A+B) / 21b92f5 (ψ_Y+settled 発見) / [S3 commit]。

### 起きたら見るべき (走らせ中の run)
- `run_s2b` (cfl4 S2, 10000): production 候補の settled 確認。
- `run_c2_s0/s2/s3` (cfl2, 6000): S0/S2/S3 比較 — S3 が S2 を更に改善するか (cfl2 で全安定)。

## 10. 変更ログ

- `2026-06-20` — 初稿。診断フラグ (`FORGE_CONTACT_1ST/BLEND/LOG`, commit `5dcdc80`) と既存切り分け結果
  (mixed-order をこの方向では反証, Gate 2/4 示唆) を反映。助言からの差分は §4.1 (1D ゲートを §4.2 実装の前に置く)。
- `2026-06-20` — **1D ゲート実施・結論 (§4.1 ステップ3)**: `case/05` に He/N2 material-contact ベンチ
  (`gen_contact_ic.py`, p,u 一様・組成のみ界面) を追加し Baseline A (2次) / Case D (1次) を比較。
  - **stationary contact**: 両スキームとも完全保持 (ε_p ~ round-off 7.8e-8, ε_u=0, ε_ΣY=0)。
  - **moving contact (u0=100)**: contact 近傍の ΔP は小さく**減衰** (377→45 Pa, ~0.4%, 1次 Case D も同程度)。
    大きな ε_p≈0.2 は **slip 壁反射の人工物** (max|dP| が x=0.94-0.98 の壁、contact x~0.53 ではない)。
  - **正しい判定 (限定的)**: 「**鋭い 1D He/N2 material-contact では現行 SLAU が大きな古典 pressure
    oscillation を発生させない**」ことのみ確認。**mixed-order 仮説は反証していない** (重要訂正):
    sharp contact 近傍では limiter が強く効き flow MUSCL も実質 1 次化 (`ψ_ρ,ψ_p,ψ_u≈0`) していた可能性が
    高く、Baseline A も contact では実質 1 次 → Case D と同程度なのは当然。**この試験は mixed-order を励起
    できていない**。
  - **再評価が必要 (§4.3 へ)**: 有限厚さ smooth composition layer (∇Y≠0, ψ>0) で **S0(現行 mixed-order)/
    S1(全 1 次)/ S2(thermo 用 Y のみ face 整合)/ S3(species も 2 次)** を比較し、face で `ΔT_f^MO=
    T_f^mixed-T_f^consistent` を直接計測して、pressure/energy-flux 振動と同位相か確認する。case/28 でも
    D2a を 1 本残す (実混合層は ∇Y≠0・ψ>0)。
  - **現時点で言えること**: double-flux を**直ちに必要とする根拠は弱まった**が、**species face 再構成 /
    thermo face 整合化が無効とはまだ言えない**。Gate 4 (高 CFL implicit defect-correction; Test4+D3) も
    並行候補だが、mixed-order の励起試験 (smooth layer S0-S3) を先に決着させる。
  - status: 1D sharp は負結果だが mixed-order 未反証。§4.3 (smooth layer 再評価) を追加し継続。
- `2026-06-20` — **smooth layer 再評価 (訂正)。mixed-order を「反証」とした前項を撤回**:
  - **sharp contact の confound を実測確認**: contact 面で `ψ_ρ(limRo)=0` (Venkat が密度再構成を抹消) →
    Baseline A も実質 1 次 → Case D と同等で当然。**sharp 試験は mixed-order を励起できていなかった**。
  - **smooth layer (δ≈4 セル, ∇Y≠0) では `ψ_ρ=1`** (再構成 active)。ここで **mixed-order の face 温度誤差
    `ΔT_f^MO = T_f(R_mix(Y_cell)) - T_f(R_mix(Y_face))` を直接計測 → 最大 30K (T~300K の 10%)**、
    \|ΔT_f^MO\|>10K の face が多数。**mixed-order 誤差は大きく、無視できない** (前回の「反証」は誤り)。
  - **S0 (mixed) vs S1 (all-1次) @ smooth moving**: S0 は ε_p 1.5e-3→**2e-4 に減衰**、S1 は ε_p~1.8e-3
    **持続**・ε_u が ~1.05 に増大。→ S0 (2次 flow) は S1 より良い (case/28 D2a が悪化したのと整合) が、
    S0 にも残差床 ~2e-4 (round-off 2.3e-7 の ~1000倍) が残り、これが mixed-order 由来かは S2/S3 で要検証。
  - **正しい現在地**: mixed-order face-state 誤差は**実在し大きい (30K)**。ただし S0(2次 flow) は all-1次より
    良いので「flow を 1 次へ落とす」方向 (D2a) は逆効果。**正しい修正方向は S2/S3 = species/組成を flow と
    **同じ face 位置へ整合再構成**して 30K の誤差を消すこと**。これは未実装・未検証。double-flux を直ちに要する
    根拠は弱いが、species face 再構成/thermo face 整合化が無効とは**言えない** (むしろ最有力候補に復帰)。
  - **次**: S2/S3 実装 (`calcGradient_d.cu` に ∇Y 追加 + Venkat ψ_Y + flux で Y_L/R 整合再構成 +
    positivity/simplex + thermo/flow/species 流束で同一 face 組成)。smooth layer で S0/S1/S2/S3 比較 +
    case/28 で D2a(off/on) と S2 を比較。診断: `FORGE_CONTACT_LOG` に `ψ_ρ,ψ_p,ψ_u` と `ΔT_f^MO` 追加済 (commit 後続)。

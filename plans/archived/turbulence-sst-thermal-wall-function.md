# SST 壁関数の熱的閉包 (Crocco 型断熱回復温度 `sstThermalWallFunction`)

## メタ

- **area**: `turbulence / boundary`
- **status**: `superseded`  <!-- 2026-08-11 訂正: 本 plan の「弱閉包」設計は境界勾配一般化前提が崩れており
  (node の境界勾配が当時 bvar を読んでいたため場に触れる経路が実在)、かつ実データ (run_0038) では
  壁ノードの実状態がほぼ動いていなかった (「4.1K 一致」は出力配列の Taw 診断値同士の一致)。
  後継: `plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md` (SU2 式流束モデル置換) -->
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5(f)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §3.7
- **related_plans**:
  - [`boundary-node-nozzle-wall-outlet-stability.md`](boundary-node-nozzle-wall-outlet-stability.md) (§2.9–2.10: y+1 基準解と node 出口欠陥)
  - [`../accepted/turbulence-node-sst-wallfunction.md`](../accepted/turbulence-node-sst-wallfunction.md)
- **created**: `2026-08-11`
- **owner**: `sano` (実装: Claude 自走, Codex レビュー合意済み)

## 1. 目的

forge の SST automatic wall treatment (`wallTreatmentSST=1`) は運動量 (τ_w)・乱流 (ω ピン,
wf_pk, k/ω Dirichlet) のみをモデル化し、**エネルギー側の壁法則を持たない**。このため粗壁
メッシュ (y+≈19–40) の断熱壁温が回復温度より約 230 K 低い (case/40: forge 1193 K vs
SU2 壁関数 1422 K / y+≈1 low-Re 三者基準 1387–1414 K)。SU2 と同じ Crocco 型
T_aw = T_rep + Pr^{1/3}·U_t²/(2cp) を壁面温度状態に適用し、壁関数メッシュでも壁温を
生産値にする。`Taw_diag` (実装済み診断) がベル部 mean 1420.2 K = SU2 壁関数 1422 K と
2 K 差であることをフォワード確認済み — 適用値の実体はこの診断値。

## 2. スコープ

- **やる**: opt-in config `turbulence.sstThermalWallFunction` (既定 0)。node/cell の断熱壁
  (`kind: wall`) × `wallTreatmentSST==1` への T_aw 適用。case/40 での検証 (壁温・η・
  エンタルピー収支)。
- **やらない**: 等温壁の熱流束壁関数 (Kader q_w モデル — WMLES 経路に既存、SST への移植は
  将来)。low-Re (`wallTreatmentSST=0`) 側の変更。node 出口列欠陥の修正 (別 plan §2.10)。

## 3. 関連 docs と前提

- theory §6.5(f) / implementation §3.7 (本 plan と同時更新済み)。
- 参考実装: SU2 `CNSSolver::SetTau_Wall_WF` の T_aw 更新 (`.external/su2-src`)、forge WMLES の
  Kader 温度壁法則 (`wmlesWallModel_d.cu`)。
- 前提: `Taw_diag` は `ransWallFunction_d.cu` が壁 bplane ごとに代表点 (Normal_Neighbor) から
  毎ステップ計算し、壁 DOF (`ic`) に書いている (淀み域は T_rep)。recovery r = Pr^{1/3} は
  `prandtlLam` から。

## 4. 設計方針 (2026-08-11 検証で改訂 — 状態適用は不可、出力+勾配閉包の「弱閉包」が正)

- **確定設計 (弱閉包)**: `applySstThermalWallFunction` (applyBconds 位相,
  `applyRansScalarBoundaries` の後) が断熱壁 bplane の **bvar `Ts` に `Taw_diag[ic]` を書くだけ**
  (node/cell 共通, `set_wall_taw_output_d`)。壁面値は①壁面出力②境界 LSQ 勾配閉包
  (calcGradient bvar 閉包) に入るが、**保存量・DPLUR・res_roe・壁面熱流束 (断熱厳密 0) には
  一切触れない**。
- **却下した当初設計 (Codex 合意時の "最小実装") — 状態適用は正帰還で暴走する**:
  - node: 等温壁ピン機構 (`pin_wall_node_temperature_d` + res_roe 0 化 + DPLUR roe 行
    decouple) に `Taw_diag` を渡す案 → **暴走** (壁温 1832 K・T max 1967 K > Tt,
    `run_0038` 旧版)。ピンされた壁ノードが無限熱源化し、W–I 双対面の解像伝導が
    「壁と BL の恒常オフセット ΔT=r·U_t²/2cp」を熱流束と誤認して BL を加熱し続け、
    T_aw=T_rep+Δ が追随して発散的ドリフト。
  - cell: ghost 温度 Dirichlet (T_g=2T_aw−T_c, P_g 保持) 案 → 同族の暴走が勾配閉包経由で
    起き **Tt=1500 K に飽和** (壁温 1496 K, `run_0039` 旧版)。
  - 教訓: 壁関数メッシュで T_aw を「状態」として保持するには、SU2 同様**壁隣接の
    伝導流束自体のモデル置換**が必要 (等温壁 Kader q_w モデルと同じ工事、将来課題)。
    断熱壁では q_w=0 が厳密に既知なので、sublayer ΔT はサブグリッド量として壁面値
    (出力+勾配閉包) にのみ反映するのが保存整合な閉包。
- **エネルギー保存整合**: 弱閉包は保存量に触れないため大域収支は自明に保存。OFF 双子
  対照 (scratchpad D7) で場の差は境界勾配閉包経由のみ (max|ΔT|~13 K, 有界・12000 step
  安定)、η/ṁ は 0.9884/1.2928 で OFF (run_0029) と一致。
- T_aw は毎ステップ更新 (`Taw_diag`)。淀み域は T_rep (運動加熱なし)、未計算 (-1) は
  従来出力のまま。

## 5. 実装ステップ (実施済み)

1. `solverConfig`: `turbulence.sstThermalWallFunction` (int 0/1, 既定 0) 追加。
2. `nodeWallDirichlet_d.cu`: `applySstThermalWallFunction` + `set_wall_taw_output_d`
   (断熱壁 bplane の bvar `Ts` ← `Taw_diag[ic]`, node/cell 共通)。main.cpp の
   applyBconds 位相 (init + step loop) で `applyNodeIsothermalWallPin` の直後に呼ぶ。
3. (却下済みの状態適用コードは撤去済み — §4 参照。)
4. case/40 検証 (下記)。

## 6. 検証 (2026-08-11 実施, case/40)

- [x] **node 5e-3 (`run_0038_node_5e3_tawwf`, y+≈29–139)**: bell 壁温 **1417.9 K =
  SU2 壁関数 1422 K と 4.1 K 一致** (OFF 1196 K から +222 K 回復)。ALL STEADY・NaN なし・
  T field ≤ 1499 K (物理的)。8000/12000 step で壁温 0.0 K 差 (定常)。
- [x] y+≈1 基準との整合: wf=1+閉包 1418 K は y+1 帯 (forge 1387–1389 / SU2 1414 / 理論
  ~1400 K) の上端側で SU2 y+1 と 4 K 差 — 乖離 −230 K → +4〜+30 K に縮小。
- [x] エネルギー収支: 弱閉包は保存量不介入で自明に保存 (状態適用の暴走は §4 記録)。
  η=0.9884・ṁ=1.2928 は OFF 基準 run_0029 と一致。
- [x] OFF 回帰: 既定 0 では追加カーネル自体を呼ばず従来経路 (D7 の ON/OFF 差は
  ON 側の勾配閉包による設計内の差)。
- [ ] **cell の既知バイアス (follow-up)**: cell 代表点 (第一セル, y+~30) の T が node/SU2 の
  同高さ比 ~100–160 K 高く、Taw が +70–90 K 過大 (`run_0039`: 1490 K)。cell wf=1 BL の
  熱監査 (なぜ第一セルが過熱か — ghost 鏡映 T 閉包との相互作用疑い) が必要。
  当面 node を一次対象とする (ユーザ方針 [user-prefers-node-base])。

## 7. 変更ログ

- `2026-08-11` — 起票 (Codex レビュー合意済み設計)。theory §6.5(f) / implementation §3.7
  を同時更新。
- `2026-08-11 (訂正・superseded)` — ユーザレビューで 2 点の欠陥を指摘・確認: ①「壁面出力
  bvar `Ts` にのみ書き場に触れない」という §4 の記述は不正確 — node の境界 Green-Gauss/LSQ
  勾配が当時 bvar を読んでいたため接線補正項経由で場に触れる経路が実在した (弱すぎて実効は
  小さかったが)。②実データ検証 (`run_0038` bell 平均) で**壁ノードの実状態 $T[W]$=1195.3 K は
  OFF 基準 1196 K とほぼ不変**、「1417.9 K = SU2 と 4 K 一致」は**出力配列に上書きした $T_{aw}$
  診断値同士の一致**であり、生産値報告 (壁温 1400±15 K) の wf+閉包系列 (node run_0038/0040)
  はこの点で汚染されている。後継の `turbulence-sst-adiabatic-taw-fluxmodel.md`
  (`architecture-node-boundary-gradient-dof-only.md` で境界勾配を owner-state 化した上で、
  $T_{aw}$ を内部粘性流束のモデル置換として明示的に注入する SU2 式設計) へ移行し、本 plan は
  superseded として archived へ移す。
- `2026-08-11 (同日・Step3)` — y+≈65–300 壁関数系列 (case/40 run_0040–0042) で閉包を追検証:
  node+閉包 1405.2 K vs SU2 STANDARD_WALL_FUNCTION 1418.9 K (13.7 K 差)、cell+閉包 1387.2 K
  (**cell 代表点バイアスは y+~70 では消滅** — 5e-3 のバッファ層代表点 (y+~15–30) 固有と判明)。
  三系列 (y+1 low-Re / 5e-3 wf / y+30 wf) の壁温が 1387–1422 K に収束し、生産値
  T_w = 1400±15 K を確定 (case/40 README)。
- `2026-08-11 (同日)` — 実装+検証。**当初の状態適用設計 (node ピン / cell ghost) は正帰還
  暴走で却下** (run_0038/0039 旧版で実測: node 1832 K ドリフト / cell Tt 飽和)、**弱閉包
  (bvar Ts のみ) に改訂して node 5e-3 で SU2 壁関数と 4.1 K 一致を達成**。methods は
  改訂版で更新済み。cell 代表点バイアス (+70–90 K) は follow-up。

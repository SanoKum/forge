# セッションプロンプト: 壁温の真値確定 3 段階 (y+1 基準解 → 熱的壁関数 → y+>30 系列)

以下を新セッションの最初のプロンプトとして使う。

---

case/40 ベルノズル (Pt 4 MPa / Tt 1500 K / 背圧 20 kPa / ε9, SLAU+SST) の**壁温と η_CF を
生産値に格上げする 3 段階**を、上から順に実施してほしい。Codex レビューで合意済みの計画である。

基準文書 (まず読む):
- [`notes/sessions/wall-temperature-three-way-analysis.md`](../../notes/sessions/wall-temperature-three-way-analysis.md)
  — 壁温問題の全記録。**§7 の訂正まで読むこと** (§3 の初期推察は棄却済み)。
- [`plans/active/boundary-node-nozzle-wall-outlet-stability.md`](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md) §2.8
- [`case/40.nozzle_design_tool/README.md`](../../case/40.nozzle_design_tool/README.md) の run 一覧

## 背景 (確定済み事実 — 再調査不要)

1. **node 断熱壁の伝導リークは修正済み** (壁温 1631→1196 K、市松 206→44 K; `run_0029_node_adiabfix`)。
2. 残る壁温 −230 K (forge 1193 K vs SU2 壁関数 1422 K) の**真因は forge SST 壁関数
   (`wallTreatmentSST=1`) に断熱回復温度の熱的閉包 (Crocco 型) がないこと**。
   `Taw_diag` 診断 (実装済み, `ransWallFunction_d.cu`) は T_aw = T_rep + Pr^{1/3}U_t²/(2cp) を
   出力し、ベル部 mean **1420.2 K = SU2 壁関数 1422 K と 2 K 差**で一致することを確認済み。
   BL 内部状態は壁関数同士なら forge≈SU2 (第一内点 T ~1000 K, μt/μ 4–7, ṁ 一致)。
3. η_CF は壁処理依存で 0.955 (low-Re 系) 〜 0.988 (壁関数系) に割れており、**真値は
   y+≈1 low-Re 基準と y+>30 壁関数系列で挟み込むまで未決**。現メッシュ (壁 frac 5e-3,
   y+≈19–40) はどちらの手法にも不適な帯。
4. `thermCondMethod: 1` (constant-Pr: k=μ(T)cp/`prandtlLam`) 実装済み (既定 0)。SU2 は
   CONSTANT_PRANDTL 0.72 なので、**SU2 と揃える比較では forge も `thermCondMethod: 1,
   prandtlLam: 0.72` を使うこと** (従来の固定 k は分子 Pr が 1.5–1.9 に漂う)。
5. 軸対称は `nodeAxisDirichlet: 1` (現 runner 既定)。単成分 CPG 限定ガードあり。
6. `check_quasisteady.py` は step パース修正済み (`^res_(\d+)\.h5$` のみ)。

## Step 1: y+≈1 low-Re 基準解 (forge cell / forge node / SU2 の三者)

- **メッシュ**: 現行 5e-3 で y+≈19–40 → y+≈1 には `wall_first_frac` ≈ **1.5e-4〜2.5e-4** が目安
  (走らせて `ypls` 出力 or 第一セル距離×u_τ で実測確認し調整)。`problem_bell_smoke.yaml` を
  コピーして専用 yaml を作る (`problem_bell_smoke_yp1.yaml` 等)。
  **必ず `check_mesh_quality.py` で AR≤1000 を確認** (近壁細分で AR が跳ねたら nj 増や
  接線方向も細分して調整。品質 FAIL のまま投入しない)。
- **forge 設定**: `wallTreatmentSST: 0` (low-Re)、`thermCondMethod: 1, prandtlLam: 0.72`。
  cell と node の両方 (ユーザ方針: 軽いケースは両系統確認)。
- **既知の罠 (node × 細壁)**:
  - 壁 frac 1e-3 は node で発散した既往があるが、真因は cross-mesh IC の壁法線 P 段差で
    **原始変数補間で根治済み** (memory: node-isothermal-wall-thin-cell-mass-source,
    「ソルバは y1+=1 等温壁で安定」)。`interp_field.py` の補間モードを確認して使うこと。
  - warm start 必須。cfl4 直投入で初手発散したら**段階起動** (soft explicit cfl0.1 1次
    3000 step → 陰解法 2次; run_0019/0024 の手順)。
  - node wf=0 は case/40 で ṁ が 4.5% 低下した未解明挙動あり (分析文書 §4)。y+1 で
    再評価し、cell/SU2 と ṁ・壁温・Cf を突き合わせる。
- **SU2 側**: `run_0028_su2_sst` の cfg をベースに同一メッシュ (gmsh で .su2 化) で low-Re。
  収束は固定 CFL 継続まで見る (`procedures/su2-cross-check.md` の CFL_ADAPT 注意)。
- **判定**: 三者の壁温 T_w(x)・Cf・ṁ・η_CF を突き合わせ。y+≈1 で SST low-Re は正当なので、
  ここでの一致が「真値」の基準になる。壁温は解析回復温度 T_aw≈1390–1420 K と比較。
- run は `run_0032_` から連番 (**run_0030 が 2 個ある衝突あり** — README の備考参照。
  以後は必ず既存 run を ls してから採番)。

## Step 2: SST 壁関数の熱的閉包 (compressible adiabatic thermal wall-function) の設計実装

**開発フロー必須**: `methods/turbulence/implementation.md` (§6.5 壁関数) を先に更新し、
`plans/active/turbulence-sst-thermal-wall-function.md` を起票してから実装する。

- **設計方針** (Codex 合意): Crocco 型 T_aw を壁状態に適用する。
  - node: 既存の等温壁ピン機構 (`pin_wall_node_temperature_d` + `iso_wall_flag` の roe 行
    decouple + res_roe ゼロ化) を**再利用**し、Tsb に固定値でなく毎ステップの `Taw_diag`
    (代表点から計算) を渡すのが最小実装。適用は `wallTreatmentSST==1 && 断熱壁` のみ。
  - cell: 壁 ghost/bvar の Ts_b を鏡映値でなく T_aw にする (壁面熱流束は断熱=0 のまま。
    ghost T の変更は勾配閉包と出力壁温に効く)。
  - **エネルギー保存との整合を plan に明記**: T_aw ピンは壁ノードのエネルギー式を
    Dirichlet 化する。大域エンタルピー収支 (入口 vs 出口の Ht 流束) を検証項目に入れ、
    ピンが偽のエネルギー源になっていないことを確認する。
  - 検証: 現行 5e-3 メッシュで壁温が SU2 壁関数 (1422 K) と数 K で一致、Step 1 の y+1 基準
    と比較して壁温・Cf の乖離が縮むこと。既定 ON にするかは結果を見て判断 (まず opt-in
    config `sstThermalWallFunction` 等)。
- 参考実装: SU2 `CNSSolver::SetTau_Wall_WF` 周辺の T_aw 更新 (`.external/su2-src`)、
  forge WMLES の Kader 温度壁法則 (`wmlesWallModel_d.cu`)。

## Step 3: y+>30 壁関数系列で SU2 と再比較

- `wall_first_frac` を 1e-2 前後に粗くした系列 (y+≈40–80) を作り、forge wf=1+熱的閉包
  (Step 2 実装後) vs SU2 STANDARD_WALL_FUNCTION (`run_0030_su2_sst_wallfunc` の cfg 流用) を
  同一メッシュで比較。
- Step 1 (y+1 基準) と合わせて **η_CF の真値を挟み込み**、壁温・η を「生産値」として
  README に確定させる。格子 3 点あれば Richardson/GCI も添える。

## ルール (AGENTS 準拠の要点)

- 実行は `solver_density_cuda/tools/run_case.sh <絶対パス>`。struct 変更後は **full rebuild**
  ([stale-build-struct-layout-trap])。
- 「収束/一致」報告には `check_convergence.py` の VERDICT 必須、派生量 (壁温・η・Cf) には
  `check_quasisteady.py` の VERDICT 必須。プラトーは「準定常」と書き「収束」と書かない。
- run 作成ごとに case/40 README の run 一覧を同期。破棄 run も所在を明記。
- 検証済みマイルストーンごとに commit し `feature/median-dual-3d` へ push。
- 誤りが見つかったら記録 (plan/memory) を先に訂正してから進む。

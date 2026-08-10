# セッションプロンプト: node 出口列欠陥の根治 (最優先) と後続 3 課題

以下を新セッションの最初のプロンプトとして使う。

---

case/40 の壁温真値キャンペーン (3 段階) は完了し、生産値は確定済み (壁温 1400±15 K /
η_CF 0.978±0.003, `case/40.nozzle_design_tool/README.md` の「生産値確定」節が正本)。
本セッションは残課題を優先順に進めてほしい。**最優先は①**。

基準文書 (まず読む):

- [`plans/active/boundary-node-nozzle-wall-outlet-stability.md`](../../plans/active/boundary-node-nozzle-wall-outlet-stability.md)
  §2.9–2.10 — y+1 キャンペーンで修正済みのバグ 2 件と、**未修正の node 出口列欠陥の切り分け状況**。
- [`plans/accepted/turbulence-sst-thermal-wall-function.md`](../../plans/accepted/turbulence-sst-thermal-wall-function.md)
  — 熱的閉包の確定設計 (弱閉包)。**§4「状態適用は正帰還で暴走」は再発防止の必読事項**。
- [`case/40.nozzle_design_tool/README.md`](../../case/40.nozzle_design_tool/README.md) run 一覧 (run_0032–0042 が今回分)。

## 確定済み事実 (再調査不要)

1. **修正済み** (commit 3cdddb9e): 双対 CV 重心/体積の float32 シューレース桁落ち
   (`gmshReader.hpp` 2D, A 相対+double 化) と `interp_field.py` の node 照会点
   (`MESH/COORD` 使用)。node y+≈1 low-Re はこれで安定 (`run_0034`)。
2. **実装済み** (commit 210ce982): `sstThermalWallFunction: 1` = Crocco T_aw の**弱閉包**
   (壁面値 bvar Ts のみ、保存量不介入)。node 壁温 = SU2 壁関数と 4–14 K 一致。
   **T_aw を状態 (node ピン / cell ghost) として課すと正帰還で暴走する** (実測済み) —
   強閉包には壁隣接伝導流束のモデル置換が必要 (課題③)。
3. **生産構成キー** (case/40): `wallTreatmentSST: 1` + `sstThermalWallFunction: 1` +
   `thermCondMethod: 1, prandtlLam: 0.72` + `turbulentPrandtl: 0.9`
   (**Prt 0.85→0.9 で壁温 +20.5 K の強感度**、SU2 既定 0.9 に整合させた)。
4. **未修正**: node 出口列欠陥 (課題①)。node の η 出口積分値は +1.3% 過大であり、
   **修正まで node の η は内部列積分か cell/SU2 を正とする** (README に注記済み)。
5. y+≈1 node の起動レシピ: warm は **node 場から** (cell 場は x 半間隔ずれで近壁 staircase)、
   段階起動 = 陰解法 cfl_pseudo 0.5・convMethod 0・nStepInner 10 で 3000 step →
   陰解法 cfl4・全域 2 次 (explicit は壁 ω~1e10–11 の剛性で不可)。
6. run 採番は **run_0043_ から** (README の表を ls と突き合わせてから採番)。

## ① node 出口列欠陥の根治 (最優先 — η_CF の生産品質に直結)

**現象** (plan §2.10 に定量記録済み): node の出口境界ノード行が内部 2 列の外挿から
P 平均 −19% (最大 −24%)・T −15% ずれる。**超音速コア (M≈4.1) でも 1 列で P −12%** =
背圧情報が届かないはずの行なので純数値欠陥。cell は ~1%、SU2 (同じ頂点中心
median-dual) は −0.2% で清浄。`thrust_metrics` が outlet iCells (=この列) で積分するため
node η が +1.3% 過大になり、**旧「node−cell η 差 = 2 次化の離散差」帰属はこの
アーティファクトだったと訂正済み**。

**切り分け済み** (再実験不要):

- 全域 1 次でも近壁行の sag −20% は残存、コア行は −12%→−4% に減少 → 2 次再構成
  (出口ノードの片側勾配×limiter) は増幅要因であって種ではない。
- 出口半割面の対流流束は `convectiveFlux_boundary_d` で F(bvar)、`outlet_statPress` の
  bvar は超音速分岐で node 自身の状態 = F(W) 純風上なので対流は自己整合。
- **残容疑 2 つ**: (a) **出口半割面の粘性流束の扱い** (sag が壁側ほど大きいことと整合 —
  W–I 面や境界半割面でエネルギー/運動量の粘性項がどう閉じているか)、(b) W–I 主ループ面の
  2 次再構成の境界側閉包 (増幅分)。

**進め方の推奨**: 収束場 (`run_0037_node_yp1_prt09/res_12000.h5` か
`run_0040_node_yp30_tawwf`) を読み、**出口ノード 1 個 (超音速コア行) の半割 CV の
離散収支を項別に再構成**する (対流/粘性/軸対称ソースを面ごとに手計算 or デバッグ出力) —
どの項が O(12%) の不均衡を支えているかを直接見るのが最短。修正後は
run_0043_ 系で node y+1 / y+30 を再取得し、出口列の外挿一致 (~1% 以下) と
η 出口積分 = 内部列積分の一致を判定基準にする。η の生産値 (0.978±0.003) が
node 出口積分でも再現されたら README の注記を解除する。

## ② SST エネルギー壁関数 (Kader q_w 流束置換) — plan 起票済み

[`plans/active/turbulence-sst-thermal-flux-model.md`](../../plans/active/turbulence-sst-thermal-flux-model.md)
に設計済み (draft)。本命は**等温壁×粗メッシュの熱負荷 q_w** (冷却ノズル壁の生産量)。
必要部品は全部ある: Kader T⁺ は `wmlesWallModel_d.cu`、流束置換は `Qw_Wall` マーカの
xor 機構 (`viscousFlux_d`)、u_τ は SST Reichardt Newton。新設は「SST から q_w を計算し
Qw_Wall に書くカーネル」のみ。検証は平板等温壁 y+ 掃引 (相関式+SU2) → case/40 等温壁変種。
着手時は plan の status を in_progress にし、methods §6.5(g) を先に書く (開発フロー)。

## ③ 3D median-dual 幾何の桁落ち除去 — plan 起票済み

[`plans/active/architecture-median-dual-3d-double-geometry.md`](../../plans/active/architecture-median-dual-3d-double-geometry.md)
に現状分析済み (draft)。3D は tetVol/重心は既に安全で、**露出は `newell()` 面積ベクトル
(絶対座標の差×和) と境界半割面蓄積のみ**。修正は 2D と同じローカル原点化+double 蓄積。
3D 薄壁 (y+~1) node を回す前に必ず実施。検証手順は plan §6 に列挙済み。

## ④ cell 5e-3 の wf=1 BL 熱監査 (低優先)

cell の熱的閉包は 5e-3 メッシュ (代表点=第一セル中心が y+~15–30 のバッファ層) でのみ
T_aw が +70–90 K 過大 (`run_0039`)。y+~70 では消滅 (`run_0041`) と確認済みなので実用上の
支障は小さいが、cell 第一セルが node/SU2 の同高さより ~100–160 K 熱い理由 (ghost 鏡映 T
閉包との相互作用疑い) は未解明。accepted plan の checklist 残項目。

## ルール (AGENTS 準拠の要点)

- 実行は `solver_density_cuda/tools/run_case.sh <絶対パス>`。`solverConfig.hpp` など
  struct 変更後は **full rebuild** (`ninja -t clean && ninja forge convertGmshToForge`)。
- 「収束/一致」報告には `check_convergence.py`、派生量には `check_quasisteady.py` の
  VERDICT 必須。プラトーは「準定常」と書く。
- run 作成/破棄ごとに case/40 README の run 一覧を同期。新 run 作成時は旧 res_*/log を
  必ず rm ([run-dir-copy-clean-res])。
- 検証済みマイルストーンごとに commit し `feature/median-dual-3d` へ push。
- 誤りが見つかったら記録 (plan/memory) を先に訂正してから進む。
- `bndFirstOrder` は正当な手段とみなさない (診断にも使う場合は明示し、レシピ・結論の
  前提にしない — ユーザ指示 2026-08-11)。

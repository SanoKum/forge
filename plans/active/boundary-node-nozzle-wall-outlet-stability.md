# node (median-dual) 2D の傾斜壁・出口コーナー不安定 (ノズル系で初手〜百 step 発散)

## メタ

- **area**: `boundary / discretization`
- **status**: `draft`  <!-- 調査完了・修正未着手 -->
- **related_docs**:
  - `methods/discretization.md` / `methods/boundary.md`
- **related_plans**:
  - [`tooling-nozzle-phase0-foundation.md`](tooling-nozzle-phase0-foundation.md) (発見元: node での再計算が全滅)
  - [`../accepted/turbulence-node-inlet-dirichlet-conserved.md`](../accepted/turbulence-node-inlet-dirichlet-conserved.md) (入口側の類似修正 — 出口版が欠けている疑い)
  - [`../accepted/diffusion-node-boundary-real-distance.md`](../accepted/diffusion-node-boundary-real-distance.md) / [`../accepted/turbulence-node-sst-wallfunction.md`](../accepted/turbulence-node-sst-wallfunction.md)
- **created**: `2026-08-03`
- **owner**: `sano` (調査: Claude 自走)

## 1. 目的 (現象)

forge_design の③ベルノズルメッシュ (構造化 quad 221×65、壁クラスタ AR≦19、品質 PASS、
**cell モードでは同一 msh が完走**) を node で回すと、あらゆる設定で序盤発散する。
2D node の検証実績 (平板 case/26・backstep case/18・case/36) は全て**軸平行壁**であり、
傾斜壁 + 超音速沿い流れ + 超音速流出コーナーを持つノズル系は node 初適用で未踏だった。

## 2. 調査結果 (2026-08-03, 全て 400 step 以下の再現 run at scratchpad)

### 切り分けマトリクス (単独では全て無効)

| 変更 | 結果 |
| --- | --- |
| SST → laminar / convMethod 1→0 / cfl 4→0.5 / axisym→平面 | いずれも step 2–4 で NaN |
| outlet_statPress→outflow / wall→slip / inlet_Pressure→速度入口 | 同上 (残差爆発値がビット同一 = 種は境界種別に非依存) |
| bndFirstOrder 0→1 / implicitRelax 1.0→0.3 / 壁ノード IC 速度ゼロ化 | 同上 |
| 双対幾何検査 (体積・面積・法線外向き) | 異常なし |

### 段階起動 (stage A: convMethod 0 + cfl_pseudo 3 + implicitRelax 0.5) で層が分離

| 構成 | 発散 step | 種の位置 |
| --- | --- | --- |
| 既定設定 (2次, relax 1, cfl 4) | 2–4 | 傾斜壁 (収縮部/ベル) の壁・壁隣接ノード。step 1 で壁 slip ノードの ΔroUx ~ O(10³) |
| stage A + SST (wf 0/1 とも) | 11–15 | **出口∩壁コーナーの outlet ノードで roOmega → 1e19–1e20** (流れ残差は健全 1.7e-4) |
| stage A + laminar | 128 | ベル中腹の壁隣接帯 (x≈47mm)。それまで残差は順調に低下 (1e-4) |
| 参考: 準直線壁 (勾配≦5°) 対照 | 5–6 (SST) | 壁隣接。傾斜が緩いと遅延 → 勾配依存性あり |

### リグレッション対照 (2026-08-03 追記)

「以前は node が動いていた」との整合確認: 検証済み node ケース case/26 平板 SST
(run_node_sst_final の入力一式、turbulence ブロックのみ新 API に翻訳) を**現行バイナリ**で
500 step 再実行 → **NaN なし・残差順調低下 (1.6e-3→2.5e-4)**。すなわち**ブランチの
リグレッションではなく**、軸平行壁で検証された node が、ノズル構成 (傾斜壁・超音速流出
コーナー) という未検証領域で破綻している。なお既存 node run の solverConfig は旧
turbulence API のため現行バイナリではそのまま再実行できない (要翻訳) 点にも留意。

### 結論 (二層の独立した node 問題)

1. **壁境界の運動量不安定 (傾斜壁×沿い超音速流)**: 設定を柔らかくすると遅延するが根治しない
   (2→15→128 step)。壁隣接ノードで成長するモードで、slip でも no-slip でも発生。既知の
   未修正バグ「node slip + 接線密度勾配の市松スプリアス流」(node-slip-spurious-flow) と
   同族の可能性。壁半割面流束の傾斜壁での取り扱い (弱形式置換・実距離 over-relax の適用域)
   を疑う。
2. **出口境界の SST スカラー爆発 (コーナー)**: roK/roOmega が出口∩壁コーナーから成長。
   入口側は [turbulence-node-inlet-dirichlet-conserved] で「Dirichlet 保存量整合+残差除外」
   が入っているが、**出口側に同等の処置が無い**。

## 2.5 追補 (2026-08-03 後半): 実証レシピ発見と原因の確定 — §2 の初期仮説を大幅更新

ユーザ指摘「以前ノズルでも node をやっていた」に基づき **case/29.bell_vs_conical の node
キャンペーン (run_0019〜0041: 層流壁レシピ・CFL 掃引・implicit 掃引・SST 実走)** を発見。
§2 の「ノズル系は node 初適用」は**誤り**だった。以後の対照実験で真因を確定:

- **リグレッション対照 2**: case/29 `run_0038_node_sst_wc` (node SST ノズル) を現行バイナリで
  verbatim 再実行 (turbulence ブロックのみ新 API 翻訳) → 3000 step 健全 (rms_ro 7.9e-5→1.45e-6)。
  **リグレッションではない**。
- **真因 = 壁第一セルの細かさ**: 同一レシピで `wall_first_frac` 1e-3 (第一セル 10–30μm,
  y+~15–50) は全設定で発散、**5e-3 (50–150μm) で完全収束** —
  `case/40/run_0006_bell_node_warm`: **VERDICT PASS / ALL PASS** (全列 3.0–4.6 桁低下,
  η_CF=0.9809 vs cell 0.9790)。case/29 の壁寄せ (Bump0.4) メッシュが粗めだったことと整合。
- **実証レシピ (node ノズル SST)**: explicit RK3 + cfl 0.1 + convMethod 0 (1 次) +
  `bndFirstOrder/nodeWallDirichlet/axisCentroidShift` + katoLaunder + **収束場からの
  warm start 必須** (冷間 IC は粗壁でも発散)。forge_design runner に既定として実装済み。
- なお §2 の「出口コーナー/傾斜壁」は症状の**発現場所**であり、独立バグとしての扱いは
  下記課題の解決後に再評価する。

### 残る解決すべきソルバ課題 (本 plan の実体・優先順)

1. **node × SST × 陰解法 (blockDPLUR)**: 収束済み node 場からの引き継ぎでも step 4 で発散
   (`run_0007_bell_node_2nd_imp` 相当)。case/29 でも implicit 成功は層流のみ (run_0035 cfl2–5
   で rms_ro 2e-9)。→ 最適化ループの評価速度に必須 (explicit cfl0.1 は ~10 倍遅い)。
2. **node × 2 次精度 (convMethod 1)**: 収束場からでも step 7–10 で発散 (cfl 0.3/0.15 とも)。
   case/29 の node SST も全て 1 次だった。→ 生産精度に必須。
3. **壁第一セル細分の頑健性**: y+~1 (low-Re) 相当の壁間隔が node で回らない。→ q_peak
   壁解像 (親計画 §4.6 ①) に将来必要。

## 3. 修正方針 (次セッションの作業項目)

1. 再現系は確立済み: `case/40` の問題 YAML + runner で 3–5 分/run。切り分け操作は
   `wall_first_frac` (1e-3 ⟷ 5e-3)・`timeIntegration` (3 ⟷ 11)・`convMethod` (0 ⟷ 1) の
   3 ノブで decisive に再現する。
2. 課題 1 (SST 陰解法): blockDPLUR の k/ω 行 (壁 Dirichlet/壁関数ノード・出口ノード) の
   Jacobian 整合を点検。scalar-DPLUR 側 ([time_integration-scalar-dplur] 系) と壁代表点修正
   ([turbulence-node-sst-wallfunction]) の併用パスが濃厚。
3. 課題 2 (2 次精度): 壁近傍ノードの勾配再構成 (Green-Gauss / [discretization-lsq-gradient])
   と `bndFirstOrder` の適用範囲を点検。発散種は常に壁隣接第一内点。
4. 課題 3 (細壁間隔): 課題 1–2 解決後に wall_first_frac を段階細分して再評価。
5. 検証: 本ノズルメッシュで node SST **2 次 + 陰解法** 12000 step 完走 + cell と η_CF 一致 ~1%。

## 4. 影響・当面の運用

- **node は実証レシピ (1 次・explicit cfl0.1・warm start・粗め壁) で運用可能** — forge_design
  runner の node 既定に実装済み。2 次+陰解法が要る生産評価 (Phase 2 最適化) は課題 1–2 の
  解決までは cell を併用。
- 関連 run: `case/40.nozzle_design_tool/` README の run 一覧 (run_0005 発散記録 /
  run_0006 node 収束 PASS / run_0007 2次・陰解法の発散記録)。
- `interp_field.py` を node 対応に改良 (入力 h5 は CELLS/centCoords 優先、node res は
  ノード座標直用)。

## 5. 変更ログ

- `2026-08-03` — 起票。調査マトリクス完了 (§2)。
- `2026-08-03` — ユーザ指摘で case/29 node ノズルキャンペーンを発見し §2.5 で原因確定:
  リグレッションでも未踏領域でもなく、**壁第一セル過細 + レシピ不一致 (warm start /
  explicit / 1 次)**。実証レシピで node SST ノズル初の VERDICT PASS (run_0006)。残課題を
  「SST 陰解法 / 2 次精度 / 細壁間隔」の 3 点に再定義。

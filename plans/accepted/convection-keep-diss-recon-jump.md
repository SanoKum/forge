# KEEP ES 散逸の再構成ジャンプ化 (解像スケール保護)

## メタ

- **area**: `convection`
- **status**: `done`
- **related_docs**:
  - `methods/convection/implementation.md` (KEEP_d の散逸レイヤ節)
- **related_plans**:
  - `plans/accepted/convection-keep-es-dissipation.md` (親: 散逸レイヤ本体)
- **created**: `2026-07-19`
- **owner**: CFD Dev

## 1. 目的

matrix ES 散逸 (keepDissType=2) は隣接対の**生ジャンプ** Δw に働く 2 次精度型のため、
よく解像された滑らかな構造からも常に KE を吸う (64³ TGV L3 で層流期 −0.6%・終値 −1.5%、
遷移期の ε\* カーブが DNS より早く立つ)。ジャンプを**再構成後ジャンプ**
$\Delta_{rec} = \Delta_{raw} - \tfrac12(g_i+g_j)\cdot d_{ij}$ に置換し、市松減衰能力を
無傷のまま解像スケールへの散逸コストを桁で下げる。

## 2. スコープ

- **やる**: `keepDissJump` config (0=生ジャンプ・既定・ビット不変 / 1=再構成ジャンプ)。
  対象は **matrix CPG 枝 (keepDissType=2, thermalMethod≠2)**。node を一次検証対象とする
  (ユーザ方針: node ベース主体)。cell も同一コードパスで有効。
- **やらない (フォローアップ)**: scalar (type1)・TP matrix 枝への展開 / リミタ付き再構成
  (LES 用途は滑らか場前提でリミタ無し κ=0 線形) / 厳密 ES 証明付き再構成
  (Fjordholm sign-property 型) — v1 は L1/L2/L3 の実測ゲートで担保する。

## 3. 関連 docs と前提

- 散逸レイヤ本体: [convection-keep-es-dissipation](../accepted/convection-keep-es-dissipation.md)
  (σ=0.02 を 64³+実 DNS で確定済)。
- 勾配は既存の `calcGradient_d_wrapper` (KEEP でも毎ステップ計算済) + node periodic 勾配 gather
  (case/35 run_0011/0012 で検証済) を使用。追加コストは KEEP_d での勾配読みのみ。

## 4. 設計方針

**鍵となる性質 (実データで事前検証済)**: 純 2Δ モード (市松) の中心勾配は厳密ゼロ
→ 再構成はジャンプを全く減らせない (フル振幅のまま散逸される)。滑らかな場のジャンプは
O(h)→O(h³)。**閾値パラメータ不要・線形 (指数減衰維持)・2Δ 選択的**。

- 事前測定 (node 実 run の VALUE/dPd\* 使用, 2026-07-19):
  散逸ドレイン比 $\Sigma\Delta_{rec}^2/\Sigma\Delta_{raw}^2$ = **市松 1.003 / TGV 層流期 0.006 /
  遷移期 0.16 / ピーク 0.45** — 市松無傷・滑らか場のコスト 2〜170 分の 1。
- 参考: Jameson 型振幅センサ案は棄却済み — 発達乱流の 2Δ 帯振幅 (p90 1.9e-3) が市松 (1e-3)
  を**上回り**分離不能、かつ振幅比例ゲートは市松減衰を指数→代数に落とし L1 破綻。
- 実装: matrix 枝の Δw を再構成 primitives
  ($\rho_{L/R} = \rho_{i/j} \pm \tfrac12 g_\rho \cdot d_{ij}$ 等, κ=0 線形・リミタ無し) から
  組み直す。面平均量 (roF, c, Ht 等) は従来どおりセル平均。
- **再構成不可の面は生ジャンプへフォールバック**: cell の ghost 面 (`ip >= nNormalPlanes`,
  ghost の勾配/重心が無効) と、再構成 ρ/p が非正になる面 (リミタ無しの正値性ガード)。
  ~~host の per-plane マスク~~ は実装時に不要と判明: node の主ループは内部双対面のみで、
  slave CV の重心はノード自身の座標・勾配は gather 済のため $d_{ij}$ も勾配も健全 (§変更ログ)。
- ES 性: 生ジャンプの ES 証明は再構成ジャンプには直接は延びない (厳密には sign-property
  再構成が必要)。v1 は「散逸形 (SPD 二次形式) を保つ + L1/L2/L3 実測」で担保し、既定 off。

## 5. 実装ステップ

1. config: `keepDissJump` (optional, 既定 0)。
2. wrapper: reconOK マスク構築 (初回のみ)・勾配ポインタと共に KEEP_d へ。
3. KEEP_d matrix CPG 枝: `keepDissJump==1 && reconOK[ip]` で Δw を再構成状態から計算。
4. full rebuild ([stale-build-struct-layout-trap])。
5. 検証 (node 一次):
   - L1 node (case/35): 市松 + matrix σ=0.05 + jump=1 → run_0024 (raw) と同等の減衰
     (共分散射影, plot_checkerboard_decay.py)。
   - L3 node 64³ (case/09): node σ=0 基準 run + σ=0.02 jump=1 → K/K0 カーブが σ=0 側へ寄る
     (予測: 層流期コスト ~0、終値ギャップ半減以上)。
   - cell: 既定 off のビット不変 + cell L1 スモーク。

## 変更ログ

- 2026-07-19: plan 作成。振幅センサ案の棄却根拠と再構成ジャンプの事前測定を記録。
- 2026-07-19: **実装・検証完了 (node 一次)**。
  - 実装: `keepDissJump` (既定 0=ビット不変) + `KEEP_d` matrix CPG 枝で Δw を再構成 primitives
    (面重心へ κ=0 線形, `FaceGeom.pcx/ccx` + `GradFields`) から構成。ghost 面 (`ip>=nNormalPlanes`)
    と再構成正値性 NG の面は生ジャンプへフォールバック (host マスク配列は不要と判明 —
    node 主ループは内部双対面のみで勾配 gather 済・d_ij 健全)。
  - **L1 (市松, case/35)**: node jump=1 は A_cb 1e-3→6.73e-6/400step = 生 (6.78e-6) と同一。
    cell も 6.29e-8 = 生 (7.09e-8) 同等以上。**市松減衰は無傷**。既定 off 回帰は 4 桁一致 (run_0041)。
  - **L3 (64³ TGV Re=1600+WALE, case/09 run_0039/0040)**: σ=0.02 の KE コストが
    **終値 −2.1%→−0.8%、層流期 −0.7%→+0.7% (σ=0 と同一)**、遷移期 ε\* の過剰立ち上がり解消、
    ピーク時刻 8.96 維持 (σ=0 単独は 9.24 と遅い)。**市松頑健性 × σ=0 並み KE 追従の両立を達成**。
  - **node LES 推奨構成を更新: `keepDissType: 2, keepDissCoeff: 0.02, keepDissJump: 1`**。
  - 残 (フォローアップ): scalar 枝・TP matrix 枝への展開 / 衝撃を含む流れでのリミタ付き再構成 /
    厳密 ES (sign-property) 化。

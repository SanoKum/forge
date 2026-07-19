# 引き継ぎプロンプト: KEEP-LES 基盤完成後の展開 (2026-07-19 セッション末更新)

> **⚠ 後継あり**: 本ファイルの続きは [iddes-phase2-session-prompt.md](iddes-phase2-session-prompt.md)
> (2026-07-19 後段セッション末)。候補1 (advGauge node)・候補2 (backstep パイロット) は完了し、
> 前処理散逸・単一スキーム DES 基盤まで進んだ。新セッションは後継を読むこと。

> 新セッションの冒頭にこのファイルを読ませる。KEEP-LES の数値基盤 (散逸設計・free-stream 保存・
> SGS 整理) が実測+DNS 裏付きで完成した。ここから先の候補タスクと、罠・運用ルールをまとめる。

## 完了済み (再調査・再実装しないこと)

1. **ES 散逸レイヤ完成形**: `keepDissType:2` (matrix) + **`keepDissJump:2`** (再構成ジャンプ
   + sign-property クリップ = **証明付き ES**)。市松減衰無傷・解像スケール保護
   (64³ TGV K/K0(10) −1.4%)。正本: [convection-keep-diss-recon-jump](../../plans/accepted/convection-keep-diss-recon-jump.md)。
   振幅センサ案は棄却済 (乱流 2Δ 帯と市松は振幅で分離不能) — 再挑戦しない。
2. **σ=0.02 を実 DNS で確定**: 参照は `case/09.Taylor-Green/ref_dns/TGV_Re1600.dat`
   (Dairay+2017 512³, 要引用; K/K0(10)=0.596, ε*ピーク 0.1029@t*=8.98)。旧近似帯 0.50-0.57 は
   誤りだった。σ は物理較正でなく「L1 を満たす最小値」。node=cell 同等検証済。
3. **移流基準差分ゲージ** (`space.roRef`+`uRef`): U∞≠0 一様流の float32 桁落ち根治
   (エネルギー 6.5 桁改善)。**スコープ: KEEP CPG × cell/node 両対応** (node 境界半割面ゲージは
   2026-07-19 後段に実装済 = 旧候補1 完了。case/33 run_0017-0020 で cell と同水準の機械精度を
   検証)。[freestream-preserving plan §8/§8.5](../../plans/active/convection-freestream-preserving-flux.md)。
4. **WALE 2 バグ修正**: ①壁なしメッシュで wall_dist≡0→不活性 (読込後 1e30 ガードで修正)
   ②Sd テンソルが成分2乗の誤式→行列2乗へ。**過去の TGV「WALE 併用」は全部 ILES だった**
   (結果は ILES として有効)。[turbulence-wale-fix](../../plans/accepted/turbulence-wale-fix.md)。
5. **σ-model 追加** (`LESmodel:2`, Nicoud 2011): 実装・検証済だが 64³ TGV では WALE より散逸的
   (−5.5%)。**解像/遷移流では静的 SGS はどちらも逆効果 = ILES+ES が最良**。σ-model の本領
   (壁乱流・未解像高Re) 評価は未実施。[turbulence-sigma-model](../../plans/accepted/turbulence-sigma-model.md)。
6. **定常局所 dt の音響過渡不安定**: unsteady:0 は閉じた周期系の波でスキーム/メッシュ/CFL 不問で
   数 step 発散 (CFL 下げても同 step = 指紋)。**過渡ケースは unsteady:1 必須**。
   case/33 の旧「発散」報告はこれとの交絡 (訂正済)。
7. (前セッションから継続) ES 散逸の CPG/TP/多成分対応・多成分 matrix プラトー根治 (Y 正規化)・
   pRef の KEEP 展開は accepted 済: [convection-keep-es-dissipation](../../plans/accepted/convection-keep-es-dissipation.md)。

## KEEP-LES 推奨設定 (確定版)

```yaml
solver: "KEEP"
keepDissType: 2          # ★ keepDiss* は top-level キー! space 配下は黙って無視される
keepDissCoeff: 0.02      # 市松頑健性優先なら 0.05
keepDissJump: 2          # 証明付き ES・市松無傷・渦保護 (1=証明なし最軽量)
keepDissPrecond: 1       # Turkel 前処理音響散逸 (2026-07-19 後段): 真の市松を c' 版の 142 分の 1 に
                         # 減衰しつつ物理コストは c' 以下。低マッハ壁付き LES の推奨。cprime を上書き
space: {pRef: <動作静圧>, roRef: <動作密度>, uRef: [<平均流速>,0,0]}
  # 非直交+平均流で必須。roRef/uRef は cell/node 両対応 (KEEP CPG のみ; SLAU/TP は off のまま)
turbulence: {LESorRANS: 0, LESmodel: 0}   # 解像/遷移流は ILES (静的 SGS は逆効果と実測済)
time: {unsteady: 1}                        # 過渡は必ず物理 dt
```

## 次の候補タスク (優先度順の提案)

1. ~~**advective gauge の node 対応**~~ **済 (2026-07-19 後段)**: 境界半割面
   (`convectiveFlux_boundary_d`) に upwind 差分形ゲージを実装し fail-fast 解除
   (plan §8.5, case/33 run_0017-0020)。「非直交メッシュ+平均流の node LES」解禁。
   node メッシュは `case/33/mesh/make_wavy_node_msh.py` で wavy.h5 から復元した
   `wavy_node.msh` を使う (旧 .msh 散逸のため)。
2. ~~**実形状 LES パイロット**~~ **完了 (2026-07-19 後段, case/18 run_0103-0124)**。収穫:
   ①σ/c' 掃引で「一様固有値スケーリングは単一トレードオフ曲線に縮退」を実測 →
   ②**Turkel 前処理音響散逸 `keepDissPrecond` を実装・検証** ([plan](../../plans/active/convection-keep-diss-lowmach-precond.md),
   真の市松 142×減衰・物理コスト c' 以下・ES 維持・cell/node 両検証・commit 56a2ede)。
   ③backstep の残存縞 (~22 Pa) は数値市松でなくリップ Δx=0.13 の**解像限界物理**と確定。
   ④SGS 3 者 (ILES 0122 / WALE 0123 / σ 0124, 共通 3D 乱流 seed) が **cfl1.0 で全て安定**・
   vis_turb 正常 (σ は WALE の 1.7 倍散逸的 = TGV 傾向と整合)。
   **定量 x_R 比較はスコープ外と judgment** (入口一様 Pt/Tt = 流入 BL 不整合、near-wall 未解像で
   壁モデル無し)。→ 定量化には (a) `inletProfile`+合成乱流、(b) **IDDES (候補3)** が前提。
   ⚠ keepDiss* は top-level キー ([memory: keepdiss-keys-toplevel-trap])。
3. **IDDES Phase 1.5 統合** ([turbulence-iddes-sst](../../plans/active/turbulence-iddes-sst.md) §4.8):
   LES 枝の基盤散逸を「ILES + keepDissType:2 σ0.02 jump2」とする設計に更新して実装
   (旧案の「WALE を LES 枝に」は今日の結果で見直し)。
4. **256³ 準 DNS** (RTX3060 で ~6h@CFL 引き上げ, メモリ ~4GB): forge 自身で DNS に肉薄し
   LES 設定を自前検証。検証インフラ (--dns 重ね描き・--node-mesh dedupe) がそのまま使える。
5. **回帰スイート整備** (小): naca_slau/rans_flat_plate_sst の baseline 陳腐化 (96c12c0 から未更新・
   変更前バイナリでも同値で FAIL = 既存問題) + ハーネスが古い `.build-native` バイナリを
   `build/` より優先する罠の解消。
6. 小粒: gauge の TP 枝 e∞ / SLAU / 種輸送 massflux / メトリック FP64 (Karp,
   [memory: fp32-metric-freestream-fix]) / node periodic dt 半減 (setDT が合併前 half-CV 体積) /
   定常局所 dt 不安定の深掘り (運用ルールで回避中) / node 多成分 L1。

## 罠・ルール (恒久)

- **SGS を使う run は vis_turb≠0 を必ず確認** (WALE 不活性を見逃した教訓)。
- 過渡ケースは `unsteady: 1`。unsteady:0 で数 step 発散したら CFL を疑う前にこれ
  ([memory: steady-localdt-acoustic-transient-instability])。
- node の KE 等の体積積分は周期 slave ミラー重複を除外 (`plot_dissipation_rate.py --node-mesh`)。
- 市松振幅は共分散射影 (`case/35/plot_checkerboard_decay.py`) でのみ測る。
- 式をコードに書く前に python で数値検証 (`tools/verify_advective_gauge.py` /
  `verify_sigma_model.py` が型)。実装前にオフラインで効果を測る (再構成ジャンプの事前測定が型)。
- config struct 変更後は full rebuild ([memory: stale-build-struct-layout-trap])。
- forge 実行は `tools/run_case.sh` 経由・新 run dir・README run 表同期。
- DNS 比較は実データ (`case/09/ref_dns/`) のみ (近似帯の誤りで判定を間違えた事故あり)。
- node dt は周期 half-CV の CFL 律速で cell の半分 (fixed-dt なら影響なし)。
- ユーザー方針: **node ベース主体**。新機能・検証は node 一次、cell は回帰対照。

## 関連 memory (自動ロードされる)

keep-es-dissipation-status / wale-inactive-fix / advective-gauge-uref /
steady-localdt-acoustic-transient-instability / user-prefers-node-base / verify-node-and-cell-both

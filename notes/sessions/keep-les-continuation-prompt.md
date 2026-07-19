# 引き継ぎプロンプト: KEEP-LES 安定化の続き (2026-07-19 セッションから)

> 新セッションの冒頭にこのファイルを読ませる。前セッションで「KEEP を LES/DES で安定に回す」基盤が完成した。
> ここから先の候補タスクと、踏んだ罠・運用ルールをまとめる。

## 完了済み (再調査・再実装しないこと)

1. **KEEP 用 entropy-stable 散逸レイヤ — 完成・accepted 済**
   [`plans/accepted/convection-keep-es-dissipation.md`](../../plans/accepted/convection-keep-es-dissipation.md) が正本 (経緯・全検証 run・教訓)。
   - config: `keepDissType` (0=off/1=scalar/2=**matrix 推奨**)・`keepDissCoeff` (σ)・`keepDissCprime` (既定1, lowMachPrecond から独立)
   - **σ の使い分け (較正済)**: 低マッハ市松対策=0.05 (既定) / **WALE 併用の解像 LES=0.02-0.03**
   - CPG/TP 単成分/多成分すべて検証済。多成分 matrix のプラトーバグは根治済 (真因: ΣρY≠ρ 共通モードノイズ → カーネル内 Y 正規化で除去)
   - 市松 null-mode の実証・診断ケースは `case/35.uniform_periodic_box` README の run 表 (L1)、TGV は `case/09.Taylor-Green` (L2/L3)
2. **pRef free-stream 修復の KEEP 展開 — 済**
   `case/33.wavy_hex_freestream` run_0004-0006: KEEP は pRef 無しで step11 発散 → `space.pRef=動作静圧` で 3.96e-12 機械精度。散逸レイヤ込みでも無害。
   plan: [`plans/active/convection-freestream-preserving-flux.md`](../../plans/active/convection-freestream-preserving-flux.md) (まだ active、残課題あり)
3. **背景サーベイ**: [`notes/investigations/convection-central-scheme-oscillation-control.md`](../investigations/convection-central-scheme-oscillation-control.md)
   (振動源4分類・検証ラダー L0-L5・優先度表。外部レビュー反映済み)

## 次の候補タスク (優先度順の提案)

1. **U∞≠0 の動く一様流保存** (激安・次のボトルネック判定):
   `case/33` に U≠0 一様流 run を追加し、移流項 (ρu·s, ρH u·s) の float32 桁落ちが KEEP でどれだけ効くか実測。
   効くなら移流項の基準差分 (ρ_ref u_ref 差し引き) を設計 — freestream-preserving plan の残課題そのもの。
2. **IDDES Phase 1.5 統合**: [`plans/active/turbulence-iddes-sst.md`](../../plans/active/turbulence-iddes-sst.md) §4.8。
   f_d ブレンドの LES 枝の基盤散逸として keepDissType=2 (σ~0.02) を使う設計に更新して実装。
3. **64³ TGV + digitize した DNS 参照で L3 定量化** (σ 最適の確定。32³ では層流期と終値を同時に満たせない=解像度律速と判明済)
4. **メトリック FP64 化 (Karp 2506.05150)**: pRef+移流対策で足りない場合の本丸。RTX は FP64 1/32 なので
   「セットアップだけ double・格納 float」等の設計検討から ([memory: fp32-metric-freestream-fix])
5. 小粒: ROE/AUSM への pRef 展開 / node periodic の dt 半減の解消 (setDT が合併前 half-CV 体積を使う、効率のみ) /
   node モードでの多成分 L1 (未実施)

## 罠・ルール (このテーマ固有)

- **市松振幅の測定は必ず共分散射影** (`case/35/plot_checkerboard_decay.py` 使用)。サブセット射影の
  sum(parity)≠0 DC バイアスで「継ぎ目散逸弱化」を誤認した事故あり (存在しなかった)
- **「強い散逸で問題が見えない」は無罪証明にならない** (scalar が多成分バグを隠した)
- **別カーネルで更新される保存量の比 (Y=ρY/ρ 等) をカーネルで使うときは正規化必須**
- 式をコードに書く前に python で数値検証する (`solver_density_cuda/tools/verify_entropy_scaling*.py` が型)。
  初出の $R^{-1}$ 形は誤りだった (正: $R|\Lambda|SR^{\mathsf T}\Delta w$)
- config struct 変更後は full rebuild ([memory: stale-build-struct-layout-trap])
- forge 実行は必ず `solver_density_cuda/tools/run_case.sh <run_dir>` 経由・新 run dir・README run 表同期
- node dt は cell の半分 (周期 half-CV の CFL 律速)。減衰を「ステップ数」で比較しない (単位時間で比較)

## KEEP-LES の推奨設定 (現時点)

```yaml
solver: "KEEP"
keepDissType: 2        # matrix ES 散逸
keepDissCoeff: 0.02    # 解像 LES (WALE 併用)。市松頑健性優先なら 0.05
keepDissCprime: 1      # マッハ混在流。単一低マッハ領域なら 0 + σ そのまま
space: {pRef: <動作静圧>}   # 非直交メッシュでは必須 (無いと発散)
turbulence: {LESorRANS: 1, LESmodel: 1}   # WALE
```

## 関連 memory (自動ロードされる)

keep-es-dissipation-status / fp32-metric-freestream-fix / backstep-unsteady-not-convergence-failure /
node-mode-periodic-and-backstep-status / verify-node-and-cell-both

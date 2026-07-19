# 引き継ぎプロンプト: KEEP-LES 安定化の続き (2026-07-19 セッションから)

> 新セッションの冒頭にこのファイルを読ませる。前セッションで「KEEP を LES/DES で安定に回す」基盤が完成した。
> ここから先の候補タスクと、踏んだ罠・運用ルールをまとめる。

## 完了済み (再調査・再実装しないこと)

1. **KEEP 用 entropy-stable 散逸レイヤ — 完成・accepted 済**
   [`plans/accepted/convection-keep-es-dissipation.md`](../../plans/accepted/convection-keep-es-dissipation.md) が正本 (経緯・全検証 run・教訓)。
   - config: `keepDissType` (0=off/1=scalar/2=**matrix 推奨**)・`keepDissCoeff` (σ)・`keepDissCprime` (既定1, lowMachPrecond から独立)
   - **σ の使い分け (較正済)**: 低マッハ市松対策=0.05 (既定) / **WALE 併用の解像 LES=0.02 (64³+実DNS で確定, 候補3参照)**
   - CPG/TP 単成分/多成分すべて検証済。多成分 matrix のプラトーバグは根治済 (真因: ΣρY≠ρ 共通モードノイズ → カーネル内 Y 正規化で除去)
   - 市松 null-mode の実証・診断ケースは `case/35.uniform_periodic_box` README の run 表 (L1)、TGV は `case/09.Taylor-Green` (L2/L3)
2. **pRef free-stream 修復の KEEP 展開 — 済**
   `case/33.wavy_hex_freestream` run_0004-0006: KEEP は pRef 無しで step11 発散 → `space.pRef=動作静圧` で 3.96e-12 機械精度。散逸レイヤ込みでも無害。
   plan: [`plans/active/convection-freestream-preserving-flux.md`](../../plans/active/convection-freestream-preserving-flux.md) (まだ active、残課題あり)
3. **背景サーベイ**: [`notes/investigations/convection-central-scheme-oscillation-control.md`](../investigations/convection-central-scheme-oscillation-control.md)
   (振動源4分類・検証ラダー L0-L5・優先度表。外部レビュー反映済み)

## 次の候補タスク (優先度順の提案)

1. ~~**U∞≠0 の動く一様流保存**~~ **済 (2026-07-19)**: 実測で「M0.1 でも KEEP/SLAU 全発散・
   エネルギー ρuH 項が支配」と確定し、移流基準差分ゲージ (`space.roRef`+`space.uRef`) を
   KEEP CPG×cell に実装・機械精度根治 (case/33 run_0007-0013)。正本:
   [freestream-preserving plan §8](../../plans/active/convection-freestream-preserving-flux.md)。
   **残**: TP 枝 e∞ / node 境界半割面 (config で fail-fast 中) / SLAU / 種輸送 massflux。
   非直交メッシュ LES では pRef に加え roRef/uRef=平均流を設定すること。
2. **IDDES Phase 1.5 統合**: [`plans/active/turbulence-iddes-sst.md`](../../plans/active/turbulence-iddes-sst.md) §4.8。
   f_d ブレンドの LES 枝の基盤散逸として keepDissType=2 (σ~0.02) を使う設計に更新して実装。
3. ~~**64³ TGV + digitize した DNS 参照で L3 定量化**~~ **済 (2026-07-19)**: 実 DNS
   (Dairay+2017 512³, `case/09/ref_dns/`) との 64³ 比較で **σ=0.02 を解像 LES 推奨として確定**
   (K/K0 誤差 ≤1.6%、σ=0 と DNS を対称に挟む)。σ=0.03 は撤回 (旧 DNS 近似帯が誤りだった)。
   σ は「L1 を満たす最小値」で選ぶ (最適は解像度と共に低下)。case/09 README 64³ 節参照。
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
keepDissJump: 1        # 再構成後ジャンプ (2026-07-19): 市松無傷のまま KE コスト 1/4 (終値-0.8%)
                       #   (plans/accepted/convection-keep-diss-recon-jump.md, node/cell 両検証済)
keepDissCprime: 1      # マッハ混在流。単一低マッハ領域なら 0 + σ そのまま
space: {pRef: <動作静圧>, roRef: <動作密度>, uRef: [<平均流速>, 0, 0]}
  # 非直交メッシュでは pRef 必須 (無いと発散)。U∞≠0 (平均流あり) なら roRef/uRef も必須
  # (M0.1 でも移流桁落ちで発散: plan §8)。CPG×cell のみ対応 (node/TP は fail-fast)
turbulence: {LESorRANS: 1, LESmodel: 1}   # WALE
```

## 関連 memory (自動ロードされる)

keep-es-dissipation-status / fp32-metric-freestream-fix / backstep-unsteady-not-convergence-failure /
node-mode-periodic-and-backstep-status / verify-node-and-cell-both

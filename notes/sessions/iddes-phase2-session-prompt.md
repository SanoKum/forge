# 引き継ぎプロンプト: DES 単一スキーム基盤完成 → IDDES Phase 2 (2026-07-19 セッション末)

> 新セッションの冒頭にこのファイルを読ませる。前身は
> [keep-les-continuation-prompt.md](keep-les-continuation-prompt.md) (KEEP-LES 基盤)。
> 本セッションで「advGauge node 対応 → Turkel 前処理散逸 → 単一スキーム KEEP+ES の
> RANS/DDES 検証 → DDES の LES 化実証」まで完了し、IDDES Phase 2 が着手可能になった。

## 完了済み (再調査・再実装しないこと; commits a648e50 → f2ea72c)

1. **advGauge の node 対応**: 境界半割面 (`convectiveFlux_boundary_d`) に upwind 差分形ゲージ。
   非直交+平均流の node が機械精度 (case/33 `run_0017`-`0020`)。
   [freestream plan §8.5](../../plans/active/convection-freestream-preserving-flux.md)。
2. **`keepDissPrecond: 1` (Turkel 前処理音響散逸)**: σ/c' の一様固有値スケーリングは「市松と
   2-4Δ 物理を分離できない単一トレードオフ曲線」に縮退することを実測 (case/18 `run_0109`-`0118`)
   → 前処理 2×2 (Δp 散逸 ∝c²/Ur 増強・ΔUn 散逸 ∝Ur 縮小) で**真の市松を c' 版の 142 分の 1 に
   減衰** (case/35 `run_0043`/`0044`)・ES 維持 (sym(K) 全域 PD)・物理コスト c' 以下。
   [plan](../../plans/active/convection-keep-diss-lowmach-precond.md)、式検証
   `tools/verify_precond_dissipation.py`。**backstep の残存縞 ~22 Pa は数値市松でなく
   リップ Δx=0.13 の解像限界物理と確定** (根治はメッシュ細分化のみ)。
3. **SGS パイロット** (case/18 `run_0122`-`0124`): ILES/WALE/σ が共通 3D 乱流 seed から
   cfl1.0 安定・vis_turb 正常 (σ は WALE の 1.7 倍散逸的)。定量 x_R はスコープ外と判断
   (入口 BL 不整合・壁モデル無し)。
4. **DES Phase 1.5 = 全域単一スキーム KEEP+ES 採用** ([iddes plan §4.8 設計更新](../../plans/active/turbulence-iddes-sst.md)):
   f_d flux ブレンドは廃 (フォールバック温存)。前提検証合格:
   - flat plate SST (case/26 `run_0012`/`0013`): KEEP+ES+implicit cfl20 安定・Cf 完全収束
     (ドリフト ≤0.04%/20k)・**Cf/Schl 0.88/0.91/0.94 = SLAU 比 −3.6〜−3.9% (既知のスキーム間
     ばらつき帯)**。u+-y+ も SLAU と重なる (`run_0013/uplus_yplus_slau_vs_keepes.png`)。
   - DDES 機能確認 (case/18 `run_0125`): 3D 857k+implicit 500step NaN 無し・f_d zoning 正常。
5. **DDES の LES 化実証** (case/18 `run_0126`, unsteady RK3 cfl1 30k): せん断層 f_d=1.000 で
   **解像 3D 乱流がゼロから発達・飽和** (Uz rms 0→~5 m/s ≈13%U∞)、モデル νt は RANS 比 79%
   (粗 Δ の Smagorinsky 相当水準)、付着 BL 遮蔽保持。**入口分布 A/B** (`run_0127`,
   `inletProfile:1`+壁法則 channel): 乱流強化 (Uz rms 8)・瞬時再循環短縮 — 統計未 (下記候補2)。

## 推奨設定 (確定版, 2026-07-19 末)

```yaml
solver: "KEEP"
keepDissType: 2          # ★ keepDiss* は top-level キー! space 配下は黙って無視される
keepDissCoeff: 0.05      # 壁付き低マッハ; 解像 LES (TGV 系) は 0.02
keepDissJump: 2
keepDissPrecond: 1       # Turkel 前処理 (cprime を上書き)。市松 142×・ES 維持
space: {pRef: <動作静圧>, roRef: <動作密度>, uRef: [<平均流速>,0,0]}  # cell/node 両対応 (KEEP CPG)
turbulence: {LESorRANS: 2, RANSmodel: 1, DESmode: 1}   # DDES (LES 単体は LESorRANS:0/1)
time: {unsteady: 1}      # DES/LES の発達は物理 dt 必須 (unsteady:0 では乱流が育たない)
```

## 次の候補タスク (優先度順)

1. **IDDES Phase 2 (本命)**: [iddes plan §3.4/§4.6](../../plans/active/turbulence-iddes-sst.md)。
   壁距離依存 $\Delta_{IDDES}$、$f_B$/$f_{dt}$/$f_e$、$l_{IDDES}$、`DESmode: 2`。
   実装前に python で関数形状を Shur 2008 と突き合わせ (verify スクリプト方式)。
   検証: DDES 回帰 (DESmode:1 ビット不変) → 周期チャネル WMLES (log 層/mismatch) →
   backstep (RANS 入口 BL + 剥離 LES + 再付着後 WMLES)。
   ★ 並行セッション作の [turbulence-wmles-wall-stress plan](../../plans/active/turbulence-wmles-wall-stress.md)
   (代数壁応力モデル) が active にある — IDDES と相補なので着手前に整合を確認。
2. **x_R の時間統計化**: run_0126/0127 の再付着比較を定量に (check_quasisteady への x_R 量追加
   or 時間平均統計)。瞬時値は激しく非定常で報告不可 (0127 は下流 x≈9-11 に二次剥離パッチも)。
3. **precond plan の残**: TGV 64³ L2 (KE cost) 実測、TP 枝展開。
   [plan §6](../../plans/active/convection-keep-diss-lowmach-precond.md)。
4. **SBLI へ伸ばす前提**: 衝撃センサ駆動 σ ブースト (JST k₂ の ES 版, Ducros×σ ランプ) を別 plan
   (現 ES レイヤは shock capturing なし。backstep 級の低マッハは不要)。
5. 小粒: 回帰ハーネスの `.build-native/relwithdebinfo/forge` 古 symlink 優先罠 /
   naca_slau・rans_flat_plate_sst baseline 陳腐化 / freestream plan 残課題 (TP/SLAU/種) /
   backstep リップ細分化メッシュ (縞の根治、Δx 0.13→0.06)。

## 罠・ルール (恒久; 今日ハマったもの含む)

- **keepDiss\* は top-level キー**。space 配下に書くと黙って無視され pure KEEP になる
  (A/B がビット同一なら設定キー位置を疑う) [memory: keepdiss-keys-toplevel-trap]。
- **restart 直後の res_0 は vis_turb=0** (restart 非忠実性)。νt は res_0 の k/ω から ρk/ω で復元。
- backstep DES メッシュ系は **visc=0.001** (νt/ν 正規化を 1.8e-5 でやると桁を誤る — 今日誤報告済)。
- flat plate の Cf 基準は **run_0009 (EWT 細メッシュ)**。run_0003 は陰解法回帰テンプレートで
  Cf 未検証 (理論+30%・x 無依存) — A/B の土俵にしない。
- **Cf/x_R の収束はドリフト実測で判定** (10k step では全然足りない: Ux 1%/5k。40-60k で
  ≤0.04%/20k まで落ちる)。check_quasisteady + スナップ間ドリフトを必ず貼る。
- SST×explicit RK3 は SLAU では発散した組だが **KEEP+ES では安定** (run_0126)。
- 入口分布は `ints:{inletProfile:1}` + `inlet_profile_<physID>.csv`
  (`tools/gen_inlet_walllaw.py` で生成)。k/ω floats は inlet 種別共通で有効。
- forge 実行は `tools/run_case.sh` 経由・新 run dir・README run 表同期。cwd はバックグラウンド
  コマンド後にリセットされることがある — run dir は絶対パスで扱う。

## 関連 memory (自動ロードされる)

keepdiss-keys-toplevel-trap / keep-es-dissipation-status / advective-gauge-uref /
user-prefers-node-base / verify-node-and-cell-both / backstep-unsteady-not-convergence-failure /
forge-sst-restart-nonfidelity

# time_integration: line-implicit v2 試作 (粘性結合 + K 凍結) — 上限確認の限定実験

- status: done
- 起票: 2026-09-02
- related_docs: [methods/time_integration/implementation.md](../../methods/time_integration/implementation.md)
- 先行: [plans/accepted/time_integration-line-implicit.md](../accepted/time_integration-line-implicit.md) (v1: 対流のみ、DDES A/B で不採用)

## 目的 (製品化ではない)

v1 の DDES A/B (case/39 ny160, run_diag_lineimp_*) で確定したのは「**現行 v1 は割に合わない**」
(サブ反復収束 +8% / step 単価 2.44 倍) であって line-implicit の到達点ではない。本 plan は
**上限確認のための一段だけ**を区切って実施し、「現実装が遅かった」のか「この離散化・ケースでは
line-implicit の上限自体が低い」のかを決着させる。

損益分岐: 現収束性能 (nSub 20→17 相当) のままなら line のサブ反復単価を ctrl の **~1.18 倍以下**に
落とす必要がある。K 凍結だけで届かなければ、粘性結合による大幅な nSub 削減が必須。

## 実験項目

1. **K 凍結 (`lineKFreeze`)**: K ブロック抽出 (単位ベクトル×5 の Jacobian 列展開 = v1 の主コスト仮説)
   を物理 step の最初のサブ反復のみ実行し、以後のサブ反復で再利用。純粋な計算コスト下限を測る。
   近似の質: defect-correction の LHS はもともと近似であり、サブ反復間の状態変化は O(ΔQ)。
2. **壁法線 line へのスカラー粘性結合**: line 面に限り K_prev/K_next へ −α_f·I を加算
   (α_f = ν_eff·δ/dcc, 既存 `viscous_diag`=2α と同じメトリック)。対角は既存の 2α を維持
   (行和優位を保つ安全側。真の対称形 [−α, 2α, −α] が line 内で完成する)。
3. **粘性 CFL 割引 (`lineViscousDtRelief` θ∈[0,1])**: on-line セルに限り `setDT` の粘性
   スペクトル半径項 2ν_eff/(ρ·dx_min) を (1−θ) 倍。θ=0 (現状) → 0.5 → 0.8 → 1.0 で
   安定限界と必要 nSub を測る。dx_min は壁近傍で壁法線=line 方向なので方向整合は近似的に成立。
4. **同一プロトコル比較**: case/39 ny160・発達場 (run_0023 res_16000)・300 step 窓・
   総 wall time で ctrl (point) と比較。測るもの: (a) K 凍結後のサブ反復単価比、
   (b) θ を上げたときの安定限界、(c) 同等残差に必要な nSub、(d) 総時間が ctrl を下回るか。

## 判定

- **勝ち筋**: (K凍結後単価) × (必要 nSub / 20) < 1 × ctrl 単価 → line-implicit は本実装
  (完全粘性 5×5 Jacobian・adaptive nSub) を検討する価値がある。
- **負け筋**: θ を上げても nSub が削れない/不安定 → 「対流 defect-correction が律速で、
  粘性 line 陰化の上限は低い」と結論し、v2 は記録して閉じる。

## 変更ログ

- 2026-09-02: 起票。
- 2026-09-02: 実装・検証完了 (**勝ち筋で決着 — 「現実装が遅かった」が正**)。case/39 ny160 DDES
  発達場・300 step 窓 (`run_diag_lineimp2_*`, IC=run_0023 res_16000, ctrl=point 345 s):
  - **コスト分解**: v1 モノリシック 843 s (2.44×) の主犯は **Thomas の毎 sweep 再分解**
    (rhs 非依存の LU + W + Kprev·W 625 積を 5 回/subiter 重複)。factor/solve 分離 (厳密) で
    547 s (1.59×)、`lineKFreeze` (サブ反復間凍結) で **456 s (1.32×) = コスト下限**。
    split/kfreeze は v1 と収束軌道が一致 (分離の厳密性・凍結の無害性を同時検証)。
  - **粘性結合・dt 割引はこのケースでは僅差** (roe +0.02 桁)。理由: 壁第一セル (Δn≈3.5e-5 m)
    の λ_visc/λ_acoustic = 2ν/(Δn·c) ≈ **0.02** — 圧縮性 pseudo-dt の壁法線律速は
    **音響 (c/Δn) であり v1 の対流・音響 line が既に陰化している**。粘性 Δn² 律速が主役に
    なるのは Δn < 2ν/c (本ケース y+≈0.02 相当) の超極薄セルのみ。θ=1.0 でも安定。
  - **本命は pseudo-CFL 引き上げ**: point は cfl_pseudo 2 で発散 (step 99)、**line は 2/4/8 全て
    300 step 安定**。cp4 で ctrl@20 品質に subiter 12-13 到達 (ctrl は 19)、cp8+nSub20 は
    roUx 17 倍深い・ωバースト最大 10.8→8.5 に低減。
  - **実測の勝ち**: `cp4 + nSub13` = **302 s (ctrl の 0.88 倍) で ctrl 同水準の step 終端残差**
    (rms_roUx 5.6e-7 vs 8.0e-7)。同時間 (nSub15) なら ctrl@20 より深い。
  - **判定**: M6 定常ノズル (streamwise 対流律速 → line 無効) と DDES dual-time (壁法線音響律速
    → line 有効) で**律速モードが違う**。DDES 生産推奨 = `lineImplicit:1 + lineKFreeze:1 +
    cfl_pseudo 4 + nSub 13-15` (粘性結合/割引は任意)。
  - **未検証の注意**: 300 step 窓のみ — nSub 削減のバースト余裕 (run_0013 の教訓: バースト最悪値は
    乱流発達と共に成長) は長時間 run で要検証。cp8×nSub 削減の組合せ、adaptive early-exit は今後。
- 2026-09-03: **方向別 dt (`lineDtDirectional: 1`) 追試** — ユーザー指摘「陰化した方向の λ は
  Δτ を縛らなくてよいはず」の実証。line 面の λ (音響込み) を CFL の max から完全除外し、
  壁セルの Δτ を off-line (streamwise) 基準 (×AR≈74, BDF 物理項が対角支配する Newton 的
  レジーム) にする。結果 (`run_diag_lineimp2_dir_cp4*`, 300 step):
  - **cp4+nSub20 で安定** (478 s)。収束は line_cp4 と同等〜微改善 (roUx 3.13→3.37 桁)、
    **ωバースト最大 9.5→6.15 に低減** (壁 Δτ 拡大が k-ω 隔離更新の突発を均す)。
  - **nSub8 は品質不足** (0.77/1.35 桁 < ctrl@20 の 1.53/2.17) — 方向別 dt でもサブ反復の
    収縮率は上がらない。**収縮律速はもはや壁擬似時間でなく off-line lag / SST segregated 結合**。
    ctrl@20 品質に必要な nSub は 12-13 のまま。
  - 運用示唆: 同品質の最速点は cp4+nSub13 (~0.87×) で不変だが、**directional はバースト余裕を
    +50% 積み増す**ので nSub 削減時の安全マージンとして併用推奨。
- 2026-09-03: **記録の補正 (Codex レビュー 3 点)**:
  1. **「壁セル Δτ が AR≈74 倍」は現実装では不成立** — `setDT` の directional 除外は
     `line_prev/next` に一致する**内部 line 面だけ**で、壁ノードの**境界半割面は CFL の max に
     残る**。壁 CV は内部 line 面と境界面が同じ V/S (実測 1.408e-5 m) なので、境界面が同じ音響
     制約を残す。**AR 倍化は壁の 1 つ内側以降の line CV のみ**。壁境界面 λ の除外可否
     (Dirichlet 行 decouple 済みなら安全かもしれない) は未検証の将来項目。
  2. ωバースト半減の帰属は「directional dt により半減 (実測)」に留める — 「壁 Δτ 拡大が均した」
     は未分離の推論 (1. の通り壁セル自身の Δτ は伸びていない)。
  3. 収縮律速の候補は **3 者**: off-line lag / segregated SST / **2次 KEEP RHS×1次 FVS LHS の
     defect-correction 不整合** (M6 の cfl×relax 飽和と同根の可能性)。
- 2026-09-03: **推奨構成の直接 A/B (directional cp4 × nSub13/15)** — 結果は下記追記。
- 2026-09-03: **推奨構成の直接 A/B 完了** (`run_diag_lineimp2_dir_cp4_nsub13/15`, 300 step):
  - `directional cp4 + nSub13`: **317 s = ctrl の 0.92 倍**、step 終端品質は同等〜良
    (roUx 4.4e-7 vs ctrl 8.0e-7, roe 2.45e-4 vs 2.13e-4)、ωバースト最大 6.15 (ctrl 10.8)。
  - `directional cp4 + nSub15`: 363 s (1.05 倍) で全量 ctrl より良い (roUx 2.0e-7, roe 1.33e-4)。
  - 非 directional nSub13 (302 s, 0.88 倍) より directional は +5% 遅いが、バースト余裕 −35%
    (9.5→6.15) を買う取引。**生産候補 = directional cp4+nSub13 (速度優先) / nSub15 (品質・余裕優先)**。
  - 残作業 (本採用の条件): 生産再開時に nSub15 で長時間 (数万 step) のバースト余裕検証を 1 本。

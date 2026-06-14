# 対流流束の free-stream 保存 (基準静圧差分)

## メタ

- **area**: `convection`
- **status**: `in_progress`
- **related_docs**:
  - `docs/convection/implementation.md` (TODO: 反映)
- **related_plans**:
- **created**: `2026-06-14`
- **owner**: CFD Dev

## 1. 目的

forge (float32) は非直交メッシュで一様場 (free-stream) を保持できず発散する。原因は
対流流束の圧力項 `Σ_f p_tilde·s_f` を float32 加算する際、非直交セルでは面ベクトル和が
厳密に相殺せず (Σs≈1e-7·面積)、大きな `p·s`(各項 ~1e3)の桁落ちが偽の運動量源
(~1e-4) になり増幅されること。基準静圧 `pRef` を差し引いて `(p_tilde − pRef)·s` で組み、
一様場で大きな `p·s` 自体を作らせず free-stream を機械精度で保存する。

検証は `case/33.wavy_hex_freestream` (歪み非直交hex 一様静止) で残差が machine zero に
落ちることを回帰とする。

## 2. スコープ

- **やる**:
  - `space.pRef` を config に追加 (既定 0.0 = 従来挙動・ビット不変)。
  - SLAU 対流流束 (内部 `SLAU_d` + 境界 `convectiveFlux_boundary_d`) の運動量圧力項を
    `(p_tilde − pRef)·s` 化。`__constant__ d_pRef` + wrapper の `cudaMemcpyToSymbol`。
- **やらない (別 plan)**:
  - 移流項 (質量 ρu·s, エネルギー ρH u·s) の基準差分 (U_∞≠0 の一様流の完全保存)。
    本件は IC が U=0 のケースでは不要 (基準=静止)。
  - ROE/AUSM/KEEP 系カーネルへの展開 (まず SLAU を検証)。
  - **inlet/outflow 流れありでの発散 (別機構)**。本 plan の範囲外。下記「残課題」。

## 3. 設計方針

- 物理的にはゲージ (定数圧) シフトで、`Σs=0` (解析的) なので解は不変。float32 桁落ちのみ除去。
- 一様基準は**全面共通の定数** pRef (= 動作/フリーストリーム静圧)。両側セルで同値なので
  内部面は等価逆符号のまま → **保存性を厳密に保つ** (局所セル値 U_c 差分は保存を崩すため不採用)。
- 重要: 「セル残差で precomputed r_c を一括減算」は float32 では不十分 (case/33 で 27% 減のみ)。
  大きな `p·s` を作って丸めた後では取り戻せない。**面ごとに引く** (本実装) のが要点。

## 4. 実装 (済)

- `input/solverConfig.hpp`: `flow_float pRef = 0.0;`
- `input/solverConfig.cpp`: `pRef = getOptionalValidatedValue<double>(space,"pRef",0.0,...)`
- `cuda_forge/convectiveFlux_d.cu`:
  - `__constant__ flow_float d_pRef;`
  - `SLAU_d` / `convectiveFlux_boundary_d` の運動量に `p_tilde_r = p_tilde - d_pRef`。
  - `convectiveFlux_d_wrapper` 先頭で `cudaMemcpyToSymbol(d_pRef, &pRef_h, ...)`。

## 5. 検証結果 (2026-06-14)

一様静止 (ρ=1.177, U=0, P=101325, 全slip) で、pRef=101325 を入れた効果:

| ケース | pRef なし | pRef=101325 |
|---|---|---|
| case/33 wavy hex 一様残差(運動量rms) | 5.9e-5 → step6 発散 | **3.9e-12, 300歩安定** |
| fan 実hex 全slip静止 | step7 発散 | **4e-14, 200歩安定** |
| StaticMixer tet+prism 全slip静止 | step3 発散 | **4.6e-12, 200歩安定** |

→ free-stream 保存崩れが(静止場)発散の原因と確定。基準圧差分で機械精度保存・発散解消。
(エネルギー残差は ~1e-6 まで線形微増 = 機械精度ドリフト。エネルギー基準差分は未実施・別 plan。)

## 6. 残課題 (別 plan 化)

- ROE/AUSM/KEEP への pRef 展開 (現状 SLAU のみ)。
- `docs/convection/implementation.md` への反映と回帰テスト登録。

## 7. 「流れありで発散」の決着 (2026-06-14) — pRef とは別問題だった

当初「inlet/outflow で流れを与えると発散」を free-stream/inlet-3D の問題と疑ったが、
**backstep の動く実メッシュに同設定を載せても発散** したことで、原因は **実行設定** と確定。
真因は3つ (いずれも pRef/コンバータ/メッシュとは無関係):

1. **explicit (timeIntegration 3) + 定常 (unsteady=0) は局所時間刻みで不安定**。
   → 定常は **implicit (11, blockDPLUR, cfl_pseudo≈1, nStepInner≈10)**、explicit なら unsteady=1。
   (静止場は流れが無いので explicit+定常でも安定 = 当初 free-stream に見えていた)。
2. **静止場に速度入口を当てると危険**。初期場の流速の向き・速さを入口とおおむね合わせる。
3. **outlet_statPress は逆流 (Un<0) 分岐で Ptb/Ttb を使う**。Ps のみ指定だと Pt=0→ρ=0→NaN。
   逆流用 Pt/Tt を設定する (StaticMixer の混合域過渡逆流での step137 発散はこれが原因)。

実証: **fan (Fluent hex) 完全収束** (case/32 run_0011_proper)、
**StaticMixer (Fluent tet+prism) 300歩健全完走** (case/31 run_0012_backflow)。
hex/tet/prism の Fluent メッシュが forge で計算可能と実証。レシピは
`.github/forge-calculation-workflow.md` の Fluent 節に記載。

## 変更ログ

- 2026-06-14: SLAU に pRef 差分を実装。case/33 (歪みhex 静止) で free-stream を machine zero 化。
- 2026-06-14: 「流れあり発散」を切り分け、原因は実行設定 (explicit+定常 / 静止IC / 逆流 Pt/Tt 欠落)
  と確定。pRef は静止場の free-stream にのみ寄与。Fluent メッシュ (fan/StaticMixer) の実行を実証。

## 変更ログ

- 2026-06-14: SLAU に pRef 差分を実装。case/33・fan・StaticMixer の静止場発散を解消 (machine zero)。
  流れありの発散は別機構として残課題に分離。

# 対流流束の free-stream 保存 (基準静圧差分)

## メタ

- **area**: `convection`
- **status**: `in_progress`
- **related_docs**:
  - `methods/convection/implementation.md` (TODO: 反映)
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

- ~~KEEP への pRef 展開~~ **済 (2026-07-19)**: `KEEP_d` の運動量圧力項 Gtilde を `(Ps−d_pRef)` 化
  (エネルギー Ptilde は −pRef∇·u で方程式が変わるため対象外 = SLAU と同判断)。検証
  `case/33...run_0004` (pRef 無 = step11 発散) / `run_0005` (pRef = **3.96e-12 機械精度・300step**) /
  `run_0006` (matrix ES 散逸込みでも 4.07e-12 = 散逸レイヤは一様場で厳密ゼロ)。
  **KEEP-LES を非直交実メッシュへ載せる前提条件が成立**。
- ROE/AUSM への pRef 展開 (現状 SLAU/KEEP)。
- ~~移流項 (ρu·s, ρH u·s) の基準差分~~ → §8 で実測・設計 (KEEP に実装)。SLAU の移流ゲージは未実装 (§8.1 の実測どおり SLAU も U∞≠0 で発散するが、upwind 切替の非線形性で因数分解が難しく別課題)。
- `methods/convection/implementation.md` への反映と回帰テスト登録。

## 8. 移流項の基準差分 (U∞≠0 の動く一様流保存, 2026-07-19)

### 8.1 実測 — 移流桁落ちは pRef 済みでも致命的 (次のボトルネック確定)

`case/33` の wavy hex を x 方向 periodic 化 (xmin/xmax は平坦・(y,z) 重心が厳密一致、
最近傍マッチングで対合成立) し、U∞=M0.1 (34.72 m/s, +x) の動く一様流を投入した:

| run | 設定 | step0 rms (ro / roUx / roe) | 結果 |
|---|---|---|---|
| `run_0007_keep_pref_u0_xperiodic` | KEEP+pRef, U=0 (periodic 化対照) | 1.2e-14 / 4.0e-12 / 3.6e-9 | 300step 安定 = run_0005 と一致 (periodic 無害) |
| `run_0008_keep_pref_uxM01_xperiodic` | KEEP+pRef, **U∞=M0.1** | 3.7e-8 / 1.3e-6 / **1.3e-2** | **step14 発散** |
| `run_0009_keep_dissmat_uxM01_xperiodic` | + matrix ES 散逸 σ=0.05 | 同上 | **step13 発散** (散逸レイヤでは救えない) |
| `run_0010_slau_pref_uxM01_xperiodic` | SLAU+pRef, U∞=M0.1 | 3.8e-8 / 1.4e-6 / 1.1e-2 | **step6 発散** (upwind でも救えない) |

- **帰属の定量確認**: メッシュの float32 セル閉じ誤差 $r_c=\sum_f \pm s_f$ (|r_c| max 1.9e-9) から
  予測した残差 rms (質量 $\rho u\, r_c$=1.7e-8 / 運動量 $\rho u^2 r_c$=5.9e-7 / エネルギー
  $\rho u H\, r_c$=5.1e-3) が実測と全成分で因子 ~2 で一致し、残差比も u (35.6≈34.7)・
  H (3.5e5≈3.0e5) と厳密に一致 → **step0 残差は純粋に移流流束×メトリック閉じ桁落ち**。
- エネルギー式が支配的 ($\rho u H \sim 1.2\times10^7$ ≫ p の 1e5)。
- Cartesian 対照は case/35 (run_0014/0016: 同 M0.1 through-flow・全周期で 400step 安定・機械ゼロ)
  が既存 → wavy メトリック起因で確定。

### 8.2 設計 — 参照一様流束の面ごと差分 (差分形因数分解)

参照状態 $(\rho_\infty, \mathbf u_\infty, p_\infty{=}\texttt{pRef})$ の厳密流束
$F_\infty(s) = (\rho_\infty U_{n\infty},\ \rho_\infty U_{n\infty}\mathbf u_\infty,\ 
\rho_\infty U_{n\infty} H_\infty)\,S$ ($U_{n\infty}=\mathbf u_\infty\!\cdot\!\mathbf n$) は
**s に線形な定数係数流束**なので、面ごとに引いても (i) 内部面で等価逆符号のまま=保存厳密、
(ii) セル和は $F_\infty(\sum_f s_f)=0$ 解析的にゼロ=解不変のゲージ。pRef の教訓どおり
「大きな積を作ってから引く」のでは桁落ちが戻らないため、**全項を (差分)×(平均) + (参照)×(差分)
に因数分解して組む** (一様場では全因子がビット単位ゼロ):

- 質量: $\tilde C - C_\infty = \overline{(\rho-\rho_\infty)}\,\overline U_n + \rho_\infty\,\overline{(\mathbf u-\mathbf u_\infty)}\!\cdot\!\mathbf n \equiv dC$
- 運動量: $dC\,\overline{\mathbf u} + C_\infty\,\overline{(\mathbf u-\mathbf u_\infty)}$ (+既存 Gtilde の (p−pRef))
- エネルギー:
  - $\tilde K$: $dC\,\overline{k} + C_\infty\,\tfrac12[(\mathbf u_0-\mathbf u_\infty)\!\cdot\!\mathbf u_1 + \mathbf u_\infty\!\cdot\!(\mathbf u_1-\mathbf u_\infty)]$
  - $\tilde I$: $dC\,\overline e + C_\infty\,\overline{(e-e_\infty)}$, $e_i-e_\infty=(p_i/\rho_i - p_\infty/\rho_\infty)/(\gamma-1)$
    ($p_\infty/\rho_\infty$ はカーネル内で同一 float 除算 → 一様場でビット零)
  - $\tilde P$: $\tfrac12[(\mathbf u_0-\mathbf u_\infty)p_1 + \mathbf u_\infty(p_1-p_\infty) + (0\leftrightarrow1)]\!\cdot\!\mathbf n$
    ((p−pRef) を再利用。$p_\infty U_{n\infty}$ の定数流束なので方程式は不変 —
    §4 で見送った「Ptilde の pRef 差分」($-p_{ref}\nabla\!\cdot\!u$ で方程式が変わる) とは別物)
- **スカラー輸送用 `massflux[ip]` は物理流束に戻して格納** ($dC + C_\infty$)。ゲージは 5 式のみ
  (係数が場に依存する $Y_f$ には $F_\infty$ が telescoping しないため。種の free-stream 桁落ちは残課題)。
- **config**: `space.roRef` (既定 0=off・ビット不変) + `space.uRef: [x,y,z]`。有効化は roRef>0。
  参照は **全面共通の定数** (pRef と同じ理由で局所参照は保存を崩すため不採用)。
- **スコープ**: KEEP CPG 枝 + cell モード (境界 ghost 面も主ループ経由で同一処方が乗る)。
  TP 枝 ($e_\infty$ に thermo 評価が要る)・node 境界半割面 (`convectiveFlux_boundary_d`)・SLAU は
  フォローアップ (ゲージは CV の**全**面に載らないと telescoping が破れるため、node は境界
  カーネル側の対応とセットで入れる)。node+roRef>0 は config で fail-fast。

### 8.3 実装・検証 (2026-07-19, 済)

- **式の事前数値検証**: `tools/verify_advective_gauge.py` — (1) 因数分解の解析等価 (float64
  rel err 2.4e-14)、(2) 一様場 float32 で全項ビット単位ゼロ、(3) wavy 内部面セル和で
  gauged 残差厳密ゼロ + orig−gauged ≡ F∞(r_c) の telescoping、の 4 検証 PASS。
- **実装**: `solverConfig.{hpp,cpp}` (`space.roRef`/`space.uRef`, 既定 0=off・roRef>0 は
  pRef>0 必須 + node で fail-fast) / `convectiveFlux_common_d.cuh` (`d_roRef`,`d_uRef*`) /
  `convectiveFlux_d.cu` (wrapper 転送) / `convectiveFlux_keep_d.inc.cuh` (advGauge 枝)。
  `massflux[ip]` は $C_\infty s$ を戻して物理流束のまま (種輸送整合)。散逸レイヤ・implicit
  Jacobian は定数流束の減算で不変のため無変更。
- **検証** (`case/33`, 詳細は README run 表):
  - `run_0011_keep_advgauge_uxM01_xperiodic`: U∞=M0.1 で step0 rms 1.7e-14/8.1e-12/5.6e-9 =
    **U=0 (run_0005) と同水準の機械精度**・300step 安定。ゲージ無し (run_0008: roe 1.3e-2,
    step14 発散) から**エネルギー 6.5 桁改善・発散根治**。ulp 級の roundtrip ずれ (roUx/ro,
    EOS) は「別の一様状態」として定数係数×s に載り telescoping されるため機械精度を壊さない。
  - `run_0013_keep_advgauge_dissmat_uxM01`: matrix ES σ=0.05 併用でも 3.3e-7・300step 安定。
  - 回帰: `run_0012` (roRef=0, atomicAdd ノイズ床 4e-12 内で run_0007 一致) /
    `case/35 run_0037` (node KEEP 市松 = run_0023 と同挙動) / tests/regression の
    implicit_bump_slau・tp_sod_n2 PASS (naca_slau・rans_flat_plate_sst の FAIL は
    変更前バイナリでも同値の**既存 baseline 陳腐化**で本件と無関係)。

## 7. 「流れありで発散」の決着 (2026-06-14) — pRef とは別問題だった

当初「inlet/outflow で流れを与えると発散」を free-stream/inlet-3D の問題と疑ったが、
**backstep の動く実メッシュに同設定を載せても発散** したことで、原因は **実行設定** と確定。
真因は3つ (いずれも pRef/コンバータ/メッシュとは無関係):

1. **静止場に速度入口を当てると危険**。初期場の流速の向き・速さを入口とおおむね合わせる
   (最大の主因。fan/StaticMixer の早期発散はこれ)。
2. **outlet_statPress は逆流 (Un<0) 分岐で Ptb/Ttb を使う**。Ps のみ指定だと Pt=0→ρ=0→NaN。
   逆流用 Pt/Tt を設定する (StaticMixer の混合域過渡逆流での step137 発散はこれが原因)。
3. **初手から難条件 (2次移流・乱流・no-slip) にしない**。段階起動 (引き継ぎ) で上げる。
   時間積分は陰解法 (blockDPLUR) が定常収束に有利だが、陽解法+定常も可 (発散主因ではない)。

実証: **fan (Fluent hex) 完全収束** (case/32 run_0011_proper)、
**StaticMixer (Fluent tet+prism) 300歩健全完走** (case/31 run_0012_backflow)。
hex/tet/prism の Fluent メッシュが forge で計算可能と実証。レシピは
`procedures/calculation-workflow.md` の Fluent 節に記載。

## 変更ログ

- 2026-07-19: **移流項の基準差分 (advGauge) を KEEP CPG×cell に実装** (§8)。U∞=M0.1 の動く
  一様流で全スキーム発散する移流桁落ちを実測・帰属確定 (run_0008-0010) し、参照一様流束の
  差分形因数分解ゲージ (`space.roRef`/`uRef`) で機械精度保存に根治 (run_0011, 6.5桁改善)。
  式は `tools/verify_advective_gauge.py` で事前数値検証。残: TP 枝 / node 境界 / SLAU / 種輸送。
- 2026-06-14: SLAU に pRef 差分を実装。case/33 (歪みhex 静止) で free-stream を machine zero 化。
- 2026-06-14: 「流れあり発散」を切り分け、原因は実行設定 (静止IC と入口流れの不整合 / 出口逆流
  Pt/Tt 欠落 / 初手から難条件) と確定。pRef は静止場の free-stream に寄与。陰解法は収束に有利だが
  陽解法+定常も可 (発散主因ではない)。Fluent メッシュ (fan/StaticMixer) の実行を実証。

## 変更ログ

- 2026-06-14: SLAU に pRef 差分を実装。case/33・fan・StaticMixer の静止場発散を解消 (machine zero)。
  流れありの発散は別機構として残課題に分離。

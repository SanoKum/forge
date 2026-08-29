# axis-Mach チェーンの粘性 δ\* 補正 (A12 = 親計画の A7): inviscid 壁 → 物理壁

## メタ

- **area**: `tooling / optimization`
- **status**: `done`  <!-- 2026-08-16 起票・同日 v1/v2 完了 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (「axis-Mach チェーン」節)
- **related_plans**:
  - 親: [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md)
    (§6 A7「着手時に後続 plan を起票して分離」の実行)
  - 資産: [`tooling-nozzle-phase3-windtunnel.md`](tooling-nozzle-phase3-windtunnel.md)
    (NS v1 = run_0028–0030 の起動レシピ・δ\*/相関比 1.04–1.09 の実測)
- **created**: `2026-08-16`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

axis-Mach チェーン (inviscid、‖ΔM‖∞ 0.240% $M_d$・設計 0.39 s) の出力壁は
**非粘性の流線**であり、実機では境界層の排除厚 $\delta^*(x)$ 分だけ有効断面が細る。
物理壁 $r_{\rm phys}(x) = r_{\rm inv}(x) + \delta^*(x)/\cos\theta_w$ を与えて
RANS (SST) でも軸 M 分布が目標に乗る状態を作る。

- **v1 (相関)**: 乱流平板相関 + Eckert 参照温度 (`feedback/deltastar.py`、実装済) を
  cplus 壁へ適用。CFD 不要。B8 系の実測では相関は真値の 4–9% 過小 (run_0030)。
- **v2 (CFD 抽出)**: v1 壁の RANS 場から $\delta^*_{\rm CFD}(x)$ を抽出し、壁を
  $r_{\rm inv} + \delta^*_{\rm CFD}$ で置き直す固定点反復。

## 2. スコープ

- **やる**: `runner_axismach` の RANS 経路 (runner_wt の `_config_sst_node`・coarse 中継・
  ω 床・段階起動レシピを流用)、δ\* 抽出ユーティリティ (質量流束欠損積分 +
  エッジ検出、探索窓はフリーストリームまで [run_0030 の教訓])、v1→v2 1 反復の実測
- **やらない**: 固定点反復の完全自動化 (収束判定・パス管理) — v2 1 反復の結果を見て判断。
  凝縮・実在気体。②以降の機種

## 3. 前提 (B8 系 NS v1 で確立済みのレシピ)

- **coarse SST 中継が必須**: Euler 場 + 一様 k/ω から y+1.5 へ直接入ると出口コーナー
  roOmega NaN (run_0028 で 3 連敗)。y+~50 (wall_first_frac 5e-3) で SST を回して
  「構造の正しい k/ω/BL 場」を作り、そこから y+1.5 へ cross-mesh restart。
- **y+1.5 メッシュ**: wall_first_frac 4.5e-5, ni 561, nj 97 (AR ≤ ~880)。
  y+1 (3e-5) は AR 1533 でゲート FAIL。
- **ω の粘性底層フロア**: 粗→細の最近傍 interp は近壁 ω が 5 桁不足のまま貼るので
  ω = 6ν/(β1 y²) を床に (runner_wt 実装済み)。
- **背圧は整合圧**: p_ambient = Pt/151.8 ≈ 6588 Pa (1000 Pa のままだと出口∩壁で吸い出し NaN)。
- **本段は cfl 1** (M4 低圧 y+1.5 固有。cfl 2/4 は出口コーナーで不可)。
- **δ\* 抽出の探索窓はフリーストリームまで届かせる** (0.35 r* 窓は下流で切れて
  非物理な「δ\* 減少」を出した — run_0030 で撤回済み)。x<8 はコア流が半径方向に
  未一様でエッジ検出が破綻するため測らない (相関にフォールバック)。

## 4. 設計方針

1. **pass0 (v1)**: cplus inviscid 壁 → `deltastar_offset` (相関) → 物理壁 →
   coarse 中継 → y+1.5 本段 → 軸 M・$\delta^*_{\rm CFD}$ 測定
2. **pass1 (v2)**: $\delta^*_{\rm CFD}$ を平滑化 (spline) し、x<8 は相関へブレンド →
   物理壁を置き直し → 同じ 2 段 (IC は pass0 の y+1.5 場を cross-mesh interp) →
   軸 M 改善と $\delta^*$ の固定点性を確認
3. ゲート: RANS 軸 M の ‖ΔM‖∞ (テスト部 x∈[x_E−2, x_E]) と、
   $\max|\delta^*_1 - \delta^*_0|/\delta^*_0$

## 4.1 実測 (2026-08-16 自走)

**v1 (相関 δ\* 壁)**: coarse 中継 `run_0069_ns_v1_coarse` (43 s, 品質 PASS・
`check_convergence --drop 2` ALL PASS [roK 3.8 桁/roOmega 3.4 桁]・NaN 0) →
y+1.5 本計算 `run_0070_ns_v1` (162 s, 48k step cfl1)。
**品質 PASS (AR≤1000)・ALL PASS (全 7 列 2.6–3.3 桁)・ALL STEADY・NaN 0**
(res_48000 の ro/P/T も物理的)。

- **軸 M (無帰還)**: ‖ΔM‖∞ = **0.533% $M_d$**、rms 0.0125、overshoot +0.30%。
  B8 系 NS v1 (run_0030) の 1.2% $M_d$ を半減以下に更新。
- **δ\*_CFD/相関比: 中央値 1.095、範囲 [1.043, 1.153]** (58/58 ステーション測定成功、
  x∈[8, 22.4])。B8 系の実測 1.04–1.09 と同傾向 (相関は過小、下流ほど乖離増)。
- 出口 $\varepsilon_M$ の全断面評価は BL 込みで 7% になる — **粘性では出口一様性の
  評価を BL 除外 (コア内) に変える必要がある** (残件)。

## 5. 完了条件

- [x] RANS 経路 + δ\* 抽出の実装、単体レベルの検算
- [x] pass0 (v1) 完走・健全性 VERDICT・δ\*/相関比の測定
- [x] pass1 (v2) 完走・固定点の確認
- [x] methods 更新・README run 表・commit

## 4.2 v2 と結論 (2026-08-16)

**v2** (`run_0071_ns_v2`: 壁 = inviscid + δ\*_CFD、IC = v1 場の cross-mesh、161 s):
品質 PASS・conv ALL PASS・ALL STEADY・NaN 0。

1. **δ\* は 1 反復で固定点**: v2 場から抽出した δ\* は入力に対し比 0.995–1.010
   (中央値 1.003)。固定点反復の自動化は不要。
2. **軸 M は v1 と不変** (‖ΔM‖∞ 0.533 vs 0.533% $M_d$、プロファイル全点 1e-3 以内で一致
   [max 0.0213 @ x=5.93 が両者で同一])。**残差は δ\* で表現できる誤差ではない** —
   Euler でも見えた設計側残差 (x≈6 の谷) + 近スロートの粘性効果で、消すには
   law 側の帰還 (A5 機構の RANS 版) が要る。
3. **実務上の結論: 相関 δ\* (v1) で十分**。δ\*_CFD との差 (~10%) は軸 M に観測可能な
   影響を持たない。生産チェーン = inviscid 設計 (0.39 s) + 相関オフセット (ミリ秒) +
   RANS 検証 1 回 (205 s)。RANS-in-the-loop は不要。
4. v1 の無帰還 0.533% $M_d$ は B8 系 NS v1 (run_0030, 1.2%) の半分以下。

**初回 v2 の失敗記録**: δ\* 測定域端 (x=22.4) の端値クリップ外挿が折れ目を作り、
壁 B-spline が下流端で非単調化 → 端勾配の線形外挿 + 全域再平滑化で解消
(`analyze_ns_deltastar.py`)。

## 残件

- RANS 軸 M の 0.5% ゲート化には law 側帰還 (A5 の RANS 版) が要る (設計側残差が支配)
- 粘性の出口一様性評価 (BL 除外・コア内評価) は未実装 — 全断面 $\varepsilon_M$ は
  BL 込みで 7% になり指標として使えない
- coarse 中継の省略可否 (SST→SST warm start なら直接 y+1.5 で立つか) は未検証

## 6. 変更ログ

- `2026-08-16` — 起票 (ユーザ指示「2→1 で自走」の 1)。
- `2026-08-18` — case/44 (M4.19 加熱空気, TP semi-perfect, r_t 0.21 m) で NS チェーンを TP×SST×node に初適用 (cfl 1 で安定)。
  相関 δ\* は本 case では**上流で 30–50 % 過大** (case/41 の 10 % 過小と逆向き; 強い加速域で平板相関が過大) → v1 出口 M −0.5 %。
  x<8 も BL 縁で ρu が極大になり δ\* が測れる (x≥3) ため、`prepare_ns(dstar_blend=)` を追加し CSV を全域採用 (x<3 は測定端の比で
  相関をスケール)。δ\* は 2 反復で固定点 (1.001)、‖ΔM‖∞ 0.52→0.31 % Md、出口コア M 4.183 ±0.002。末端の抽出勾配 (出口 BC 影響) は
  トレンド勾配で外挿する (`case/44 build_dstar_v3_va.py --x-hi-trust`)。残差 (近スロート −0.2〜−0.3 %) は δ\* 非表現で law 側帰還が要る
  (残件のまま)。詳細: [case/44 README](../../case/44.vitiated_air_wt/README.md) NS 節。
- `2026-08-30` — case/42 M5 (R3/L_U9/L_c14, Tt 1550 K) に適用。v1 相関は中流 1.29 倍過大→出口 0.94 倍
  過小で出口コア −0.71 %; v3 (build_dstar_v3_m5.py = va 版流用) 1 反復で −0.29 %・固定点 1.001。
  **抽出信頼下限は case 依存を再確認**: M5 は x_lo=12 必須 (x=3 で 10 倍の異常抽出。case/44 の x≥3 は
  「コア一様性が十分」の条件付きだった)。run_0101 で「生抽出 CSV 全域採用→壁波打ち (軸 M −5 % ディップ)」
  という本 plan 記録済みの v2 失敗を再演した — **反復は必ず build_dstar_v3 系 (平滑化+外挿) を使うこと**。
  残件に追加: 本フローの forge_design 本体への配管 (analyze/build スクリプトが case コピーのまま) と
  x_lo の自動判定 (抽出/相関比の sanity ゲート)。
- `2026-08-30` — case/42 M6 (R2/L_U12/L_c50, Tt 1550 K) にも適用 (M5 と同日・同フロー)。v1 出口コア −0.57 % →
  v3 −0.43 %・固定点 0.998。抽出信頼下限は **x_lo=14** (M5 の 12、case/44 の 3 に続き case 依存を再々確認)。
  収束レベルが M5 (−0.29 %) より浅いのは M6 の BL が厚く (δ99 大) δ* 一次補正の限界が早いためとみられる。
  **メッシュの罠 (要申し送り)**: 第一セル 2.4 µm (wf 4.5e-5 × r_t 53.75 mm, ni1700) で median-dual 変換が
  閉性破綻 (`dual faces not closed` normalized 0.099, 壁 bcond 面積 3 % 欠損)。~3.5 µm では健全 →
  node y+1 双対幾何の float32 桁落ち修正 (既往) でも守れない下限がある。回避 = wall_first を上げ y+~1.4 で運用。

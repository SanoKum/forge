# 平衡凝縮の EOS 拘束形 (`condEquilibrium=2`): 湿り度を従属変数計算で代数的に決める

## メタ

- **area**: `condensation`
- **status**: `done`
- **related_docs**:
  - [`methods/condensation.md`](../../methods/condensation.md) (「平衡凝縮 — EOS 拘束形」節)
- **related_plans**:
  - [condensation-equilibrium.md](../accepted/condensation-equilibrium.md) (緩和形 `condEquilibrium=1`、本計画の前身)
  - [condensation-nonequilibrium.md](../accepted/condensation-nonequilibrium.md) (二相 EOS・rog 輸送の枠組み)
- **created**: `2026-08-19`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

緩和形の平衡凝縮は onset 直後 ~2 $r_t$ で θ 律速の遅れ (S≤1.18、過冷却 2 K; case/44 `run_0097`) が残る。
ユーザ指摘 (2026-08-19)「過冷却度がゼロとなるように β を収束計算させる」に対応し、液相分率 $g$ を**輸送量ではなく
状態量**として従属変数計算の中で $(T,g)$ 同時反転で決める EOS 拘束形を追加する。完了時: 飽和セルで各ステップ厳密に
$S=1$、下流の M/T/g は緩和形と一致 (固定点同一)、onset 帯の遅れが消える。

## 2. スコープ

- **やる**: `condEquilibrium: 2` (dependentVariables で $(T,g)$ 反転 → `rog` 射影 → 原始量同期; source kernel は診断のみで
  `res_rog`=0)。TP carrier (H₂O in 擬似種) / TP pure / CPG pure の 3 経路。case/44 va3old (旧条件×va3 形状) で緩和形と比較。
- **やらない**: 氷の飽和線、平衡音速で NS ブロック $\kappa$ を置換 (loose coupling のまま)、多凝縮種 (`nCondSpecies>=2` はエラー)、
  二温度、亜音速出口 ghost の二相整合。

## 3. 関連 docs と前提

- 二相 EOS $e=e_{gas}(Y,T)+g(R_wT-L)$、$p=\rho T(R_{mix}-gR_w)$、SLAU 面温度の $R_{gas}=R_{mix}-gR_w$ 補正 (methods §5)。
- SLAU 流束・出力は原始量 `g_0` (= `rog_0`/ρ) を読む → `rog` を射影したら `condensationPrimitive_d_wrapper` で同期が要る。
- 反復順序 (main.cpp `assembleResidual`): update → speciesPrimitive/condPrimitive → **dependentVariables** → bconds → 勾配 →
  流束 → condTransport/condSource → 時間積分 → condPointImplicit → condPrimitive。射影は dependentVariables 内で行うので
  流束・BC・出力は同一ステップ内で $g_{eq}$ を見る。

## 4. 設計方針

- **反転**: `cond_equilibrium_Tg_carrier(sp,nSp,Y,e_in,rho,Yw,Rw,cprops,T_guess,g_guess,&T,&g)`
  1. $T_0$=`thermo_T_from_e`(g=0)。$\rho Y_wR_wT_0\le p_{sat}(T_0)$ (or $Y_w\le0$) → $g=0,T=T_0$。
  2. 飽和: $G(g)=\rho(Y_w-g)R_wT(g)-p_{sat}(T(g))$、$T(g)$=`cond_T_from_e_carrier` (前回 T を warm start)。
     括弧 $[0,Y_w]$ ($G(0)>0$, $G(Y_w)<0$)、初期 $g$=前ステップ値をクランプ、括弧付き Newton (解析 $dG/dg$) + 二分法退避、
     収束 $|G|<10^{-8}p_{sat}$ or 括弧幅 $<10^{-10}Y_w$、最大 40 反復。
  - pure (TP/CPG): $p_v=\rho(1-g)R_vT$、$Y_w\to1$、$T(g)$=`cond_T_from_e_onetemp`/`cond_T_from_e_cpg`。
- **射影**: `rog[0][ic] = ρ g_eq` (dependentVariables kernel 内, nCells_all)。wrapper 末尾で `condEquilibrium==2` なら
  `condensationPrimitive_d_wrapper` を呼び `g_0` を同期 (realizability クランプも通る)。
- **輸送凍結**: source kernel の `eq==2` 分岐で診断 (`condS`,`condTsat`) を書いた後 `res_rog[ic]=0`, `sj_g[ic]=0`, `return`
  → point-implicit 更新 δ=0、`rms_rog` 恒等 0。$Q_0$–$Q_2$ は緩和形と同じ (輸送のみ)。
- **config**: `condEquilibrium` 0/1/2 のバリデーション、`nCondSpecies>=2 && condEquilibrium==2` はエラー。
- ビット互換: `condEquilibrium!=2` では dependentVariables の分岐に入らない (既存 0/1 経路不変)。

## 5. 実装ステップ

1. `input/solverConfig.{hpp,cpp}`: `condEquilibrium` の範囲 0–2、`nCondSpecies>=2` との排他。
2. `cuda_forge/condensationEOS_d.cuh`: `cond_equilibrium_Tg_carrier` / `_pure_tp` / `_pure_cpg`。
3. `cuda_forge/dependentVariables_d.cu`: kernel 引数 `condEquilibrium` 追加、TP/CPG 分岐で eq==2 のとき $(T,g)$ 反転・
   `rog` 射影 (以降の $e_{mix}$/$P$/$H_t$ 計算は既存式に $g_{eq}$ を流す)。wrapper で primitive 同期。
4. `cuda_forge/condensationSource_d.cu`: `eq==2` 分岐 (診断のみ・残差 0)。
5. `methods/condensation.md` (済) / 本計画 / `plans/README.md`。

## 6. 検証

- **ビルド**: `solverConfig.hpp` を触るので **full rebuild** ([[stale-build-struct-layout-trap]])。
- **回帰**: case/44 `run_0091` (dry) と `run_0097` (緩和形 eq=1) を新 binary で再実行 → dry はビット同一、eq=1 は変化なし。
- **本命**: `run_0097` と同条件で `condEquilibrium: 2` (`run_0098_va3old_Lc8_eq2`)。判定:
  - onset 以降の軸で $|T_{sat}-T|$ max ≪ 1 K (緩和形 2.0 K)、S∈[0.99,1.01] (緩和形 1.18)。
  - onset+2 $r_t$ 以降の軸 M/T/g が緩和形と一致 (|ΔM|≲1e-3, |ΔT|≲0.5 K)。
  - `check_convergence` / `check_quasisteady` / NaN / メッシュ PASS、軸 $h_0$ 保存 ≤5e-4。
  - 段階起動 (cfl 0.5→1→2) が緩和形と同じく完走 (不安定なら cfl を下げて記録)。
- 0-D 単体: host で `cond_equilibrium_Tg_carrier` を旧条件 (ρ,e) で呼び $S=1$、$g$ が緩和形の `cond_equilibrium_delta` 反復
  の固定点と一致することを確認 (任意)。

## 7. 影響範囲

- `condensation` ブロック (`condEquilibrium 2` のときのみ)。`dependentVariables_d` kernel の引数 1 個追加 (0/1 で無害)。
- docs: `methods/condensation.md` (追記済)、`plans/README.md`。

## 8. 完了条件

- [x] `methods/condensation.md` 更新
- [x] 実装・full rebuild・回帰 (dry / 緩和形 eq: 新旧 binary 差 ≤2e-6 相対 = node の run 間ノイズ床)
- [x] `run_0098` で判定基準を満たす、case/44 README に run 追記
- [x] `status: done`、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-08-19` — 起票 (ユーザ選択 (c): EOS 側 β Newton)。
- `2026-08-19` — 実装完了 (`condEquilibrium: 2`)。`cond_eq_solve_g` (括弧付き Newton+二分法) を `condensationEOS_d.cuh` に、
  carrier TP / pure TP / pure CPG の 3 経路。`cond_equilibrium_delta` (緩和形) も同ヘッダへ移動 (source kernel と単体テストで共用)。
  単体 `tests/unit/test_cond_equilibrium_eos.cpp` (nvcc host): 未飽和→g=0、飽和→S=1 & e 保存 (2e-10 J/kg)、
  緩和形 Δ 反復の固定点と g/T 一致 (g 0.00428, T 257.60 K)、初期値非依存、pure CPG N2 も S=1 → ALL PASSED。
  case/44 `run_0098_va3old_Lc8_eq2` (va3 L_c8 形状 × 旧条件, 段階起動 cfl0.5→1→2 そのまま完走, PASS/NaN 0/STEADY):
  **onset (x=6.21 r_t) 以降 全軸で S=1.0000・|T_sat−T|=0.000 K** (緩和形 `run_0097` は S≤1.18・2.0 K が ~2 r_t 続く)。
  下流 (x>8.5) は緩和形と |ΔM|≤1.0e-3・|ΔT|≤0.14 K・|Δg|≤1e-4 で一致 (固定点同一)、出口 g/Y 12.8 %・T 259.0 K・M 4.144 同値。
  onset 帯 (5<x<8.5) だけ ΔM −0.015・ΔT +1.8 K (θ 遅れの分)。`rms_rog` 恒等 0、軸 h0 保存 5.7e-4。回帰: 新 binary の
  dry (`run_0099`) / 緩和形 (`run_0100`) は旧 binary と相対 ≤2e-6 (node の run 間ノイズ床)。図 `figs/va3_eq2_vs_eq1_axis.png`。
  残件: 氷飽和線オプション、平衡音速の NS ブロック結合 (今回は不要だった)。

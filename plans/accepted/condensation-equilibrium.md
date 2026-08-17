# 平衡凝縮モデル (局所瞬時平衡) と飽和温度診断 `condTsat`

## メタ

- **area**: `condensation`
- **status**: `done`
- **related_docs**:
  - [`methods/condensation.md`](../../methods/condensation.md) (「平衡凝縮」節、モデル切替表)
- **related_plans**:
  - [condensation-nonequilibrium.md](../accepted/condensation-nonequilibrium.md) (4 モーメント非平衡、EOS 結合)
  - [condensation-evaporation.md](../accepted/condensation-evaporation.md) (蒸発、θ 律速の規約)
  - [tooling-nozzle-tp-split-h2o-condensation.md](../accepted/tooling-nozzle-tp-split-h2o-condensation.md) (split_h2o 経路)
- **created**: `2026-08-18`
- **owner**: `sano` (実装: Claude 自走)

## 1. 目的

風洞ノズル (case/44 加熱空気 M4.19–5.0) で「凝縮し得る熱力学的上限」を評価するため、核生成・成長を経ない
**局所瞬時平衡** の凝縮モデルを追加する。併せて過冷却度 $T_{sat}(p_v)-T$ を評価できるよう飽和温度を h5 に出力する。

## 2. スコープ

- **やる**: `condEquilibrium` フラグ (source kernel の平衡分岐、モーメントソース 0)、`condEqRelax`、診断 `condTsat_<s>`、
  0-D 検算 (平衡到達で S=1)、case/44 で非平衡 vs 平衡の比較。
- **やらない**: 氷の飽和線 (現状は過冷却液)、平衡液滴径、二温度、NS。

## 3. 関連 docs と前提

理論・式は methods の「平衡凝縮」節。二相 EOS (e = e_gas + g(R_w T − L))・蒸気分圧 $p_v=\rho(Y_w-g)R_wT$・
θ 律速 (Δg≤5e-3, ΔT≤1 K/step)・src_jac の規約は非平衡実装をそのまま使う。

## 4. 設計方針

- kernel 冒頭で $p_v$, $p_{sat}(T)$, $T_{sat}(p_v)$ を計算して診断に書き、`eq` なら平衡分岐へ (蒸発分岐・核生成より前)。
- $\Delta$ は $[-g, Y_w-g]$ で二分法 (40 反復) — $F$ は単調 (Δ↑で $p_v$↓, $T$↑→$p_{sat}$↑) なので一意。
- $S_g=\alpha\rho\Delta/\Delta t$、θ 律速は非平衡と共通、`sj_g = αθ/Δt`。$Q_0$–$Q_2$ ソース 0 (輸送のみ)。
- $T_{sat}$: $\ln p_{sat}(T)=\ln p_v$ を Clausius–Clapeyron 勾配 $L/(R_wT^2)$ の Newton で反転 (15 反復、$p_v<10^{-6}$ Pa は 0)。

## 5. 実装ステップ

1. `input/solverConfig.{hpp,cpp}`: `condEquilibrium` (int, 0), `condEqRelax` (double, 1.0)。
2. `cuda_forge/condensationProperties_d.cuh`: `cond_Tsat(cprops, pv)`。
3. `cuda_forge/condensationSource_d.cu`: `cond_equilibrium_delta(...)`、kernel に `eq, eqRelax, diagTsat` 引数、平衡分岐。
4. `variables.cpp`: 診断 `condTsat_` を登録 (h5 出力)。
5. `methods/condensation.md` 更新 (済)、case/44 で検証。

## 6. 検証

- ビルド (full rebuild — solverConfig.hpp 変更 [[stale-build-struct-layout-trap]])。
- 回帰: case/44 `run_0024` (非平衡, g=0) と同条件を新 binary で再実行し軸 M 一致。
- 平衡: case/44 M4.19 で `condEquilibrium 1` → S=1 の飽和線上に乗る (condS→1、Tsat−T→0)、g>0 が x≈6.4 r_t (S=1 到達) から立ち上がる。
  CEA 平衡 (氷基準、出口 8 % 析出, T 261.8 K) と定性一致 (液基準なので少なめ)。
- 判定: NaN 0、品質 PASS、非平衡 (M4.75 で onset) より上流で凝縮が始まること。

## 7. 影響範囲

- `condensation` ブロックのみ (既定 off でビット不変)。診断変数 1 つ追加 (凝縮 on のときのみ確保)。

## 8. 完了条件

- [x] 実装・ビルド・回帰一致 (case/44 run_0030: P 差 2e-6)
- [x] case/44 M_d 6 点の平衡 run と表 (run_0033/36/39/42/45/48)
- [x] methods / plans/README 同期、accepted へ移動

## 変更ログ

- `2026-08-18` — 起票 (ユーザ指示: 加熱空気風洞の平衡凝縮評価)。
- `2026-08-18` — 実装完了。`condEquilibrium`/`condEqRelax`/`condEqDTmax` (10 K)/`condEqDgMax` (0.05) と診断 `condTsat_<s>`。
  非平衡の θ 律速 (ΔT≤1 K/step) を平衡分岐に流用すると急膨張 (M_d 5) で S 67 が残るため平衡専用の上限を追加。
  case/44 va2 (Pt 1.137 MPa/Tt 1058 K) 6 点: onset 直後 ~1 r_t を除き S=1.00 (Tsat−T=0)、出口凝縮率 9 % (M4.19) → 45 % (M5.0)、
  M4.19 の平衡出口凝縮 8.4 % は CEA 平衡 (氷基準 8 %) と同程度。詳細は case/44 README「新条件」節。
- `2026-08-18` — **エネルギー非保存バグ修正** (methods「対流流束の面温度と二相圧力の整合」): SLAU TP 面温度が全水分気相の $R_{mix}$
  で再構成され凝縮帯で $h_0$ +0.6 %・T +4–5 K。出口ノードの「跳ね」はこの非整合の顕在化 (出口値が正しかった)。修正後 case/44 va2 の
  凝縮 12 run を `run_0067`–`0078` (va2fx) で再計算。`outlet_statPress` 超音速流出の ghost 全量外挿も同時に (結果不変)。

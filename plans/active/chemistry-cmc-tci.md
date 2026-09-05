# 乱流‐化学相互作用 (TCI): 1st-order CMC の導入

## メタ

- **area**: `thermophysics / chemistry / turbulence`
- **status**: `in_progress` (2026-09-05 起票・設計確定、Phase A 着手)
- **related_docs**:
  - `methods/chemistry_cmc.md` (現在仕様: 理論・実装設計)
  - `methods/chemistry.md` (有限速度化学の基盤、§5 TCI)
  - `notes/investigations/cabra-liftoff-model-fidelity-survey.md` (方式選定の文献根拠)
- **related_plans**:
  - `plans/active/chemistry-finite-rate-h2.md` (親。§5.1 P0-1/P0-4 の決定で本 plan を起票)
- **created**: `2026-09-05`
- **owner**: `CFD Dev`

## 1. 目的

laminar chemistry (セル平均で Arrhenius を評価) では自着火安定化火炎の着火遅れを表現できず、Cabra (case/48) で
リップ付着・BK (case/47) で早着火になる (機構交換・PaSR では直らない — 親 plan §9 2026-09-05)。
混合分率に条件付けた組成分布で反応を評価する 1st-order CMC を導入し、**Cabra の浮き上がり高さ H/d≈10 とその
コフロー温度応答、BK の着火位置**を再現できる状態にする。加熱器設計 (火炎位置・壁熱負荷まで) の前提。

## 2. スコープ

- **やる**: Bilger 混合分率の診断、分散輸送 + $\tilde\chi$ (SST)、$\eta$ 格子上の条件付きスカラー ($n_s$ 種 + $h$) の輸送
  (物理空間 1 次風上 + 乱流拡散、$\eta$ 拡散陰的 (AMC)、化学は各 $\eta$ 点で点陰解)、β-PDF 積分によるソース結合
  (`tci: 2`)、Cabra $H(T_c)$ 応答曲線と BK 着火位置での検証。定常 (擬似時間) と dual-time の両方で動く。
- **やらない**: 2nd-order CMC (条件付き分散)、密度重み付き CMC 形、多流 (3 流以上) 混合分率、LES-CMC、
  transported PDF (粒子/stochastic fields)、条件付き速度の高次モデル、$\tilde Y$ を PDF 積分で上書きする full coupling。

## 3. 関連 docs と前提

- 理論・設計: `methods/chemistry_cmc.md` (本 plan と同時に起草)。
- 既存資産: `scalarTransport_d` (1 変数の移流+拡散残差と時間積分)、`chemistry_d` (`chem_source`, 種ブロック LU,
  `jacobianInterval`)、`ransTransport_d` (SST の $k,\omega$ → $\varepsilon/k=\beta^*\omega$)、`thermo_d` (NASA-9, $h\leftrightarrow T$)。
- 前提: 全化学種の拡散係数を同一 ($Sc_t$ 共通) とする現行仕様のもとで元素質量分率は保存スカラー → $\tilde\xi$ は診断量。
- 制約: 条件付きスカラーは node 当たり $(n_s+1)N_\eta$ 個 (9 種 × 41 = 410 double)。Cabra 60k ノードで 25M 値 (200 MB) は可。
  化学評価は node 数 × $N_\eta$ 倍になる ($\approx40\times$) ので `jacobianInterval` と「$\tilde P(\eta_k)$ が閾値未満の $\eta$ 点は
  化学をスキップ」で抑える。

## 4. 設計方針 (差分のみ。本文は `methods/chemistry_cmc.md`)

1. **source coupling**: 流れ側は $\tilde Y,\tilde h$ を従来どおり輸送し、`chemistry_source_d` の $\dot\omega(\tilde Y,\tilde T)$ を
   $\sum_k w_k\tilde P(\eta_k)\dot\omega(Q(\eta_k))$ に置換する (`tci: 2`)。既存の陰解法配管を変えない。
2. **$\tilde\xi$ は診断** (Bilger、元素組成は機構 YAML の `composition`)、**分散だけ輸送** (`scalarTransport_d` 汎用拡散、
   生成 $2\mu_t/Sc_t|\nabla\tilde\xi|^2$、散逸 $\bar\rho\tilde\chi$ は `src_jac` で点陰解)。
3. **$\eta$ 格子**: $N_\eta$=41 (config)、$\xi_{MR}$ と $\xi_{st}$ 付近を tanh で密に。$\eta$ 端は Dirichlet (燃料/酸化剤流の状態)。
4. **$\eta$ 拡散 + 化学の陰解**: node 毎に「三重対角 ($\chi$ 項) + 各 $\eta$ 点の種ブロック LU」を分離反復 (先に $\eta$ 拡散を陰的に、
   次に化学を点陰解)。分離誤差は擬似時間収束で消える。
5. **PDF**: β 分布を $\eta$ 格子上で解析正規化して離散化 (台形則)。$\widetilde{\xi''^2}<\epsilon$ はデルタ扱い。
6. **AMC** で $\langle\chi|\eta\rangle$。$C_\chi=2$ (config)。
7. **精度**: $Q$・化学は double、輸送残差は flow_float。
8. **出力**: `xi`, `xiVar`, `chi`, `cmc_dY`。$Q_T(\eta)$ は指定 node 群のみ (診断)。

## 5. 実装ステップ

| Phase | 内容 | 主要ファイル | 状態 |
| --- | --- | --- | --- |
| **A** 混合分率インフラ | 元素組成読込 (`composition`)、Bilger $\tilde\xi$ 診断カーネル、分散 `roXiVar` の登録・輸送・ソース、$\tilde\chi$、出力。Cabra 混合場で検証 | `chemistry_mech_io.hpp`, `chemistry_d.cuh`, `cuda_forge/cmc_d.{cuh,cu}` (新規), `variables.cpp`, `main.cpp`, `input/solverConfig.*` | in_progress |
| **B** 条件付きスカラー | `Q` ストレージ・$\eta$ 格子・初期化 (混合線)・物理空間輸送 (scalarTransport ループ)・$\eta$ 拡散陰解 (AMC)・各 $\eta$ 点の化学点陰解・境界条件 | `cmc_d.{cuh,cu}`, `scalarTransport_d`, `chemistry_d` | todo |
| **C** ソース結合 | β-PDF 積分で $\bar{\dot\omega},\bar{\dot Q}$, 対角 Jacobian を置換 (`tci: 2`)、`cmc_dY` 診断、`jacobianInterval`/PDF 閾値スキップ | `chemistry_d.cu`, `speciesTransport_d.cu` | todo |
| **D** 検証 | 凍結混合整合 → Cabra $H(T_c)$ 応答曲線 (5 点) → BK 着火位置 → 回帰 (`tci 0/1` ビット不変) | `case/48`, `case/47`, `tests/` | todo |

### 5.1 残作業表 (優先順)

| 優先 | 項目 | 状態・次アクション |
| --- | --- | --- |
| 1 | Phase A: 元素組成読込 + Bilger $\tilde\xi$ 診断 | in_progress (2026-09-05) |
| 2 | Phase A: 分散輸送 + $\tilde\chi$ + 出力、Cabra 混合場 (`run_0067`) で $\tilde\xi$–$Y_{H_2}$ 整合と分散の妥当性 | todo |
| 3 | Phase B: $Q$ ストレージ・$\eta$ 格子・混合線初期化・凍結整合 (化学 OFF で $Q$ が線形) | todo |
| 4 | Phase B: $\eta$ 拡散陰解 (AMC) + 化学点陰解、境界条件 | todo |
| 5 | Phase C: ソース結合 (`tci: 2`)、`cmc_dY` | todo |
| 6 | Phase D: Cabra $T_c$ 1015/1030/1045/1060/1075 K の $H$ 応答曲線 (Y_OH>2e-4) vs 実験・文献、中心軸/半径 T ±30 K | todo |
| 7 | Phase D: BK 着火位置、`tci 0/1` 回帰、性能 (PDF 閾値スキップ) | todo |
| 8 | ドキュメント: `procedures/solver-settings.md` に `cmc` キー、`methods/index.md` | todo |

## 6. 検証

- **単体 / ビルド**: Bilger $\xi$ が純燃料/純酸化剤で 1/0、混合線上で線形 (ホストテスト)。β-PDF の正規化・モーメント
  ($\int P=1$, $\int\eta P=\tilde\xi$, $\int(\eta-\tilde\xi)^2P=\widetilde{\xi''^2}$) を $\eta$ 格子上で 1e-3 以内。
- **検証ケース**:
  1. 凍結混合 (case/48 化学 OFF): $Q_\alpha(\eta)$ = 混合線、`cmc_dY` < 1 %。
  2. Cabra H₂/N₂ (case/48): $H(T_c)$ 応答曲線 (5 点)。合格: 1045 K で H/d が 10 ± 3、傾きが実験 (10 K で約 2 倍) と同符号・同桁、
     $T_c\pm30$ K のどこかで H/d=10 を通る。中心軸・半径 T の平均差 ±30 K。判定は `check_quasisteady --quantity ignx` STEADY の場で。
  3. Burrows–Kurkov (case/47): 着火位置 18–25 cm (現状 4–6 cm)。
  4. 回帰: `tci: 0` / `tci: 1` の既存 run が ビット不変。
- **判定基準**: 上記。加えて質量・元素・全エンタルピー収支 (`check_quasisteady` の exit 量) が反応 ON でも閉じること。

## 7. 影響範囲

- 新規: `cuda_forge/cmc_d.{cuh,cu}`, `tools/test_cmc.cpp` (ホストテスト), `methods/chemistry_cmc.md`。
- 変更: `chemistry_mech_io.hpp` (composition), `chemistry_d.{cuh,cu}` (`tci: 2` 分岐), `variables.cpp` (`xi`, `roXiVar`, `Q`),
  `main.cpp` (呼び出し), `input/solverConfig.*` (`cmc` キー), `procedures/solver-settings.md`。
- 既定 (`tci: 0/1`) では全経路ビット不変。

## 8. 完了条件

- [x] `methods/chemistry_cmc.md` 起草
- [ ] Phase A–C 実装
- [ ] Phase D 検証 (§6) — Cabra 応答曲線と BK 着火位置
- [ ] `status: done` + §9 変更ログ、`plans/active/` → `plans/accepted/`、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-05` — 初稿。ユーザ決定「分布で考える TCI を入れて検証する」(親 plan §5.1 P0-1/P0-4) を受け、方式を 1st-order CMC
  (source coupling, Bilger 診断 + 分散輸送, AMC, β-PDF) に確定。根拠: 文献で Cabra $H(T_c)$ を再現した実績が最も強いのが
  transported PDF / CMC で、決定論的な CMC は定常 (擬似時間) ソルバと GPU の既存スカラー輸送・種ブロック陰解に乗る
  (stochastic fields は擬似時間定常と相性が悪い)。Phase A に着手。

# DES/DDES/IDDES 用 低散逸対流 flux 設計 技術調査（圧縮性 LES/DES）

<!-- ファイル名規約: turbulence-<short-slug>.md / 本文書は調査専用・コード変更なし -->

## メタ

- **area**: `turbulence`
- **status**: `draft`（**調査専用 / コード変更なし**）
- **related_docs**:
  - [`docs/turbulence/theory.md`](../../docs/turbulence/theory.md)
  - [`docs/turbulence/implementation.md`](../../docs/turbulence/implementation.md)
  - [`docs/convection/theory.md`](../../docs/convection/theory.md)
  - [`docs/convection/implementation.md`](../../docs/convection/implementation.md)
- **related_plans**:
  - [`turbulence-iddes-sst.md`](../../design/active/turbulence-iddes-sst.md) — 本調査の直接の動機（§4.8 低散逸 flux 設計を本調査で再評価・改訂）
  - [`turbulence-des-wmles-survey.md`](turbulence-des-wmles-survey.md) — DES/WMLES 手法選定サーベイ（SBLI は IDDES-SST+WENO が主流）
  - [`convection-pep-scheme-survey.md`](convection-pep-scheme-survey.md) — KEEP/PEP 系中心スキーム調査（中心枝の保存性側）
- **created**: `2026-06-22`
- **owner**: `CFD Dev`
- **調査実施**: deep-research ワークフロー（98 エージェント, 16 ソース取得, 25 クレーム検証: 24 確認 / 1 棄却）

> **位置づけ**: 本文書は実装計画ではなく **技術調査レポート**。forge に SST-DDES/IDDES を追加する際、LES 領域の対流 flux をどう設計するかを判断する材料（一次文献 + SU2 公式実装 + forge 現状との突き合わせ）をまとめる。実装に進む場合は本調査を基に [`turbulence-iddes-sst.md`](../../design/active/turbulence-iddes-sst.md) §4.8 の設計（既に本調査を反映して改訂済み）に従うこと。

---

## 0. エグゼクティブサマリ（結論先出し）

1. **標準設計は「LES 域=低散逸中心 / RANS・粗格子=風上」の連続ブレンド** `weight=(1-σ)w₁+σw₂`（Travin-Shur-Strelets-Spalart 2000 → OpenFOAM `DEShybrid` / SU2 `ROE_LOW_DISSIPATION`）。forge の `KEEP_SLAU` 中心↔SLAU ブレンドはこの系譜に正しく乗る。
2. **ブレンド係数の駆動量は「乱流解像の物理」（壁距離・速度勾配/渦度・渦粘性 = f_d/NTS）であり、衝撃センサとは直交した別物**。両者を混ぜてはいけない。
3. **forge 案 `duc = max(1-f_d, 1-ψ_limiter)` のうち `1-f_d` 項は SU2 `FD`（σ_FD=max(0.05,1-f_d)）と完全一致で妥当**。一方 **MUSCL limiter を衝撃センサに流用する設計は一次 DES/LES 文献に一例も無い**。
4. SU2 は衝撃検出に**専用 Ducros センサ**を別途用意（`FD_DUCROS`/`NTS_DUCROS`, σ=σ_Ducros+σ_NTS−σ_Ducros·σ_NTS）。ただし **`FD` 単体（f_d のみ）も正規の選択肢**で、Ducros は「実衝撃と剥離せん断層が共存するケース（=SBLI/ノズル/ピントル）」のベストプラクティスであって必須ではない（「両方必須」の過強主張は検証 0-3 で棄却）。
5. **推奨: 短期=`1-f_d` 主導で limiter は positivity フロアに格下げ、中期=改良 Ducros センサを OR ブレンドで追加（SBLI/ノズル/ピントル向け）、長期=LES 枝の LD2 型勾配再構成**。

---

## 1. 一次文献・実装事例の表

| source | solver / case | turbulence model | convective scheme | shock/discontinuity treatment | relevance to forge |
| --- | --- | --- | --- | --- | --- |
| **Travin, Shur, Strelets & Spalart (2000)** [doi:10.1007/0-306-48383-1_16](https://doi.org/10.1007/0-306-48383-1_16) | 原典 DES hybrid | DES (SA系) | LES 域=4次中心 / RANS・外部非回転=風上偏重(5/3次) の連続ブレンド | ブレンドは**乱流解像物理で駆動**、衝撃捕捉は別機構 | `KEEP_SLAU` 中心↔風上ブレンドの直接の原型 |
| **OpenFOAM `DEShybrid`** [API source](https://www.openfoam.com/documentation/guides/latest/api/DEShybrid_8H_source.html) | OpenFOAM 非構造FV | DDES | `weight=(1-σ)w₁+σw₂`, w₁=linear(低散逸), w₂=upwind-biased | σ は wall distance・velocity gradient・eddy viscosity で構築。**Ducros も MUSCL limiter も参照しない** | 写経できる連続中心↔風上ブレンドの実装テンプレ |
| **Molina et al. (SU2), AIAA 2017-4284** [arc.aiaa](https://arc.aiaa.org/doi/10.2514/6.2017-4284) | SU2 / BSCW 遷音速バフェット (AePW-2 Case 3, M=0.85, Re=4.49M) | DDES | Roe + adaptive σ | **σ=σ_Ducros+σ_NTS−σ_Ducros·σ_NTS**: 壁＋衝撃で σ=1（フル散逸）、せん断層/wake で σ→0.05 | **圧縮性遷音速 SBLI/バフェットの実 DES 例**。f_d+NTS+Ducros 共存の根拠 |
| **Löwe, Probst, Knopp, Kessler (DLR), AIAA J. 2016, J054956** [doi:10.2514/1.J054956](https://arc.aiaa.org/doi/abs/10.2514/1.J054956) | DLR TAU, **非構造FV 圧縮性** | hybrid RANS/LES | **LD2**（low-dissipation low-dispersion 2次）, 勾配情報で face 値 | 中心低散逸 + 選択的散逸 | **forge と同じ解像クラス（非構造FV・2次）で成熟**。長期の LES 枝再構成候補 |
| **Sciacovelli et al., Comp.&Fluids 230 (2021) 105134** [arXiv 2103.16426](https://arxiv.org/pdf/2103.16426) | implicit LES / 極超音速乱流BL | wall-resolved | 9次中心 + Jameson 型人工散逸 | **shock detector を「選択性向上・誤発火回避」へ改訂** | 「滑らか域=低散逸中心、衝撃のみ散逸」を高速圧縮性乱流で実証。**センサ選択性に投資する**のが定石 |
| **Hendrickson, Kartha & Candler, AIAA 2018-3710** [arc.aiaa](https://arc.aiaa.org/doi/10.2514/6.2018-3710) | 圧縮性 LES (衝撃あり) | — | — | **改良 Ducros センサ**。原 Ducros は「低渦度域で発散に過敏 → 解を著しく改変」する失敗モードを明記 | forge が Ducros を入れるなら**この改良版**を使うべき根拠 |
| **Sensor-restrained artificial shear viscosity, PRF 2025** [arXiv 2505.00158](https://arxiv.org/html/2505.00158v1) | 圧縮性乱流 | scale-resolving | LAD (Ducros×Heaviside(−θ)) | 過剰人工散逸が coherent structure・乱流統計を破壊 | **過散逸=解像乱流殺しのリスク**の一次根拠。低散逸要件の動機 |

凡例の検証票は §5 に記載（全て 3-0 確認、SU2 FD 定義のみ 2-1）。

---

## 2. 主要知見（一次文献ベース）

### 2.1 連続ブレンドの正準形（Travin / DEShybrid）

- ブレンドは `weight=(1-σ)w₁+σw₂`、`0≤σ≤σ_max`。w₁=低散逸中心（例: linear/4次中心）、w₂=風上偏重（頑健）。
- **σ の駆動量は乱流解像の物理（壁距離・速度勾配/渦度・渦粘性）。これは RANS↔LES の選択子であって、衝撃センサでも MUSCL limiter でもない**（DEShybrid の σ 式は Ducros・limiter・圧力/密度勾配のいずれも参照しない）。→ **forge が limiter をセンサ流用する設計の評価上、最も load-bearing な点**。

### 2.2 衝撃/不連続は専用センサで別処理

- 圧縮性 LES の標準的衝撃検出は **Ducros センサ** `σ_Ducros=(div u)²/((div u)²+ω²)`（LAD では負ダイレーションの Heaviside `H(-θ)` と併用）。
- **原 Ducros の失敗モード**（Hendrickson-Kartha-Candler 2018）: 低渦度域で発散に過敏になり、不要な散逸で解を著しく改変。改良版は低渦度感度を抑えつつ真の不連続では散逸を保つ。
- Sciacovelli 2021 は 9 次中心 + Jameson 型人工散逸で、**shock detector の選択性を上げて滑らか乱流域での誤発火を避ける**。**= 専用センサに投資するのが定石**。
- **MUSCL/slope limiter を衝撃センサに使う一次 DES/LES 文献は皆無**。limiter は「あらゆる急勾配（せん断層・乱流極値・接触面）」で発火し、衝撃と乱流を区別しない → DES が守りたい解像乱流をまさに過散逸させるリスク。

### 2.3 過散逸は解像乱流を殺す（低散逸要件の動機）

- Hendrickson 2018: 「散逸は安定化に必要な場所にだけ加える。過剰だと小スケール構造が解から消える」。
- PRF 2025: 「無制限の人工拡散は coherent 構造を劣化させ乱流統計を損なう。特に乱流 BL に沿って人工せん断粘性を過剰適用すると統計が崩れる」。Garnier 1999 / Johnsen 2010 も衝撃捕捉の小スケール圧縮性乱流ダンピングを裏付け。

---

## 3. SU2 専用節（compressible solver, `INC_*` 除外）

### 3.1 公式ドキュメント

- [Convective Schemes (docs_v7)](https://su2code.github.io/docs_v7/Convective-Schemes/) — `CONV_NUM_METHOD_FLOW`:
  - **Central**: `JST`, `JST_KE`, `JST_MAT`, `LAX-FRIEDRICH`
  - **Upwind**: `ROE`, `L2ROE`, `LMROE`, `TURKEL_PREC`, `AUSM`, `AUSMPLUSUP`, `AUSMPLUSUP2`, `SLAU`, `SLAU2`, `HLLC`, `MSW`
  - **`ROE_LOW_DISSIPATION` のサポートは `ROE` と `SLAU`/`SLAU2` のみ**（L2ROE/TURKEL_PREC/AUSMPLUSUP[2]/HLLC は非対応）。
- [config_template.cfg](https://github.com/su2code/SU2/blob/master/config_template.cfg): `ROE_LOW_DISSIPATION= FD ... for Hybrid RANS/LES simulations` → **DES 既定推奨 = `FD`（f_d ベース）**。

### 3.2 ソースファイルと config 名

- [option_structure.hpp](https://github.com/su2code/SU2/blob/master/Common/include/option_structure.hpp): `enum ENUM_ROELOWDISS { NO_ROELOWDISS=0, FD=1, NTS=2, NTS_DUCROS=3, FD_DUCROS=4 }`
  - `FD`="DDES F_d function", `NTS`="Travin and Shur", `NTS_DUCROS`="Travin+Ducros", `FD_DUCROS`="F_d+Ducros"
- [CNumerics.cpp `GetRoe_Dissipation`](https://github.com/su2code/SU2/blob/master/SU2_CFD/src/numerics/CNumerics.cpp):
  - `FD`: `Dissipation_ij = max(0.05, 1 - 0.5*(Diss_i+Diss_j))`、ノード値は f_d（[CNSVariable.cpp](https://github.com/su2code/SU2/blob/master/SU2_CFD/src/variables/CNSVariable.cpp) `SetRoe_Dissipation_FD`: `Roe_Dissipation=1-tanh((8·r_d)³)`）
  - `NTS`: `Roe_Dissipation = σ_max·tanh(A^ch1)`
  - `NTS_DUCROS`/`FD_DUCROS`: `Dissipation_ij = max(Min_Dissipation, φ₁+φ₂-φ₁·φ₂)`（確率OR 合成, σ_Ducros=(div u)²/((div u)²+ω²)）
- [roe.cpp](https://github.com/su2code/SU2/blob/master/SU2_CFD/src/numerics/flow/convection/roe.cpp): `F = ½(F_i+F_j)·n − σ·(½ P|Λ|P⁻¹ (U_i−U_j))` — **散逸項にだけ σ を乗じる**形（中心 flux は不変）。
- [ausm_slau.cpp](https://github.com/su2code/SU2/blob/master/SU2_CFD/src/numerics/flow/convection/ausm_slau.cpp): `CUpwSLAU[2]_Flow` のみ `val_low_dissipation` を受け `GetRoe_Dissipation` を呼ぶ。AUSM 系は非対応。

### 3.3 examples / 実例

- **BSCW（AePW-2 Case 3）遷音速バフェット DDES**（AIAA 2017-4284）: Roe + `NTS_DUCROS`。壁と衝撃で σ=1、せん断層/wake で σ≈0.05。平均 Cp と衝撃後方の coherent 構造で実験良一致。

### 3.4 SU2 の実際の推奨と制約

- **推奨パス: `CONV_NUM_METHOD_FLOW= ROE`（または `SLAU2`）＋ `ROE_LOW_DISSIPATION= FD`**（DES 既定）。実衝撃が共存するなら `NTS_DUCROS`。
- **制約**: 低散逸機構は Roe/SLAU(2) 限定。σ は**散逸テンソルへの乗数**（フルフラックスブレンドではない）。forge の `duc` がどちらの形か（散逸乗数 vs フルフラックス線形ブレンド）で SU2 σ 式の移植容易性が変わる（→ §6 未解決点）。

---

## 4. forge 評価（案 A–E）

- **案A（全域 SLAU/SLAU2 のまま）**: DES の意味を成さない。LES 域で風上散逸が解像乱流を殺す（§2.3）。**不可**。
- **案B（SLAU 散逸に σ を掛けて LES で弱める）**: SU2 `ROE_LOW_DISSIPATION` と同型で最小改修。ただし KEEP 中心枝を使わないため KE 非保存で渦保存に劣る。**代替**。
- **案C（KEEP 中心基本 + RANS/衝撃/positivity で SLAU ブレンド）**: 正準形（Travin/DEShybrid）そのもの。forge `KEEP_SLAU` で実装済の器に f_d を注入。**本命**。
- **案D（WENO/TENO 高次）**: 非構造 FV・GPU・float32 と非親和、実装侵襲大。**長期の選択肢外**（forge 方針と不整合）。
- **案E（SU2 推奨に寄せる）**: `FD`（f_d）主導 + 実衝撃で Ducros = **案C の中期形そのもの**。望ましい収束先。

**結論: 案C を本命に、中期で案E（Ducros 追加）へ拡張**。

### 4.1 forge 案 `duc = max(1-f_d, 1-ψ_limiter)` の判定

- **`1-f_d` 項**: SU2 `FD`（σ_FD=max(0.05,1-f_d)）と一致 → **妥当・採用**。フロア 0.05 も踏襲推奨。
- **`1-ψ_limiter` 項（limiter をセンサ流用）**: 一次 DES/LES 文献に支持なし。さらに forge 内部観測（[ducros-limiter-inert-with-slau] / [convection-ducros-limiter]）で **limiter on/off が SLAU 場を <0.1% しか変えない** = limiter はそもそもセンサとして不感。よって `1-ψ_limiter` は実質 `duc≈1-f_d` に縮退する公算大。**limiter は衝撃センサにせず positivity/robustness フロアに格下げ**。
- **SBLI/ノズル/ピントル（実衝撃＋剥離せん断層共存）**: SU2 BSCW がまさに `NTS_DUCROS` を要した状況 → **中期で改良 Ducros センサ（Hendrickson 2018）を明示追加**し、SU2 と同じ OR ブレンド `σ = σ_Ducros + (1-f_d) − σ_Ducros·(1-f_d)`。
- **float32 / 非構造FV / GPU リスク**: σ ブレンドは局所量から face ごとに計算する有界・平滑な乗数（フロア 0.05）で GPU 並列・float32 安定。**真のリスクはブレンドでなく既知の [forge-freestream-nonorthogonal]（非直交 free-stream 桁落ち）**。中心化で散逸が減るほど顕在化しうるので、LES 域の中心比率を上げる前に free-stream preservation を確認。

---

## 5. 検証票（deep-research 25 claims, 24 確認 / 1 棄却）

- 3-0 確認: DEShybrid 連続ブレンドの正準性 / ブレンド係数は乱流解像物理で駆動（衝撃センサと直交）/ SU2 `ROE_LOW_DISSIPATION` は ROE・SLAU(2) 限定 / σ は散逸テンソル乗数 / `ENUM_ROELOWDISS` の 5 値 / BSCW DDES の σ=σ_D+σ_NTS−σ_D·σ_NTS / 原 Ducros の低渦度過敏 / 改良 Ducros / 過散逸が解像乱流統計を破壊 / DLR LD2 が非構造 FV で成熟 / Sciacovelli 9次中心+選択的 shock detector。
- 2-1 確認: SU2 `FD` の定義 σ_FD=max(0.05,1-f_d)。
- **0-3 棄却**: 「SU2 の機構は『正準 DES がシールド関数 **と** Ducros の両方を必須とする』ことを**証明**する」── SU2 は両者を**層化して提供**するが `FD` 単体も正規の標準オプション。Ducros は実衝撃共存時のベストプラクティスであって**厳密必須ではない**。

---

## 6. 不確実な点・追加で読むべき一次文献

1. **`KEEP_SLAU` の `duc` の入り方**: SU2 同様「散逸テンソルに σ 乗算」か、「フルフラックスを中心↔風上で線形ブレンド」か。前者なら SU2 σ 式を直接移植可。→ [convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) `KEEP_SLAU_d`（L2210–, blend L2452–2456）を実確認。
2. **SBLI/ノズル/ピントルで `FD` 単体で足りるか**: SU2 BSCW は共存ゆえ `NTS_DUCROS`。forge は SU2 クロスチェック（[`.github/forge-su2-cross-check.md`](../../.github/forge-su2-cross-check.md)）で判定。
3. **`ψ_limiter` はセンサとして反応するか**: <0.1% 不変なら `1-ψ_limiter` 項は無寄与で `duc=1-f_d` に縮退。要計測。
4. **IDDES（vs DDES）での flux 係数**: f_d（DDES）か IDDES の f̃/f_dt か。本調査の検証ソース範囲外 → **Shur, Spalart, Strelets, Travin (2008, IJHFF, [doi:10.1016/j.ijheatfluidflow.2008.07.001](https://doi.org/10.1016/j.ijheatfluidflow.2008.07.001))** / **Gritskevich et al. (2012)** を一次で読む。
5. 全文 403 だった AIAA URL（2017-4284, J054956）は SU2 master ソースで内容裏取り済（同等以上の証拠）。

---

## 変更ログ

- `2026-06-22`: deep-research（98 エージェント / 16 ソース / 25 claim 検証）で新規作成。結論を [`turbulence-iddes-sst.md`](../../design/active/turbulence-iddes-sst.md) §4.8 に反映（短期=f_d 主導・limiter は positivity フロア / 中期=改良 Ducros を OR ブレンド追加）。

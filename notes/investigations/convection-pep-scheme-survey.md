# PEP (Pressure-Equilibrium-Preserving) 系スキーム技術調査と forge 実装方針

<!-- ファイル名規約: convection-<short-slug>.md -->

## メタ

- **area**: `convection`
- **status**: `draft` (**調査専用 / コード変更なし**)
- **related_docs**:
  - [`docs/convection/theory.md`](../../docs/convection/theory.md)
  - [`docs/convection/implementation.md`](../../docs/convection/implementation.md)
  - [`docs/thermophysics/theory.md`](../../docs/thermophysics/theory.md)
- **related_plans**:
  - [`convection-multispecies-contact-pressure.md`](../../design/active/convection-multispecies-contact-pressure.md) — 多成分 TP 接触面 limit-cycle の 1D ゲート切り分け (本調査の直接の動機)
  - [`su2-nemo-contact-thermo-investigation.md`](su2-nemo-contact-thermo-investigation.md) — SU2-NEMO の contact 処理棚卸し (double-flux/PEP 補正は SU2 にも無いと確認済み)
  - [`time_integration-general-eos-jacobian.md`](../../design/accepted/time_integration-general-eos-jacobian.md) — 一般 EOS 固有系 Jacobian (PEP の陰解法側で再利用)
  - [`thermophysics-multicomponent-tpgas.md`](../../design/accepted/thermophysics-multicomponent-tpgas.md) — NASA-9 TP gas 基盤 (face 内部エネルギー補正の素材)
- **created**: `2026-06-20`
- **owner**: (未定)

> **位置づけ**: 本文書は実装計画ではなく **技術調査レポート**。forge への PEP 系スキーム導入を判断するための材料 (文献調査結果 + forge 現状コードとの突き合わせ + 筆頭候補とリスク) をまとめる。実装に進む場合は、本調査を基に [`convection-multispecies-contact-pressure.md`](../../design/active/convection-multispecies-contact-pressure.md) の Gate 判定 (Gate 2 = 保存形 PEP 不足) を確定させてから、別 plan を `_template.md` で起こすこと。

---

## 0. エグゼクティブサマリ (結論先出し)

- **PEP 系は 2 系統**あり、**2023〜2026 年に融合**した。
  - **(A) 多成分・接触面系** (double-flux: Abgrall & Karni / quasi-conservative WENO: Johnsen & Colonius / DG-PEP: Ching et al.) — γ・組成ジャンプで保存形風上/DG が出す **spurious pressure oscillation** を抑える。
  - **(B) KEEP / split convective form 系** (Kuya-Totani-Kawai / Shima-Kuya-Wada / Jain-Moin) — 運動エネルギー・エントロピー保存の低散逸中心型で、分割形が **PEP 特性**も満たす。
  - 最新の **完全保存 PEP** (Coppola et al. 2026 / DeGrendele et al. 2025) は (B) の split-form を土台に、**質量流束に EOS 微係数 $\alpha=(\partial(\rho e)/\partial\rho)_p$ で結びつけた内部エネルギー流束補正**を載せ、(A) の圧力平衡を **厳密保存のまま** 達成する。多項式 $c_v(T)$ + 実在気体 departure / 任意 EOS・任意化学種数に対応 = **forge の NASA 多項式 TP と素直に整合**。
- **forge にとって最重要の評価軸 = 非構造 FV で成立するか**。PEP/KEEP のほぼ全研究は **構造格子 FD か高次 DG**。例外は 2 件のみ:
  - **Kuya-Okumura-Sawada (JCP 2023)**: 真に **非構造 FV** の KEEP (ただし **cell-vertex**, 多成分 TP-PEP ではない)。
  - **Badrkhani-Karpowski-Hasse (2025)**: **OpenFOAM (非構造 cell-centered FV)** の Entropy-Stable + Double-Flux (ただし double-flux 由来で **非保存**)。
- **筆頭候補**: **完全保存 PEP の「face 内部エネルギー流束補正」(APEP/APEC 型) を、forge 既存の `KEEP` convMethod に接ぎ木する**。
  - 理由: **厳密保存を維持** (衝撃捕捉・flux Jacobian と両立)、**EOS 一般** (TP/NASA 多項式ネイティブ)、**非保存な cell-state freezing 不要**。
  - 最大の注意点: **完全保存 PEP の厳密導出は構造 FD**。forge の **非構造 cell-centered FV へ face 補正を再導出**する必要がある (テンプレは Kuya 2023)。加えて中心型ゆえ **衝撃散逸/ハイブリッド**と、**陰解法 LHS のロバスト性**は文献に存在せず外挿。
- **重要な釘刺し** (Ranocha & Gassner 2021): **PEP ≠ 線形安定**。圧力振動を消しても entropy-conserving 系の局所線形不安定 (anti-diffusion) は残る。**PEP 単体では blockDPLUR は安定化しない** — 明示的散逸/limiter と LHS 設計は別途必須。

---

## 1. 背景: なぜ PEP が要るのか (forge の動機)

多成分/実在気体流で、**比熱比 $\gamma$ や組成が異なる気体の接触不連続・拡散界面**を、**完全保存形の風上/DG** で離散化すると、**spurious な圧力 (および速度) 振動**が発生する。これは数値スキームの欠陥であり、物理ではない。

- 根本原因 (Karni 1994, Abgrall 1996): $\gamma$ が face を跨いで変化するとき、**エネルギー方程式 / EOS 反転の離散熱力学的不整合**が生じる。
- 理論的帰結 (Ching-Johnson-Kercher, arXiv 2501.12532): **保存形 DG は velocity-equilibrium-preserving だが pressure-equilibrium-preserving ではない**ことが証明されている。

**forge での現れ方** (→ [`convection-multispecies-contact-pressure.md`](../../design/active/convection-multispecies-contact-pressure.md)):

- forge の SLAU は現状 **mixed-order face state**: $\rho, p, u$ は 2 次 MUSCL、組成 $Y$ はセル 1 次。face 温度を $T_f = P_f^{(2)} / (\rho_f^{(2)} R(Y_{\text{cell}}^{(1)}))$ と次数混在で再回復している ([convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) L60-70)。
- He/空気接触面の高 CFL 擬似時間で残差床 ~4e-6 の **limit-cycle** が出る。これが「mixed-order 由来 (Gate 1)」か「保存形 SLAU の PEP 本質不足 (Gate 2)」かを 1D ベンチで切り分ける、というのが現行 plan の主題。
- **本調査は Gate 2 (完全 PEP が要る) と確定した場合の選択肢を先回りで整理するもの。**

---

## 2. 主要スキーム比較表

非構造対応を最優先軸に、forge の評価軸 (密度ベース陽/陰・多成分 TP 整合・保存性・衝撃捕捉・実装コスト) で整理する。

| スキーム / 系統 | 発案者・年 | 元の枠組み | 非構造 FV | 保存性 | 多成分 TP | 衝撃捕捉 | 実装侵襲度 | forge 適性 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **double-flux** | Abgrall 1996 / Abgrall & Karni 2001 / Billet & Abgrall | 構造 FV/FD | △ (理論上可) | **× 非保存** (cell 内 $\gamma^*,\rho e_0$ 凍結) | 〇 (γ 凍結で対応) | 〇 (風上ベース) | 中 | △ 保存性放棄が難点 |
| **quasi-conservative WENO/HLLC** | Johnsen & Colonius 2006 (JCP 219) | 構造 高次 FV | △ | **準保存** (γ 量に非保存移流式追加, 誤差「極小」) | 〇 | 〇 (WENO) | 高 (WENO) | △ 高次 WENO は forge 非親和 |
| **KEEP (split form, FD)** | Kuya-Totani-Kawai 2018 (JCP 375) | **構造 FD** | × | 〇 厳密 | 拡張可 | × 中心型→要ハイブリッド | 中 | 〇 理論的土台 (forge に既存 `KEEP`) |
| **PEP split form** | Shima-Kuya-Wada 2021 | 構造 FD | × | 〇 | 拡張可 | × 中心型 | 中 | 〇 PEP 条件の出典 |
| **KEEP-FV (非構造)** | **Kuya-Okumura-Sawada 2023 (JCP)** | **非構造 FV** | **◎ 実証** (cell-vertex) | 〇 | (TP-PEP 明示せず) | × 中心型→要補正 | 中〜高 | ◎ **非構造化のテンプレ** (ただし cell-vertex) |
| **Entropy-Stable + Double-Flux** | **Badrkhani-Karpowski-Hasse 2025** (arXiv 2506.13231) | **OpenFOAM 非構造 FV** | **◎ 実装済** | **× 非保存** (double-flux 継承) | 〇 (γ, $H_t$ を補助変数) | 〇 (ES + 風上) | 中 | △ 非構造実証だが非保存 |
| **DG-PEP (pressure-evolution)** | Ching-Johnson-Kercher 2025 (arXiv 2501.12532 / AIAA SciTech) | 高次 DG | × | 半離散保存 (Abgrall 補正) | 〇 | 〇 (DG+limiter) | 高 | △ DG 前提で forge 非親和 |
| **完全保存 PEP (EPEP/APEP-RG)** | **Coppola-Aiello-De Michele 2026** (arXiv 2605.03617) | **構造 FD** | × (要再導出) | **◎ 厳密保存 + PEP** | **◎** ($c_v(T)$ 多項式 + 実在気体 departure) | × 中心型→要ハイブリッド | 中 | **◎ 筆頭候補の理論** |
| **完全保存 PEP + 高次 APEC** | **DeGrendele et al. (NASA Ames) 2025** (arXiv 2512.04450) | **構造 FD** | × (要再導出) | **◎ 厳密保存 + PEP** | **◎ 任意 EOS・任意種数** | × 中心型→要ハイブリッド | 中 | **◎ 筆頭候補の理論** |
| **two-phase 5-equation KEEP** | Jain & Moin 2022 (JCP) | 構造 FD | × | 〇 (相毎質量/運動量/全エネ) | (2 相 stiffened-gas, TP 非直結) | 〇 (diffuse interface) | 高 | △ 機構は参考、熱力学は別物 |

凡例: ◎=強く該当 / 〇=該当 / △=条件付き / ×=非該当。

---

## 3. 各系統の理論的要点

### 3.1 (A) 多成分・接触面系: 保存性 vs 圧力振動抑制のトレードオフ

- **double-flux** (Abgrall & Karni): 界面の圧力・速度が保たれる **必要十分条件**は、各セル内で **$\gamma^*$ と $\rho e_0$ を時間ステップ間で一定に保つ**こと。このため face flux を左右セルの EOS で **2 回計算** → 項がテレスコープせず **厳密保存を失う**。「圧力振動抑制のために保存性を犠牲にする」という古典トレードオフの原型。
- **quasi-conservative WENO** (Johnsen & Colonius): 保存 Euler 系に **γ/EOS 量の非保存移流式を 1 本追加** → 保存誤差は「極めて小さい」が厳密ではない。
- **DG-PEP** (Ching et al.): 圧力発展方程式を離散化し、Abgrall 型補正項で半離散エネルギー保存を回復。ただし厳密 PEP と厳密保存の両立にはさらなる調整が要る。

### 3.2 (B) KEEP / split convective form 系: 中心型・低散逸の理論的背骨

- **Kuya-Totani-Kawai 2018**: 質量・運動量の対流項を **split convective form** に書き換えると、**支配方程式間の解析的関係**が全エネルギー流束を **決定**する。この解析関係を **離散レベルで満たす**ことが、運動エネルギー↔内部エネルギー交換の正しさ = エントロピー保存に必須。
- **Shima-Kuya-Wada 2021**: どの分割形・どの平均が **追加で圧力振動を防ぐ (PEP 特性)** かを示した。**= forge が参照すべき PEP 条件の一次出典**。
- KEEP/PEP は **中心型・非散逸** → 衝撃には **人工散逸 / WENO / 風上ハイブリッド**が必須。

### 3.3 融合: 完全保存 PEP の本質 (筆頭候補の中身)

2025〜2026 の到達点は、**(B) の split-form を土台に (A) の圧力平衡を厳密保存のまま実現**:

- **PEP + 保存の同時条件** (Coppola et al.): 質量流束 $\delta F_\rho$ と内部エネルギー流束 $\delta F_{\rho e}$ を
  $$\delta F_{\rho e} = \alpha\,\delta F_{\rho}, \qquad \alpha = \left(\frac{\partial(\rho e)}{\partial \rho}\right)_{p}$$
  で結びつける。**EPEP-RG** は厳密な (一様状態で特異化しやすい) 平均を、**APEP-RG** は算術平均で代替し特異性を回避 — **両者とも全エネルギー厳密保存を維持**。
- **任意 EOS への一般化** (DeGrendele et al., APEP/高次 APEC): face の $\rho e$ を **中心平均 + 補正 $\beta$** で評価。$\beta$ は **face 組成 $\rho Y_i$ に関する $\rho e$ の偏微分 $\epsilon_i=(\partial \rho e/\partial \rho Y_i)_p$** を使う 2×2 (Moore-Penrose) 系から決まる。**非保存更新・過剰決定・EOS 固有設計が一切不要。**
- **TP/実在気体の扱い** (Coppola et al. §4): $c_v(T)=\sum_{k} c_k T^k$ (5/7/9 係数 = **NASA 多項式と同クラス**)、$e(\rho,T)=e_{TP}(T)+\Delta e(\rho,T)$ (van der Waals / Peng-Robinson departure)。
- **forge への含意**: 「**組成・温度・内部エネルギーの face 評価をどう組むか**が完全平衡保持を決める」という forge の問い (→ §1) に **直接答える**。鍵は **face の $\rho e$ を組成整合に補正評価する**こと。

### 3.4 反証された通説 (採用しないこと)

- 「**PEP のために再構成を primitive 変数で行うのが設計の鍵**」という主張は **0-3 で棄却**。primitive 再構成は PEP の機構ではない。**真の機構は EOS 微係数で結んだ内部エネルギー↔質量流束の結合** (§3.3)。

---

## 4. forge 現状コードとの突き合わせ

(forge コード調査の結論。所在は当該ファイル参照。)

### 4.1 再利用できる既存資産

1. **既存 `KEEP` convMethod の器**: [convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) に `KEEP_FVS_d` (L1174-)・`KEEP_SLAU_d` (L2210-) が実装済 (現状ラッパで無効化中)。`KEEP_SLAU_d` は `Ktilde` (運動エネ) ・`Itilde` (内部エネ) を持つ中心型 + SLAU 散逸ブレンド = **完全保存 PEP 補正を載せる自然な母体**。
2. **NASA-9 TP 熱力学基盤**: [thermo_d.cu](../../solver_density_cuda/cuda_forge/thermo_d.cu) に多項式 $h(T),c_p(T)$・Newton 温度反転・混合則。**face の $\rho e$ 補正と $\epsilon_i=(\partial\rho e/\partial\rho Y_i)_p$ 評価の素材が揃っている**。
3. **per-cell 物性配列** (`gamma[ic]`, `cp[ic]`, `Rmix[ic]`): $\alpha,\epsilon_i$ の face 評価に即流用可能。
4. **一般 EOS 固有系 Jacobian** ([time_integration-general-eos-jacobian.md](../../design/accepted/time_integration-general-eos-jacobian.md)): block-DPLUR が CPG 専用 EOS 仮定を脱し、実 $H_t,\kappa,\chi$ で TP 整合済 → **PEP の陰解法側の足場**。frozen-coefficient なので **新フラックス追加で LHS を作り直さず defect-correction で対応可能**。
5. **face state 再構成パイプライン**: 勾配→limiter→MUSCL が段階化済。species 2 次化はこのパイプラインへの参入で済む。

### 4.2 障壁・ギャップ

1. **非構造 FV への再導出 (最大リスク)**: 完全保存 PEP (Coppola/DeGrendele) は **構造 FD/method-of-lines**。face 内部エネルギー補正を **forge の非構造 cell-centered FV (任意多面体)** に **再導出**する必要があり、文献に前例なし。**テンプレは Kuya-Okumura-Sawada 2023 (非構造 KEEP-FV)** だが、それは **cell-vertex** を採用しており、forge の cell-centered とは離散点が異なる (cell-vertex は KEEP 性・2 次精度維持のため。cell-centered は安定性・精度が劣ると報告)。→ **cell-centered で成立させられるか / cell-vertex 化が要るか**は未解決の核心。
2. **陰解法 blockDPLUR との相互作用**: PEP 文献は **全て陽解法 RK**。$\alpha$ 結合・2×2 $\beta$ 補正の **flux Jacobian / 対角優位性 / 高 `cfl_pseudo` 安定性は完全に未文書**。forge では外挿になる。
3. **GPU コスト**: face 毎の小規模線形解 (2×2 / pseudoinverse) が、既存 Roe/SLAU/KEEP フラックスに対し **大規模 GPU でどれだけのオーバーヘッド**になるか要計測。
4. **衝撃ハイブリッド**: 中心型 PEP/KEEP を forge 既存の風上 (Roe/SLAU) と **どう切り替えるか** (平滑/界面=PEP-KEEP、衝撃=SLAU/Roe)。既存 Ducros センサ ([convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) L189) が下地。
5. **線形安定 ≠ PEP** (Ranocha-Gassner 2021): PEP 化しても entropy-conserving 系の **局所線形不安定 (密度波の anti-diffusion) は残る**。明示散逸/limiter と LHS 設計は別途必須。

---

## 5. forge への筆頭候補と推奨

### 5.1 筆頭候補

**完全保存 PEP の face 内部エネルギー流束補正 (APEP-RG / 高次 APEC 型) を、既存 `KEEP` convMethod に接ぎ木する。**

採用理由 (forge の制約に対し):

- **非構造**: 完全保存 PEP 自体は構造 FD だが、**Kuya 2023 の非構造 KEEP-FV を再導出テンプレ**にできる。double-flux 系 (非保存) や DG 系 (高次 DG 前提) より forge の有限体積骨格に近い。
- **密度ベース陽/陰**: **厳密保存** = 既存 flux Jacobian・defect-correction 陰解法と整合。非保存な cell-state freezing が無いので blockDPLUR の保存則整合を壊さない。
- **多成分 TP**: $c_v(T)$ 多項式 = NASA-9 と同クラス、任意種数対応。forge の [thermo_d.cu](../../solver_density_cuda/cuda_forge/thermo_d.cu) 基盤に直接乗る。
- **保存性**: 厳密保存を維持したまま PEP → 衝撃捕捉とも両立。

### 5.2 次点・比較対象

- **double-flux (OpenFOAM 実証ライン)**: 非構造 cell-centered FV での **動作実証**としては最も近い (Badrkhani 2025) が、**非保存**を受け入れる必要。Gate 2 確定後、完全保存 PEP の非構造再導出が難航した場合の **フォールバック**として位置づけ。
- **非構造 KEEP-FV (Kuya 2023) そのもの**: cell-vertex 前提。forge の median-dual 化 ([discretization-median-dual.md](../../design/active/discretization-median-dual.md)) が進めば node-centered = cell-vertex 的な土俵になり、**親和性が上がる**可能性 → 中長期で再評価の価値。

### 5.3 推奨する段階導入の道筋 (実装着手時)

1. **Gate 確定が先**: [`convection-multispecies-contact-pressure.md`](../../design/active/convection-multispecies-contact-pressure.md) の 1D He/空気ベンチで「Gate 1 (mixed-order 再構成主因) か Gate 2 (保存形 PEP 不足)」を確定。Gate 1 なら **species 同次再構成 (低コスト)** で済み、本格 PEP は不要。
2. **Gate 2 確定なら本 plan を `_template.md` で起こす**。まず **1D 構造的セットアップで APEP face $\rho e$ 補正を検証** (Coppola/DeGrendele 再現)。
3. **非構造 cell-centered での face $\rho e$ 補正を再導出** (Kuya 2023 を参照)。陽解法 RK で 1D/2D 接触面の圧力振動消失と保存則を確認。
4. **衝撃ハイブリッド** (既存 Ducros + SLAU/Roe フォールバック) を組む。
5. **blockDPLUR との整合**を defect-correction で確認 (frozen-coefficient、収束解不変を回帰基準)。
6. cell-centered で精度/安定が出ない場合、**cell-vertex (median-dual) ルート**を [discretization-median-dual.md](../../design/active/discretization-median-dual.md) と合流して再検討。

---

## 6. open questions (文献で未解決 = forge 側で詰める論点)

1. 完全保存 PEP の face 内部エネルギー補正を、**非構造 cell-centered FV (任意多面体)** で一貫導出できるか。cell-vertex (Kuya 2023) が必須なら、その精度/ロバスト性ペナルティは?
2. **blockDPLUR LHS** との相互作用 ($\alpha=(\partial\rho e/\partial\rho)_p$ 結合・2×2 $\beta$ 補正の Jacobian、対角優位性、高 `cfl_pseudo` 条件数)。全文献が陽解法 RK のため未文書。
3. forge の多成分 TP (per-species NASA 多項式・温度依存 $c_p$・混合 EOS) における **$\epsilon_i=(\partial\rho e/\partial\rho Y_i)_p$ の最も簡潔な厳密形**と、その face 毎コスト (小線形解) の GPU スケール比較。
4. 中心 PEP/KEEP と既存風上 (Roe/SLAU) の **最良ハイブリッド戦略** (Ranocha-Gassner の線形不安定を踏まえ、平滑/界面 PEP・衝撃 SLAU の切替閾値)。

---

## 6.5 検証ケースの提案 (1D → 2D ラダー)

PEP スキームは「**一様圧力を一様に保てるか (PEP)**」「**厳密保存**」「**衝撃捕捉との両立**」「**線形安定 (Ranocha-Gassner)**」「**非構造で成立するか**」を**独立に**潰す必要がある。1 ケースで全部見ようとせず、**易しい素過程から**積み上げる。各ケースは forge 規律に従う:

- 投入時・まとめ時に `run_*` パスをリポジトリルートからの相対で明示、各 case の `README.md` に「## 計算 run 一覧」表を維持。
- 全 run で **step 0 から NaN/Inf 監視** (`detectNaN=1`)、`residual_history.png` を残す。
- 定常を主張する場合は `solver_density_cuda/tools/check_convergence.py` の VERDICT を必ず添える。非定常は「発達/定常化」を中間 `res_*.h5` で確認 (→ [develop-flow-before-reporting])。
- **回帰基準**: PEP off (`convMethod=SLAU` 等) で既存ケースとビット一致を保つ (opt-in 実装)。

### ラダー一覧

| # | 次元 | ケース | 主眼 (何を分離して見るか) | 既存 case 対応 | 合否ライン |
| --- | --- | --- | --- | --- | --- |
| **V0** | 1D | **静止等温 材料接触面** (He\|N2, 一様 $p$, $u{=}0$) | **PEP の素過程・ゲート** | 新規 (or [`convection-multispecies-contact-pressure.md`](../../design/active/convection-multispecies-contact-pressure.md) の 1D) | $\max\|p-p_0\|/p_0$: PEP $\lesssim10^{-12}$ / 現行 SLAU は有限振動 |
| **V1** | 1D | **移流 材料接触面** (低マッハ $M{\sim}0.01$–$0.1$) | 低マッハ PEP・`lowMachPrecond` 相互作用 | 同上を $u{\neq}0$ 化 | 圧力振幅 vs $M$ が増大しない |
| **V2** | 1D | **多成分 shock tube** (He/N2 両側異種) | 衝撃+接触+膨張、**保存性**、衝撃ハイブリッド | [`case/05.sod_shock_tube`](../../case/05.sod_shock_tube/) (既存 `*_slau_tp_n2_*` 拡張) | 厳密 Riemann/CEA と接触圧 overshoot 比較、$\sum$ 保存量 機械精度 |
| **V3** | 1D | **密度/組成波 伝播** (一様 $p,u$, 正弦波長時間移流) | **線形安定** (entropy-cons の anti-diffusion) | 新規 (周期 1D) | 振幅が成長しない (制御散逸内) |
| **V4** | 1D | **厳密保存検算** (任意 Riemann の総和監視) | テレスコープ=機械精度保存 | V2 と同 run で計測 | $\frac{d}{dt}\sum(\rho,\rho u,\rho E)\approx0$ ($10^{-13}$) |
| **V5** | 2D | **斜め移流 材料界面** (周期箱, 格子非整列, 三角/歪みメッシュ) | **非構造 cell-centered で PEP 成立か** (§6-1 の核心) | 新規 ([`discretization-median-dual`](../../design/active/discretization-median-dual.md) のメッシュ流用可) | $\max\|p-p_0\|$ 機械精度、界面が歪まず移流 |
| **V6** | 2D | **shock-bubble interaction** (Haas-Sturtevant: He or R22 バブル) | 多成分+衝撃+界面不安定の総合 | 新規 (or [`case/06`](../../case/06.mach3_wind_tunnel/) 系) | 実験/既往計算と界面形状・波系を定性照合 |
| **V7** | 2D | **Taylor-Green vortex** (非粘性, 単/多成分) | **KEEP 性=運動エネ保存・低散逸** | 新規 (周期箱) | 非粘性で KE 散逸が upwind より桁で小 |
| **V8** | 2D | (発展) **多成分 Riemann / RM 不安定** | 実運用近接の総合 | 新規 | 既往参照解と定性一致 |

### 各ケースの要点

**V0 (最優先・ゲート)**: 一様 $p$・$u{=}0$ で中央に組成ジャンプ (He\|N2 など TP)。温度は $p$ と組成から決まり**密度ジャンプを伴う**。厳密 PEP なら $p,u$ は機械精度で不変、非 PEP は接触面で spurious 圧力・速度振動。**これが [`convection-multispecies-contact-pressure.md`](../../design/active/convection-multispecies-contact-pressure.md) の Gate 1 (mixed-order 再構成主因) / Gate 2 (保存形 PEP 不足) を切り分ける実体**。現行 `KEEP`/`SLAU` を baseline に、§7.5 の `Itilde` 差し替え版で振動が機械精度まで消えることを確認する最初の関門。

**V1**: V0 を一定速度で移流。低マッハ漸近で圧力振動が出ないか、既存 `lowMachPrecond` (`c'` 散逸前処理) と干渉しないかを見る。[`time_integration-lowmach-preconditioning.md`](../../design/accepted/time_integration-lowmach-preconditioning.md) の知見と接続。

**V2**: 両側異種ガスの shock tube。**衝撃捕捉ハイブリッド** (平滑/界面=PEP, 衝撃=既存 `duc`-blend で SLAU) が、接触面の偽圧力を消しつつ衝撃で振動しないかを検証。厳密 Riemann (CEA 組成) と照合。既存 `case/05` に `*_slau_tp_n2_*` の run があるので拡張が早い。

**V3 (Ranocha-Gassner 釘刺しの実体)**: 一様 $p,u$ で密度 or 組成の正弦波を長時間移流。entropy-conserving 中心流束の **anti-diffusion 由来の局所線形不安定**が出ないかを確認。**PEP 化しても安定化しない**ことの実証で、必要散逸量の下限を決める。

**V5 (非構造の要)**: 周期箱で組成界面を**格子に対し斜め**に一定速度移流。**三角形・歪み・非直交メッシュ**で実施し、構造 FD 由来の PEP 補正が **forge の非構造 cell-centered FV で成立するか** (§6-1) を直接確認。cell-centered で精度/PEP が崩れるなら cell-vertex (median-dual) ルート判断材料。forge の free-stream 非保存問題 ([forge-freestream-nonorthogonal]) とも切り分ける。

**V7 (KEEP 性)**: 非粘性 Taylor-Green vortex で運動エネルギー時間履歴を見て、中心型 PEP-KEEP の低散逸性を upwind (Roe/SLAU) と定量比較。多成分版は組成パッチを入れて PEP と KEEP を同時に見る。

### 推奨実施順序

`V0 → (V4 保存検算は V2 に同梱) → V2 → V3 → V1` で 1D を固め、**V0 で Gate 2 が確定**したら 2D へ。2D は **V5 (非構造 PEP) を最優先**、続いて V7 (KEEP 性) → V6 (総合) → V8。
V0・V5 が PEP 実装の **可否を決める 2 大関門**、V3 が **安定散逸量の下限**、V2・V6 が **衝撃両立**、V7 が **低散逸の御利益**を担保する。

---

## 7.5 付録: 完全保存 PEP face 補正の具体的な中身

> 原著 (Coppola arXiv:2605.03617 / DeGrendele arXiv:2512.04450) の式を抽出し、forge 実装に落ちる形で書く。
> **記号の衝突に注意**: Coppola の $\alpha=(\partial\rho e/\partial\rho)_p$ (EOS 微係数) と DeGrendele の $\alpha$ (組成再配分の補正スカラ) は **別物**。以下では Coppola 系を $\alpha_{\rho},\lambda$、DeGrendele 系を $\epsilon_i,(\alpha_Y,\beta)$ と書き分ける。

### A. 全体像 (一言で)

**「中心型 KEEP 流束 + face 内部エネルギーの EOS 整合補正」**。upwind 流束ではない。

- 質量・運動量・運動エネルギーは **split form (算術平均積)** の中心流束。
- **内部エネルギー流束だけ**に、「**face を跨いで変化する保存量 (密度 or 組成) に対する $\rho e$ の感度**」で決まる補正を足す。
- この補正で「一様圧力なら離散圧力も一様に保たれる ($\mathrm{d}p/\mathrm{d}t=0$)」= **PEP** が、**厳密保存 (流束テレスコープ) を壊さずに**成立する。

### B. Coppola APEP-RG (単成分・実在/TP gas, 最も簡潔な形)

セル $i,i{+}1$ 間の face。バー $\bar a=(a_i+a_{i+1})/2$ は算術平均。

**EOS 微係数 (各セルで評価):**
$$\lambda=\left(\frac{\partial e}{\partial\rho}\right)_p,\qquad \alpha_\rho=\left(\frac{\partial(\rho e)}{\partial\rho}\right)_p=e+\rho\lambda$$

**流束:**
$$\mathcal F_\rho=\bar\rho\,\bar u \quad(\text{質量}),\qquad
\mathcal F_{\rho e}=\bar{\alpha_\rho}\,\mathcal F_\rho-\bar u\,\overline{\rho^2\lambda}\quad(\text{内部エネ, PEP 補正本体})$$
$$\mathcal F_{\rho u}^{tot}=\mathcal F_\rho\bar u+\bar p,\qquad
\mathcal F_{\rho e}^{tot}=\mathcal F_{\rho e}+\mathcal F_\rho\frac{u_i u_{i+1}}{2}+\overline{up}$$
(運動エネは KEEP 流儀で $u_i u_{i+1}$ の **積**を使う点が肝。$\overline{up}=(u_ip_{i+1}+p_iu_{i+1})/2$。)

**CPG 健全性チェック**: CPG では $\rho e=p/(\gamma-1)$ が $\rho$ に依らない → $\alpha_\rho=0$, $\rho^2\lambda=-p/(\gamma-1)$ → $\mathcal F_{\rho e}=\bar u\,\overline{p/(\gamma-1)}$ = 通常の内部エネ移流。**$\alpha_\rho$ 補正が「$\gamma$/組成が face で変わるとき」だけ効く。**

- **EPEP-RG (厳密)** は $\mathcal F_\rho=\bar\rho^\lambda\bar u$ に特殊平均 $\bar\rho^\lambda=\delta^+(\rho^2\lambda)/\delta^+\alpha_\rho$ を使うが、**一様状態で $\delta^+\alpha_\rho\to0$ → 特異**。
- **APEP-RG (実用)** はこれを **算術平均 $\bar\rho$ に置換**して特異性回避。厳密 PEP ではないが「全シミュレーションで使える近似 PEP」。**forge はこちらを採る。**
- **保存性**: 流束差分 $\mathrm d(\rho E)_i/\mathrm dt=-\frac1h\delta^-\mathcal F_{\rho e}^{tot}$ がテレスコープ → 質量・運動量・**全エネルギーを機械精度で厳密保存**。

### C. DeGrendele APEC/PEP (多成分 = forge が要る形)

face 値そのものを「中心平均 + 補正」にする。2 成分例 (Eq.75-77):
$$\rho Y_1|_{f}=\overline{\rho Y_1}+\alpha_Y,\quad \rho Y_2|_{f}=\overline{\rho Y_2}-\alpha_Y\quad(\text{組成を保存しつつ再配分})$$
$$\rho e|_{f}=\overline{\rho e}+\beta-\sum_i\Big[\epsilon_i|_{(m+1)/2}\,(\rho Y_i|_{m+1}-\rho Y_i|_f)-\epsilon_i|_{m/2}\,(\rho Y_i|_f-\rho Y_i|_m)\Big]$$

**組成感度 (各セルで EOS から):**
$$\epsilon_i=\left(\frac{\partial\rho e}{\partial(\rho Y_i)}\right)_{\rho Y_{j\ne i},\,p}$$

**離散 PEP 条件 (Eq.30, これを満たすよう $(\alpha_Y,\beta)$ を決める):**
$$\rho e|_{m+1/2}-\rho e|_{m-1/2}=\sum_i\left(\frac{\partial\rho e}{\partial\rho Y_i}\right)_p\big(\rho Y_i|_{m+1/2}-\rho Y_i|_{m-1/2}\big)$$

**$(\alpha_Y,\beta)$ の決め方**: 左右 face での PEP 条件 (Eq.32-33) から 2×2 線形系 (Eq.80) を立て、**Moore-Penrose 擬似逆**で解く。条件数の逆数 $r$ が閾値 $r_g$ 未満 (≈一様状態) なら **$\alpha_Y=\beta=0$** にフォールバック。

### D. forge の TP 混合での $\epsilon_i$ 具体形 (本調査で導出)

forge の NASA-9 TP 混合 $\rho e=\sum_k(\rho Y_k)e_k(T)$, 拘束 $p=\big(\sum_k\rho Y_k R_k\big)T$ ($R_k=R_u/W_k$) から、$c_i\equiv\rho Y_i$ として:
$$\epsilon_i=\left(\frac{\partial\rho e}{\partial c_i}\right)_p=e_i(T)-\frac{c_{v,\mathrm{mix}}}{R_{\mathrm{mix}}}\,R_i\,T$$
($e_i(T),c_{v,\mathrm{mix}}$ は NASA 多項式から、$T$ は face/セルの温度。)
**単成分極限で $\epsilon=e-c_vT=0$** → C と一貫。**forge の [thermo_d.cu](../../solver_density_cuda/cuda_forge/thermo_d.cu) で評価済みの $e_i(T),c_v,R$ から追加コストほぼ無し**で計算できる。datum (生成エンタルピ基準, `thermoHrefTemp`) を $e_i$ で一貫させること。

### E. forge カーネルへの落とし込み

| 項目 | 現状 ([convectiveFlux_d.cu](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu)) | PEP 化の差分 |
| --- | --- | --- |
| 質量流束 `Ctilde` | KEEP 中心 + SLAU 散逸 | そのまま (split form 中心) |
| 運動量 `Mtilde` | `Ctilde·ū + p̄` | そのまま |
| 運動エネ `Ktilde` | `Ctilde·½(u_L·u_R+…)` | そのまま (既に $u_Lu_R$ 積) |
| **内部エネ `Itilde`** | `Ctilde·½(P_L/ρ_L+P_R/ρ_R)/(γ−1)` (**単一 γ 仮定**) | **ここを PEP 補正形に差し替え**: B の $\bar{\alpha_\rho}\mathcal F_\rho-\bar u\,\overline{\rho^2\lambda}$ か C の補正 $\rho e\|_f$ |
| 化学種 face 値 | セル 1 次 | C の再配分補正 $\pm\alpha_Y$ (or まず同次 2 次) |
| EOS 微係数 | — | per-cell `alpha_rho[ic]`, `eps_i[ic]` を D で前計算 |
| 衝撃散逸 | `duc`-blend で KEEP↔SLAU | **そのまま流用** (平滑/界面=PEP, 衝撃=SLAU) |

**非構造化**: 1D の $\delta^\pm$ を face 法線流束 $\mathcal F\cdot\mathbf S$ に置換。split form 積 ($\bar\rho\bar U_n$, $u_Lu_R$, $\bar p$) と内部エネ補正は **2 セル状態のみのローカル量**なので cell-centered でそのまま試せる。ただし KEEP 性・2 次精度の厳密保証は Kuya 2023 が cell-vertex を要したので、cell-centered での精度/安定は要検証 (→ §4.2-1, §6-1)。

### F. 1 face あたりのアルゴリズム (擬似コード, 多成分)

```
入力: 隣接 2 セル状態 (ρ, ρY_i, u, p, T), face 法線 S
1. 各セルで EOS 微係数を評価:
     T_ic, c_v,mix, R_mix, R_i, e_i(T)  ← NASA 多項式 (既存 thermo_d)
     ε_i[ic] = e_i(T) − (c_v,mix/R_mix)·R_i·T        (式 D)
     (単成分なら α_ρ[ic] = e + ρλ で式 B 経路)
2. 算術平均: ρ̄, Ū_n=S·ū, p̄, (ρe)̄, (ρY_i)̄
3. PEP 補正を解く:
     2×2 系 (Eq.80) を ε_i の左右差から組み, 条件数 r を評価
     if r > r_g:  (α_Y, β) = pinv · rhs       else: α_Y=β=0  (フォールバック)
4. 補正 face 値: ρY_i|_f = (ρY_i)̄ ± α_Y,  ρe|_f = (ρe)̄ + β − Σ ε_i 補正項
5. 流束 (法線):
     F_ρ      = ρ̄ · Ū_n
     F_ρY_i   = (ρY_i|_f) · Ū_n
     F_ρu     = F_ρ · ū + p̄ · S
     F_ρe^tot = (ρe|_f)·Ū_n + F_ρ·(u_L·u_R/2) + (u p)̄·|S|
6. 衝撃センサ duc で SLAU 流束と blend (既存機構)
7. 両セルへ ±F を atomicAdd
```

**コスト**: 追加は per-cell の $\epsilon_i$ (NASA 評価は既存) と per-face の 2×2 擬似逆のみ。Roe/SLAU の 5×5 固有分解より軽い見込みだが GPU 実測要 (→ §6-3)。

---

## 7. 一次情報源 (引用)

すべて deep-research で fetch + 敵対的検証 (3 票) を通した一次文献 (24/25 claim confirmed)。

| # | 文献 | 種別 | 役割 |
| --- | --- | --- | --- |
| 1 | Johnsen & Colonius, *JCP* 219 (2006) — quasi-conservative WENO multicomponent [ResearchGate 222580824](https://www.researchgate.net/publication/222580824_Implementation_of_WENO_Schemes_in_Compressible_Multicomponent_Flow_Problems) | primary | 問題の定義・準保存 WENO |
| 2 | Ching, Johnson & Kercher (2025) — DG-PEP [arXiv:2501.12532](https://arxiv.org/abs/2501.12532) | primary (preprint+AIAA) | 保存形 DG が PEP でない証明・DG-PEP |
| 3 | Badrkhani, Karpowski & Hasse (2025) — Entropy-Stable/Double-Flux in OpenFOAM [arXiv:2506.13231](https://arxiv.org/html/2506.13231) | primary (preprint) | **非構造 FV での double-flux 実証** |
| 4 | Kuya, Totani & Kawai, *JCP* 375 (2018) — KEEP scheme [S0021999118305916](https://www.sciencedirect.com/science/article/abs/pii/S0021999118305916) | primary | KEEP / split convective form の原典 |
| 5 | Kuya, Okumura & Sawada, *JCP* (2023) — **KEEP FV on unstructured meshes** [S0021999123006162](https://www.sciencedirect.com/science/article/abs/pii/S0021999123006162) | primary | **非構造 KEEP-FV (再導出テンプレ, cell-vertex)** |
| 6 | Jain & Moin, *JCP* (2022) — KEEP two-phase 5-equation [S0021999122003692](https://www.sciencedirect.com/science/article/abs/pii/S0021999122003692) | primary | 保存 + KE 保存機構 (2 相, TP 非直結) |
| 7 | Ranocha & Gassner, *CAMC* (2021) — PEP ≠ linear stability [Springer s42967-021-00148-z](https://link.springer.com/article/10.1007/s42967-021-00148-z) | primary | **PEP 単体で安定化しない (釘刺し)** |
| 8 | Coppola, Aiello & De Michele (2026) — EPEP/APEP-RG 完全保存 PEP for real/TP gas [arXiv:2605.03617](https://arxiv.org/html/2605.03617v2) | primary (preprint) | **筆頭候補の理論 (TP/実在気体)** |
| 9 | DeGrendele et al. (NASA Ames, 2025) — 完全保存 PEP + 高次 APEC, 任意 EOS [arXiv:2512.04450](https://arxiv.org/html/2512.04450) | primary (preprint) | **筆頭候補の理論 (任意 EOS・任意種数)** |

**ソース品質の留保**:
- 筆頭候補の中核 (Coppola 2605.03617 / DeGrendele 2512.04450 / Badrkhani 2506.13231 / Ching 2501.12532) は **2025-2026 の arXiv preprint** で査読未確定のものを含む (Ching は AIAA SciTech 2025 に対応論文あり)。基盤系統 (Johnsen-Colonius 2006 / Kuya 2018 / Jain-Moin 2022 / Kuya 2023 / Ranocha-Gassner 2021) は **査読済 JCP/Springer**。
- **最大の未充足ギャップ**: 厳密 PEP 完全保存スキーム (Coppola, DeGrendele) は **すべて構造 FD 由来**で、非構造 FV 定式化を **提供していない**。非構造化は forge 側の主要リスク (→ §6-1)。
- Jain-Moin と 5-equation 系は **2 相 stiffened-gas** で、多成分 TP の drop-in ではない (機構は参考)。
- 陰解法整合は文献に皆無 (全て陽解法 RK)。

---

## 8. 変更ログ

- `2026-06-20` — 初稿 (調査専用)。deep-research (97 agents / 24 confirmed claims) + forge コード調査を統合。筆頭候補 = 完全保存 PEP の face 内部エネルギー補正を既存 `KEEP` に接ぎ木。実装着手は [`convection-multispecies-contact-pressure.md`](../../design/active/convection-multispecies-contact-pressure.md) の Gate 2 確定後。
- `2026-06-20` — §7.5 付録追加。Coppola APEP-RG / DeGrendele APEC の原著式を抽出し具体形を記載。forge TP 混合の $\epsilon_i=e_i(T)-(c_{v,\mathrm{mix}}/R_{\mathrm{mix}})R_iT$ を導出、`Itilde` 差し替え点・1 face アルゴリズム擬似コードを追記。
- `2026-06-20` — §6.5 検証ケース提案 (1D→2D ラダー V0–V8) 追加。V0 静止材料接触面 (PEP ゲート) ・V5 斜め移流非構造界面 (非構造 PEP の要) を 2 大関門、V3 を線形安定下限、V2/V6 を衝撃両立、V7 を KEEP 性として配置。既存 case 対応と実施順序を明記。

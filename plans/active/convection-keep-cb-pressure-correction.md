# KEEP: 高周波圧力欠陥駆動の mass-flux 補正 (Rhie–Chow 型市松キラー)

## メタ

- **area**: `convection`
- **status**: `draft`
- **related_docs**:
  - `methods/convection/implementation.md` (KEEP + ES 散逸・precond・fdblend/opBlend の現在仕様)
- **related_plans**:
  - `plans/active/convection-keep-diss-lowmach-precond.md` (Turkel 前処理散逸 = 本補正の既存版・σ_min 制限付き)
  - `plans/accepted/convection-keep-diss-recon-jump.md` (keepDissJump 再構成ジャンプ = δp^HF の抽出機構)
  - `plans/active/turbulence-iddes-sst.md` §4.8 (fdblend / opBlend / 圧力モード汚染の経緯)
- **created**: `2026-07-22`
- **owner**: sano + Claude

## 1. 目的

case/39 周期丘 DDES (本番格子, 丘頂 Δz/y1≈70 の強異方性帯) に残る **spanwise 2–3Δz
定在圧力鋸歯** (±150–300 Pa ≈ 0.1–0.3 q_crest, 定在・飽和・KEEP 系固有 = SLAU 比 200 倍,
run_0014 診断 2026-07-22) を、**解像渦の散逸 (σ_min=0.05) を一切変えずに**減衰させる。
手段は、再構成不整合として抽出した高周波圧力欠陥 δp^HF だけを入力とする独立な
mass-flux 補正レイヤ (圧力ベース法の Rhie–Chow / all-speed flux の圧力差項と同族)。

## 2. スコープ

- **やる**: KEEP_d カーネルへの追加補正項 (config で opt-in, 既定 off=ビット不変)、
  係数掃引による較正、丘頂鋸歯・物理場コスト・収束性の定量評価、採用時の本番反映。
- **やらない**:
  - 既存 precond 2×2 行列の σ_p/σ_u 列別スケール化 (**第二候補**として保留。列別係数は
    sym(K) 正定値性が自動保証されず、entropy quadratic form・M→0 漸近・M≥1 復帰の
    事前検証が必要なため、まず独立レイヤ案を試す)
  - 厳密 Rhie–Chow (運動量対角 a_P 由来の面係数)。dual-time の pseudo 状態
    (nStepInner/implicitRelax/LHS 近似) に空間 flux が依存してしまうため不採用
  - 方向別 (face metric 依存) 係数、LHS (block-DPLUR) への剛性反映 (効果不足時の将来課題)
  - SLAU / Roe 側への同補正 (自前の散逸で不要)

## 3. 関連 docs と前提

- 現状の KEEP 散逸レイヤ構成は [methods/convection/implementation.md](../../methods/convection/implementation.md):
  ES 行列散逸 (keepDissType=2) + 再構成ジャンプ (keepDissJump=2) + Turkel precond
  (keepDissPrecond=1) + f̃_d 駆動 σ (keepDissFdBlend=1) + opBlend。
- **既存 precond の Δp 散逸 (∝c²/Ur, 市松キラー) は実装済みだが σ_min=0.05 に制限**
  されている点が本計画の出発点。「Δp/ΔU の分離」を新設するのではなく、
  **圧力結合の強度だけを LES 用 σ_min から独立させる**追加改良である。
- 実装に着手する前に `methods/convection/implementation.md` へ本補正の仕様
  (§KEEP 散逸レイヤに追記) を反映する (開発フロー step 1)。
- 前提となる観測 (case/39 run_0014, 2026-07-22 診断):
  - 鋸歯は x/h≈8.3–8.7 帯に限定、HP 成分のスナップ間相関 0.70–0.84 (2–8 ms) = 定在
  - 50k step で振幅トレンドなし (飽和)、旧構成 run_0013 step5000 でも同振幅 = opBlend 起因ではない
  - SLAU RANS (run_0012) では 0.8 Pa と 1/200
  - 機構の安全な表現: 「Δz/y1≈70 の強い異方性 + 低散逸 KEEP/jump 再構成の組合せで、
    壁法線方向に拘束された圧力モードの spanwise 減衰が不足」(面積極小→結合弱い、の
    面積単独論は体積割りで相殺されるため因果としては不採用)

## 4. 設計方針

### 4.1 補正の定義

面 f (ic0→ic1) ごとに、**高周波圧力欠陥**を再構成不整合として抽出する:

$$p_L^f = p_0 + \nabla p_0\cdot(\boldsymbol x_f-\boldsymbol x_0),\quad
  p_R^f = p_1 + \nabla p_1\cdot(\boldsymbol x_f-\boldsymbol x_1),\quad
  \delta p_f^{HF} = p_R^f - p_L^f$$

- 線形圧力場では厳密ゼロ (勾配再構成が線形場に厳密) → **滑らかな Cp・free-stream を触らない**
- 2Δ 市松では中心勾配 ≈0 → 生の Δp がそのまま残る = 狙い撃ち
- keepDissJump≥1 の分岐で `pL`/`pR` は**既に計算済み** ([convectiveFlux_keep_d.inc.cuh:236-238](../../solver_density_cuda/cuda_forge/convection/convectiveFlux_keep_d.inc.cuh#L236))。
  sign-property クリップ (jump=2) は特性射影 rd_k に掛かるもので、δp^HF はクリップ前の
  素の pR−pL を使う (市松では両者一致)。
- **有効面の限定 (生 Δp へフォールバックしない)**: 補正は「keepDissJump≥1 かつ内部面かつ
  両側再構成が成功 (pL,pR,roL,roR>0)」の面のみ有効とし、それ以外 (ghost / 正値性 fallback)
  は **δp^HF=0** (補正なし)。生 Δp を使うと滑らかな物理圧力勾配まで拡散し市松狙い撃ちの
  設計原則と矛盾するため。対象の 2Δ 市松では再構成圧力は正なので効力は失わない。
  config 読込時に `keepDissCbCoeff>0 && keepDissJump<1` をハードエラーにする。
- **float32 ゲージ安全化 (実装時に発覚・確定)**: δp^HF を絶対圧力の pR−pL で組むと
  GG 勾配のメトリック閉合ノイズ (+p~1e5 の量子化) が 1/Ur 増幅と掛かり、**非直交 free-stream
  (case/33) を 5 step で破壊** (step0 roe 2e-2 = 基準の 6 万倍)。対策は 2 段:
  ① **差分形** dpRec = (Ps1−Ps0) + (g1·d1 − g0·d0) (量子化除去)、
  ② **sign-property minmod**: dpHF = minmod(dpRec, dpRaw) — 生 Δp と符号反転なら 0、一致なら
  小さい方。一様場では dpRaw が厳密 0 → 勾配ノイズを構造的に完全消去 (jump2 の rd クリップが
  free-stream 安全である理由と同一原理の移植)。市松は dpRec=dpRaw でフル振幅維持。
  検証: case/33 cell/node とも jump2 単独と同一水準 (roe ~1e-6, 300step 安定) に復帰、
  case/35 L1 減衰は不変 (cell C_cb=0.02 で 5 step 700 倍減衰)。

補正 mass flux (all-speed インピーダンス 1/(ρUr) 型、密度ベース版 Rhie–Chow 相当):

$$\delta\dot m_f = -\tfrac12\, C_{cb}\, S_f\, \frac{\delta p_f^{HF}}{U_{r,cb}},\qquad
  U_{r,cb} = c\,\min\!\bigl(1,\ \max(M,\ \epsilon_{cb})\bigr) = \min\!\bigl(c,\ \max(|\boldsymbol u|,\ \epsilon_{cb}c)\bigr)$$

M≥1 で Ur=c に上限 (all-speed 整合)。実装は既存 `lowMachUr()` を ε だけ変えて再利用する。
符号は「高圧側→低圧側へ質量を逃がす」向き (既存 precond d11 項と同符号系。実装時に
2 セル市松の単体検算で確認する)。

### 4.2 全方程式への整合適用 (mass flux にだけ定義し、全移流で共通利用)

密度残差だけを動かすと $p=(\gamma-1)(\rho E-\tfrac12\rho|u|^2)$ の熱力学整合が壊れ
市松が減衰しないため、補正 flux を面平均量で全保存量へ配る:

$$\delta\boldsymbol F_f = \delta\dot m_f\,[\,1,\ \bar u,\ \bar v,\ \bar w,\ \bar H_t\,]^T$$

- 面平均は既存の `uxF, uyF, uzF, Ht` (算術平均) をそのまま使う
- 丘頂の z 面では ū_z≈0 → 運動量への影響は自動的に微小、連続+エネルギーで圧力モードのみ減衰
- **スカラー (roK/roOmega/species/受動スカラー) は自動整合**: KEEP_d は
  `massflux[ip] = res_ro_temp` ([同:521](../../solver_density_cuda/cuda_forge/convection/convectiveFlux_keep_d.inc.cuh#L521))
  で散逸込み質量流束を書き出し、`scalarTransport_d` / `speciesTransport_d` /
  SST 移流がこれを読むため、δṁ を res_ro_temp に含めれば追加作業なし
- **advGauge との整合**: δp^HF は定数 pRef・線形ゲージに不変 (定数は Δ と勾配外挿で相殺)。
  `massflux = res + CinfS` の物理化も現行コードのまま正しい
- 運動量の圧力力 p·n·S には触らない

### 4.3 係数の役割分担 (同定可能性)

- **ε_cb**: 物理カットオフ Mach。壁の M=0 ではなく「問題領域の動的に意味のある速度」で決める。
  丘頂 q≈480 Pa → U≈29 m/s → M≈0.085 なので **ε_cb=0.10 に固定** (掃引しない)。
  既存 `precondEps=0.15` とは独立キー (既存 precond 散逸の較正を動かさないため)。
- **C_cb**: 減衰強度。**0 (off, 既定=ビット不変) / 0.05 / 0.10 / 0.20 を掃引**。
  既存の実効圧力結合 σ_min/ε = 0.05/0.15 ≈ 0.333 に**加算**されるため、追加分と総量を
  区別する: 追加強度 C_cb/ε_cb = 0.5/1.0/2.0、**総強度 = 0.333+C_cb/ε_cb = 0.833/1.333/2.333
  = 現行比 2.5×/4×/7×**。C_cb=1 (総 ~31×) は初手では過剰として除外。
- 既存 σ_min×precond の Δp 項とは**重複加算を許容**する (総結合 = σ_min/0.15 + C_cb/0.10)。
  分離較正できる形 (独立レイヤ) にしておくことが仮説検証・較正根拠の保存になる。

### 4.4 config

- `keepDissCbCoeff` (C_cb, 既定 0.0 = off・ビット不変), `keepDissCbEps` (ε_cb, 既定 0.10)。
  **トップレベルキー** (他の keepDiss* と同じ。space 配下に書くと黙って無視される既知の罠に注意)。
- 有効条件: `keepDissType==2` のみ。fdblend/opBlend/precond の有無とは直交。

### 4.5 理論チェック (実装前)

- 追加項は保存形 (面 flux の反対称配布) → 保存性は構造的に OK。
- **エントロピー散逸性は sign gate で面ごとに保証する** (負寄与の許容はしない)。
  δF_f = −k_f δp^HF **g**_f (**g**_f=[1,u_f,v_f,w_f,H_f]ᵀ, k_f=½C_cb S_f/Ur>0) は rank-1 なので、
  entropy 変数ジャンプ Δ**w** との内積 r_f = Δ**w**ᵀ**g**_f を使えば entropy 寄与は
  Δ**w**ᵀδF_f = −k_f δp^HF r_f と閉じる。よって

  $$\delta p_{f,\mathrm{used}}^{HF} = \begin{cases}\delta p_f^{HF} & (\delta p_f^{HF}\, r_f > 0)\\ 0 & (\text{otherwise})\end{cases}$$

  で**面ごとに散逸的 (entropy 寄与非正) を構造的に保証**できる。純粋な圧力市松では
  δp^HF と r_f の符号が一致するので gate はそのまま通す。Δ**w** は既存の dw0..dw4
  (生ジャンプ版 dwr でなく再構成版) を再利用できるため追加コストは内積 1 本。
- 理論チェックスクリプトの検証項目 (「負寄与を測って許容判断」ではない):
  ① gate なしの符号分布 (δp^HF·r_f < 0 の頻度と状態依存性)
  ② gate 後に負寄与が厳密ゼロであること
  ③ 純圧力市松の減衰振幅が gate で削られないこと
  ④ gate 発火率 (乱流的ランダム状態でどの程度の面が落ちるか)

## 5. 実装ステップ

1. **methods 更新**: `methods/convection/implementation.md` の KEEP 散逸レイヤ節に
   本補正の仕様を追記 (開発フロー step 1)。
2. **理論チェックスクリプト**: `solver_density_cuda/tools/verify_cb_pressure_correction.py`
   (§4.5 の sign gate 検証 ①–④ + 2 セル市松の減衰方向検算)。
3. **カーネル**: `convectiveFlux_keep_d.inc.cuh` — keepDissType==2 分岐末尾で
   δp^HF (再構成 pL/pR を再利用、有効面のみ) → sign gate (r_f=Δwᵀg_f) → δṁ →
   res_ro/roUx/roUy/roUz/roe へ加算。massflux への反映は既存の代入 (res_ro_temp 経由) で自動。
4. **config**: `input/solverConfig.hpp/.cpp` に `keepDissCbCoeff`/`keepDissCbEps` 追加
   (`CbCoeff>0 && keepDissJump<1` はハードエラー)、カーネル引数へ配線 (cell/node 両モード)。
   **フルリビルド必須** (struct 変更の既知の罠)。
5. **掃引**: §6 の A/B 掃引を実行、判定。
6. **採用**: 最小有効 C_cb で本番 (run_0015 として run_0014 最新 res から引き継ぎ)。

## 6. 検証

- **単体 / ビルド**: フルリビルド。C_cb=0 で本番構成とビット一致 (res_*.h5 バイナリ比較 100 step)。
- **掃引テストベッド**: 本番格子 (粗格子は AR_z が異なり鋸歯を再現しない恐れがあるため使わない)。
  `run_0014_prod_ddes_opblend/res_56000.h5` (発達場・鋸歯あり) を IC に、**二段階**で回す:
  - **第 1 段 (スクリーニング)**: 2000 step × 4 run (C_cb = 0 基準 / 0.05 / 0.10 / 0.20、各 ~25 分)。
    減衰の方向性・NaN・pseudo 残差健全性のみ判定 (4 ms・スナップ 2–3 枚では統計判定不能)。
  - **第 2 段 (判定)**: off と最良 1 本を **6000–10000 step へ延長**、診断出力
    (`res_ylo` スライス) を **250–500 step 間隔**にして 10–20 スナップ以上確保。
    下記指標 2・3・準定常性はこの段で判定する (`check_quasisteady.py` に丘頂鋸歯量を
    追加するか、同等の専用 VERDICT スクリプトを run に残す)。
  - run 名は case/39 の連番規約に従う: `run_0015_keep_cb_off` / `run_0016_keep_cb_c005` /
    `run_0017_keep_cb_c010` / `run_0018_keep_cb_c020`、採用本番 = `run_0019_keep_cb_prod`。
    **投入時から README の run 表を同期する**。
- **判定指標** (基準 run との A/B、専用後処理スクリプトを run に残す):
  1. 丘頂帯鋸歯: x/h∈[8.3,8.7] 壁線の node-diff RMS と 5 点移動平均除去後の HP-RMS。
     **目標: 基準の ≤1/5 (diff-RMS ≤30 Pa)** [第 1 段から傾向確認可]
  2. 定在性: HP 成分のスナップ間相関 (Δt=0.5–1 ms × 多対) **< 0.3** [第 2 段]
  3. 物理場コスト: z 平均 **Cp** の L2 差 **≤ 0.02** (無次元 Cp 同士で比較、単位を混ぜない)、
     LES 域の解像速度変動 (≥10 スナップ平均からの rms) の変化 **≤ 数 %** [第 2 段]
  4. 箱モード非再発: 壁 ΔP ≤ 3.5 kPa 維持
  5. 収束健全性: pseudo 残差の 1 物理 step 内低下量が基準と同等、NaN/Inf なし
     (`check_convergence.py` VERDICT 併記)
- **副作用検証** (採用候補 C_cb で各 1 run、**共有対流カーネル変更なので node/cell 両モード**):
  - `case/35` 2 セル市松 L1 試験: node/cell (理論スクリプトと実カーネルの符号・減衰率を接続)
  - TGV `case/09` 64³ KE 減衰 A/B: node/cell (滑らか場コスト、期待 <1%; DNS 参照 `case/09/ref_dns`)
  - free-stream `case/33` 非直交: node/cell (機械ゼロ維持)
- **効かなかった場合**: C_cb=0.2 でも目標未達なら、(a) Δz 半減の局所スラブ試験で
  波長が 2–3Δz のまま追随するか (数値モードの決定的証明)、(b) 第二候補
  (precond 行列の σ_p/σ_u 分離、§2 の事前検証つき) へ進む。

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/convection/convectiveFlux_keep_d.inc.cuh` (cell/node 共用)
- `solver_density_cuda/input/solverConfig.hpp/.cpp`、カーネル呼び出し側の引数配線
- `solver_density_cuda/tools/verify_cb_pressure_correction.py` (新規)
- `methods/convection/implementation.md` + `methods/index.md` 整合確認
- 既定 off のため既存ケースへの影響なし (C_cb=0 ビット一致で担保)

## 8. 完了条件

- [ ] `methods/convection/implementation.md` 更新済み
- [ ] 理論チェック (§4.5) 実施・結果記録
- [ ] C_cb=0 ビット一致確認
- [ ] 掃引完了・§6 判定クリア・最小有効 C_cb 決定
- [ ] 本番 run_0015 反映 (鋸歯消失を本番格子で確認)
- [ ] status: done → `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-07-22` — 初稿。run_0014 の丘頂 2–3Δz 定在鋸歯診断と、σ_min 独立の圧力結合強化
  (独立高周波補正レイヤ vs 行列列別スケール) の比較検討を経て、独立レイヤ案を第一候補に採用。
- `2026-07-23` — **⚠ 撤回 (ユーザー指摘による)**: 本計画の動機だった「丘頂 spanwise 2-3Δz 定在鋸歯」は
  **抽出アーチファクトだった**。xcluster により丘頂の Δx/h≈0.031-0.033 に対し、スパン線抽出の
  バンド幅 |x/h−8.5|<0.03 (幅 0.06) が **隣接 2 列 (x/h=8.483, 8.515) を混在**させ、丘頂の急な
  x 方向圧力勾配 (~100 Pa/列) が z ソートで交互に現れてジグザグを偽装 (「定在」も定常な平均勾配の
  アーチファクト)。**単一列で再解析すると丘頂は再付着域 (物理対照) と同質の滑らかな解像乱流**
  (lag-1 corr +0.58〜+0.82, 統計量同水準)。SLAU 比 200 倍も無意味な比較 (定常 RANS は z 一様)。
  訂正図: case/39 `sawtooth_retraction.png`。**教訓: スパン線抽出は必ず列数を確認する** (バンド
  選択ではなく x 値の厳密一致で列を取る)。
  cb 補正自体の実装・検証 (L1 純 2Δ に有効・free-stream 安全・entropy gate) は有効なまま opt-in で
  残置。keepDissOpBlendRaw も散逸ゾーニングの忠実化として理屈は残る (本番採否は保留)。
  case/39 に既知の数値アーチファクトは**箱モード (opblend で解決済み) のみ**となった。
- `2026-07-23` — (撤回前の記録) 本番格子掃引: 不採用と判断していた。run_0015 (off) / run_0016
  (C_cb=0.05) / run_0017-0018 (0.1/0.2): 0.05 は鋸歯 diff-RMS が off と同一 (124.6 vs 124.6 Pa
  @2000step)、0.1/0.2 は implicit でも step≤100 で NaN。**実測原因** (res_76000 丘頂 z 面 6300 面):
  entropy gate は 100% 通過 (無罪)、**δp^HF 実効振幅が生 Δp の中央値 6%** — 鋸歯は純 2Δ でなく
  2-3Δz で中心勾配に「見える」ため、再構成ベース抽出の保護帯に入ってしまう。本補正は純 2Δ 市松
  専用 (L1 では基準の 10 倍速) として機能は正しいが、本ターゲットには構造的に届かない。
  **派生知見**: jump2 フィルタは全散逸経路 (opblend 標準 Roe 側含む) に掛かるため、fdblend で
  σ→1 にしても実効ジャンプ 6% で散逸ゾーニングが無効化されていた。これを直す
  `keepDissOpBlendRaw` (標準 Roe 側=生ジャンプ, DES 定石) を実装・検証したが、**鋸歯は不変**
  (run_0019: 143.3 vs 124.6 Pa, 自然変動幅内) — 鋸歯の主体は fd~0.99 の LES 帯 (σ_f=σ_min=純
  precond, 生ジャンプ増分の適用外) にあり、覆われる grey 帯だけ強めても再インプリントされる。
  **残る選択肢**: (i) 現状容認 (有界・局所・0.1-0.3 q_crest, 丘頂近傍 Cp に注意書き)
  (ii) メッシュ: nz 60→90-120 (コスト比例増) or 丘頂 y1 緩和で AR_z 低減
  (iii) 異方性選択的散逸: セル最小寸法 h_min の事前計算配列を追加し、面法線方向スケール
  h_n=V/S_f との比で z 面音響散逸を局所増強 (新設計・要計画更新)。
- `2026-07-23` — 実装完了 (カーネル+config+理論スクリプト全 PASS+ビット一致=run-to-run ノイズ範囲)。
  実装時の重要知見 2 件: ① 絶対圧力での δp^HF は float32 非直交 free-stream を破壊 → 差分形+
  sign-property minmod で根治 (§4.1)。② explicit RK3 CFL0.5 では実効係数 C_cb/ε≦~0.2 が安定限界
  (C/ε=1.0 で市松モード自体が overshoot 発散 ×3.7/step, case/35)。implicit 本番は掃引で確認。
- `2026-07-23` — 外部レビュー (Codex) の 6 指摘を反映: ① Ur を M=1 で上限 (`lowMachUr` 再利用)
  ② 生 Δp フォールバック廃止 (無効面は δp^HF=0、`CbCoeff>0 && keepDissJump<1` はハードエラー)
  ③ 負寄与の許容判断を廃し **sign gate (δp^HF·r_f>0 のみ有効) で面ごとに entropy 散逸性を保証**
  ④ 係数表を追加分/総量 (現行比 2.5×/4×/7×) に訂正 ⑤ 掃引を二段階化 (2000 step スクリーニング +
  off/最良の 6000–10000 step 延長・出力 250–500 step 間隔)、Cp 指標の単位統一
  ⑥ run 命名を連番規約 (`run_0015_keep_cb_off` 等) に修正、副作用検証に case/35・node/cell 両モードを追加。

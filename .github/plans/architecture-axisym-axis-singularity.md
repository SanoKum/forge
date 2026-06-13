# 軸対称 近軸の数値問題 (軸中心 k スパイク) の根本原因特定

## メタ

- **area**: `architecture`
- **status**: `in_progress`  <!-- 根本原因を特定。安価な根治は未発見 (global double / explicit が実証済の対処) -->
- **related_docs**:
  - `docs/turbulence/theory.md` (§7.5)
  - `docs/time_integration/` (block DPLUR)
- **related_plans**:
  - `architecture-axisymmetric.md` (B 流儀 r 重み実装)
  - `time_integration-implicit-stable-cfl.md` (block DPLUR)
  - `turbulence-kato-launder.md` (← 本件の対処としては**誤り**。下記参照)
- **related_docs(手順)**: `.github/forge-su2-cross-check.md`
- **created**: `2026-06-12`
- **updated**: `2026-06-13`
- **owner**: `CFD Dev`

## 1. 結論 (根本原因 = float32 の陰解法が近軸を収束させきれない)

case 29 軸対称 RANS の「軸中心 k スパイク」の根本原因は、**float32 の陰解法 (block-DPLUR) が
近軸第一セルの半径方向速度 `u_r` を収束させきれない**こと。乱流モデル (フープ項 / Kato–Launder)
は無関係で、それらは下流の**対症療法**にすぎない。

### 確定した因果連鎖

1. **平均流の欠陥**: float32 **陰解法**では、軸隣接の第一セルの `u_r` (forge の `Uy`) が物理値
   (±15〜30 m/s) ではなく **≈0 に張り付く**。隣の第二セルは正しい値に跳ぶため、第一⇄第二セル間に
   巨大な偽 `∂u_r/∂r` が生じる。
2. **乱流への波及**: その偽ひずみが SST 生産 `P_k=μ_t S²` を駆動し、軸中心 k がスパイク
   (k_axis 12〜16 vs k_core 3〜4)。フープ項 `u_r/r` を生産から外しても **6% しか減らない**=フープは主犯でない。
3. **SU2 との対比** (同一メッシュ・同一 BC): SU2 (頂点中心・**軸上に節点**・倍精度) は `u_r` が 0 から
   滑らかに立ち上がり k は軸で**最小**(スパイク無し)。forge (セル中心・第一セル中心が r=Δr/2) で起きる。

### 決定的証拠 (precision を切り分け)

x=40mm 第一セル `Uy` の収束履歴 (laminar conical):

| 構成 | 5k step | 20k step | 判定 |
| --- | --- | --- | --- |
| float32 **explicit** (RK) | — | **−14.9** ✓ | 正しい |
| **global double** (全倍精度) implicit | **−15.1** | −15.1 ✓ | 即収束・正しい |
| float32 **implicit** (block-DPLUR) | −0.2 | **−0.6** | 約1000倍遅く這う=固着 |

- **収束解そのものは正しい** (explicit float = double)。問題は float32 **陰解法**が近軸第一セルを
  ~1000倍遅くしか収束させないこと (1run では動ききらず固着して見える)。
- 残差で言うと、`rms_ro`≈3e-5 でも `rms_roUy`≈1e-2 停滞・`rms_roK` 増大 = 近軸が未収束。
  → **収束は `rms_ro` だけで判断してはならない** ([AGENTS.md] 収束ルールへ反映済)。

## 2. 試して失敗した対処 (すべて評価済・効果なし)

「桁落ちは陰解法のどこか1箇所」という仮説で部分倍精度を試したが、**いずれも固着のまま**:

| 倍精度にした箇所 | 結果 |
| --- | --- |
| 5×5 ブロック solve (`solve_5x5`) | 効果なし |
| 軸対称ソースヤコビアン除去 | 効果なし |
| FVS 面ヤコビアン `build_jacobian_split` (`a_plus`/`k_off`) | 効果なし (float の固有積は元々正確) |
| **DPLUR 反復まるごと倍精度** (diag/rhs/neighbor/solve、入出力のみ float) | **効果なし** |
| `cfl_pseudo` 引き上げ (2→4→8) | 効果なし (over-damping/CFL ではない) |
| `cfl_pseudo` 引き下げ (2→0.5→0.3) | 効果なし (−0.64→−0.58→−0.53、固着のまま) |
| scalar DPLUR (`blockDPLUR=0`) | **全発散** (cfl2 で step~9k に全場 NaN)。block DPLUR が安定性に必須 |
| `nStepInner` 増 (20→100→200 sweep) | **全く不変** (−0.64 でビット一致)。線形系は既に収束済=**Krylov(FGMRES)でも直らない** |
| `doubleResidual` (残差を double バッファに蓄積) | **効果なし** (−0.64)。桁落ちは atomic sum ではなく**float 状態由来の per-face 値**にある |
| **global double** (状態も double) | **効く** |

**SU2 との対比 (su2.log で確認)**: SU2 = 後退Euler + **FGMRES + ILU(0)** + **double**。forge = 後退Euler defect-correction
+ **block-DPLUR (LU-SGS sweep 固定回数)** + **float**。`nStepInner` テストで forge の LU-SGS は線形系を既に
解ききっている (sweep 増で不変) と判明したので、**SU2 の優位は FGMRES という手法ではなく double 精度**。
`doubleResidual` テストで残差蓄積を倍精度化しても直らず、**効くのは状態 (P,ρ 等) の倍精度**のみ
= float 状態から計算した per-face 流束が近軸で桁落ちし、陰解の `D⁻¹` がそれを増幅する。陽解は `R·dt/V` で許容。
→ **部分倍精度・手法変更では根治せず、状態の倍精度 (global double) が必要**。SU2 の `USE_MIXED_PRECISION`
(本体 double・前処理のみ float) も同じ結論。

**決定的比較 (同一 cfl=0.5・同一 20000 step)**: 陽解法 (`timeIntegration=3`) は第一セル Uy=−14.9 (正)、
陰 block-DPLUR (`timeIntegration=11`) は −0.58 (固着)。cfl も step 数も同一なので、差は
「**陽の点更新 vs 陰の block-DPLUR defect-correction そのもの**」。乱流なし (laminar) でも生じる
= **平均速度場の問題**であり、block-DPLUR (float) が近軸第一セルを構造的に収束させられない。

- 切り分けの帰結: double-D (DPLUR反復) も double-R蓄積 (doubleResidual) も単独では**効かず**、
  **状態 (P,ρ等) を double にする global double のみ効く**。近軸の悪条件な陰解 Jacobian `D⁻¹` が、
  **float 状態から計算した per-face 流束の桁落ち**を増幅するため。**陽解法は `R·dt/V` の
  スケールでこれを許容する** (だから float でも正しい)。→ 根治には**状態の倍精度**が要る。
- 乱流モデル側の対処はすべて**無関係 (マスク)**: フープ除去 (6% のみ)、`dilatationCorrection`
  (トレース除去では `u_r/r` は偏差成分に残り消えない)、Kato–Launder (軸上 `Ω→0` でたまたま生産が
  落ちるだけ)。`turbulence-kato-launder.md` の「軸スパイクの対処」という位置づけは**誤り**。

## 3. 実用上の指針 (確定)

1. **推力・平均流は妥当**: case 29 の推力 (λ=0.984, mdot チョーク整合) は検証済。スパイクは
   **近軸限定の乱流過大評価**で、大域の推力・平均流の妥当性は損なわない。
2. **近軸の k を信用しない / 過細分を避ける**: 軸対称 implicit では軸第一セルの乱流量は精度限界。
   軸セルを極端に細分しない (悪化する)。
3. **近軸の乱流精度が要る場合**: **global double**(実証済・正しい。ただしメモリ2倍・低速、
   かつ double ビルドでの RANS 安定性は要確認) か、**explicit**(平均流は float で正しい。
   ただし SST は剛体ソースで cfl 要調整・warm-start 推奨)。

## 4. 根治の候補 (要検討)

部分倍精度 (DPLUR反復 / 残差蓄積) はすべて**検証して無効** (§2)。安価な根治は存在しないと判明。残る方向:

- **状態の倍精度 (global double / 状態だけ double)**: 唯一実証済。近軸が要するのは float 状態から
  計算する per-face 流束の精度なので、最低でも状態 (P,ρ,roU) を double 化する必要がある。
  flux/勾配も double にしないと per-face 値が float 精度のままなので、実質 global double。
  消費者GPU では FP64 が遅い (RTX3060 で ~1/32) のが難点。double ビルドでの RANS 安定性は要確認
  (本調査では double SST が発散した)。
- **近軸セルだけ陽的更新するハイブリッド**: ユーザ却下 (ハック的で不可)。
- **軸セルの定式化変更**: 第1セル中心を軸上に置く (SU2 流の頂点中心) / axis-cell の特別扱い。大規模・要設計。

## 5. 検証 (根治実装時)

- laminar conical の第一セル `Uy` が float implicit で −15 近傍に即収束 (現状 −0.6 固着)。
- SST conical の軸中心 k が SU2 同様に「軸で最小」(現状スパイク)。
- 全軸対称回帰 (case 26 flat_plate / case 27 CEA / case 29 推力 mdot/λ) で悪化なし・速度低下が許容範囲。

## 6. 変更ログ

- `2026-06-12` — (旧説。誤りを含む) strain アノマリー / planar 面積 / グリッド発散 / 近軸チェッカーボード説。
- `2026-06-13` — **根本原因を確定**。SU2 クロスチェック (同一メッシュ・BC) と precision 切り分けで:
  ①フープ項は主犯でない (除去で 6% のみ)。②flux スキーム非依存 (SLAU=ROE 両方で発生)。
  ③軸ゴースト・勾配は正常 (軸面流束は r 重みで実質ゼロ、勾配は planar で正しい)。
  ④**float32 陰解法 (block-DPLUR) が近軸第一セルの `u_r` を収束させきれない**のが真因。explicit float
  と global double はいずれも正しい。Kato–Launder は無関係なマスクと確定
  (`turbulence-kato-launder.md` の位置づけを訂正)。
- `2026-06-13` (追試) — 安価な根治を網羅的に検証し**すべて無効**と確定:
  部分倍精度 (5×5 solve / source-Jac / FVS-Jac / DPLUR反復全体)、`cfl_pseudo` 上下、`nStepInner` 増
  (線形系は収束済→FGMRES でも不可)、scalar DPLUR (発散)、**`doubleResidual` (残差の倍精度蓄積)**。
  **global double (状態の倍精度) のみ有効**。SU2 (FGMRES+ILU+double) との対比から、**SU2 の優位は手法
  ではなく double 精度**であり、桁落ちは残差蓄積ではなく**float 状態から計算する per-face 流束**に由来する
  と判明。→ 根治は状態の倍精度 (≈global double) が必要、近軸ハイブリッド陽的はユーザ却下 (§4)。

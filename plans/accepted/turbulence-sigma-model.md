# σ-model SGS (Nicoud 2011) の追加

## メタ

- **area**: `turbulence`
- **status**: `done`
- **related_docs**:
  - `methods/turbulence/implementation.md`
- **related_plans**:
  - `plans/accepted/turbulence-wale-fix.md` (WALE 修正で判明した静的モデルの遷移前倒し問題が動機)
- **created**: `2026-07-19`
- **owner**: CFD Dev

## 1. 目的

修正後の WALE は遷移流 (64³ TGV) で ν_t を早く立てすぎ、ILES より DNS 一致が悪い
(ピーク t\*=7.84 vs DNS 8.98)。静的モデルの中で誤検知が最も少ない **σ-model**
(Nicoud, Toda, Cabrit, Bose, Lee, PoF 2011) を `LESmodel: 2` として追加し、
未解像高 Re 乱流用の明示 SGS の第一候補を整備する。

## 2. スコープ

- **やる**: `SIGMA_d` カーネル (LESorRANS=1, LESmodel=2)。$\nu_t=(C_\sigma\Delta)^2 D_\sigma$,
  $D_\sigma=\sigma_3(\sigma_1-\sigma_2)(\sigma_2-\sigma_3)/\sigma_1^2$, $C_\sigma=1.35$,
  $\Delta=V^{1/3}$。特異値は $G=g^Tg$ の固有値の平方根 (3×3 SPD の三角関数閉形式,
  カーネル内 double 評価)。wall_dist 不要 (壁減衰は $D_\sigma$ の性質で自動)。
- **やらない**: dynamic 化 / C_σ 再較正 (論文推奨 1.35 のまま) / WALE の置換 (併存)。

## 3. 設計方針

- **モデルの設計性質 (python で事前検証)**: 2 成分流 (σ3=0)・純せん断 (σ2=σ3=0)・
  剛体回転 (σ1=σ2, σ3... D_σ=0)・軸対称伸長の一部で ν_t=0。**TG 初期場は w=0 の
  2 成分流なので ν_t≡0 = 遷移前倒しへの構造的耐性** (WALE は初期から ν_t~0.26μ)。
- 固有値閉形式 (SPD): I1=tr G, I2=(tr²G−tr G²)/2, I3=det G →三角関数式。float32 の
  丸めに対し α1<0・acos 引数の [−1,1] クランプ・λ の非負クランプで正則化。
- 検証: TGV 64³ node で (a) σ-model 単独 (b) +ES σ0.02 jump2。期待: 層流期 ν_t≈0・
  ピーク時刻が WALE (7.84) より DNS (8.98) に寄る。

## 変更ログ

- 2026-07-19: 起票。
- 2026-07-19: **実装・検証完了**。
  - 実装: `SIGMA_d` (turbulent_viscosity_d.cu, LESorRANS=1 & LESmodel=2)。固有値はカーネル内
    double の三角関数閉形式+ソート+クランプ。`tools/verify_sigma_model.py` で SVD 一致
    (rel 1.5e-10)・設計ゼロ性質 (純せん断/剛体回転/等方・軸対称伸長/TG 2成分場) ALL PASS。
  - 検証 (node 64³ TGV Re=1600, run_0044/0045): 層流期 (t*=0.28) の ν_t は WALE 比 25% 小さい
    (設計どおり) が、圧縮性 TG は早期に 3 次元化するためゼロ性質の効く時間が短く、発達後は
    C_σ=1.35 の散逸が勝って **K/K0(10) −5.5% (WALE −3.4% より悪い)・ピーク t*=7.84 (前倒し同等)**。
  - **結論**: 64³ TGV 級の解像/遷移流では静的 SGS はどちらも逆効果で、**ILES + matrix ES
    (σ=0.02, keepDissJump=2) が最良 (−1.4%) のまま**。σ-model は `LESmodel: 2` として整備済み、
    本領 (壁乱流・回転流・未解像高 Re) での WALE との比較は将来のケースで行う。

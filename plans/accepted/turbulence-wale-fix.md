# WALE の不活性バグ (wall_dist≡0) と Sd テンソル式の修正

## メタ

- **area**: `turbulence`
- **status**: `done`
- **related_docs**:
  - `methods/turbulence/implementation.md` (WALE 節)
- **related_plans**:
  - `plans/accepted/convection-keep-diss-recon-jump.md` (L3 検証が本バグの影響を受ける)
- **created**: `2026-07-19`
- **owner**: CFD Dev

## 1. 目的

ユーザー指摘「WALE にしているのに vis_turb が 0 では?」の調査で 2 つのバグを確認した。

1. **不活性バグ**: `WALE_d` の長さスケールが `Ls = min(κ·wall_dist, Cw·V^{1/3})`。壁の無い
   メッシュ (周期 box 等) では converter が `wall_dist≡0` のまま → Ls=0 → **vis_turb≡0**。
   これまでの TGV L3「WALE 併用」run (32³/64³, cell/node すべて) は実際には **SGS なしの
   ILES (KEEP+分子粘性+ES散逸)** だった。
2. **Sd テンソル式の誤り**: WALE (Nicoud & Ducros 1999) の
   $S^d_{ij}=\tfrac12(\bar g^2_{ij}+\bar g^2_{ji})-\tfrac13\delta_{ij}\bar g^2_{kk}$
   ($\bar g^2_{ij}=g_{ik}g_{kj}$ = **行列 2 乗**) に対し、実装は成分ごとの 2 乗
   ($\tfrac12(g_{12}^2+g_{21}^2)$ 等) になっており本来の WALE 演算子ではない。

## 2. スコープ

- **やる**: ① 壁なしメッシュの wall_dist ガード (solver 読込後、`max(wall_dist)<=0` なら
  1e30 で充填・ログ出力。メッシュ再変換不要・壁ありメッシュはビット不変)。
  ② `WALE_d` の Sd を Nicoud-Ducros の行列 2 乗形へ修正 (コード前に python で数値検証)。
  ③ 本物の WALE で 64³ TGV node を再検証し、ILES との差・σ 推奨を再確認。
  ④ 既存文書の「WALE 併用」表記を「ILES (WALE 不活性)」へ訂正。
- **やらない**: WALE 定数 Cw の再較正 (標準値 0.325 のまま) / 壁近傍 κd クリップの是非
  (壁ありケースの挙動変更は別途)。

## 3. 設計方針

- wall_dist ガードは **variables.cpp の VALUE 読込直後** (host)。判定は「全 CV で ≤0」
  = 壁情報なしの converter 出力。壁ありメッシュは壁面 CV (node の壁ノード d=0) があっても
  max>0 なので触らない。
- Sd 修正は $\bar g^2 = g\cdot g$ を陽に組み、対称化+トレース減算。検証は乱数勾配行列で
  numpy の行列積と突き合わせ (verify_wale_sd.py)。
- 影響範囲: WALE が実際に動いていた既存ケースは「壁あり×LESorRANS=1」のみ (該当 run は
  要棚卸し)。TGV 系は不活性だったため②の修正で過去結果は変わらない (0 のまま→修正後有効)。

## 変更ログ

- 2026-07-19: 起票。真因確認 (wall_dist≡0 実測・Sd 式差異のコード確認)。
- 2026-07-19: **実装・検証完了**。
  - 実装: ① `variables.cpp` 読込後ガード (max(wall_dist)≤0 → 1e30 充填+ログ; 壁ありメッシュ不変)。
    ② `WALE_d` の Sd を行列 2 乗形へ (python 検証 PASS: numpy 行列積と一致・純せん断で Sd=0)。
  - 検証 (node 64³ TGV Re=1600): vis_turb が活性化 (mean ~0.26×分子粘性 @層流期)。
    **本物 WALE は遷移先回り散逸でピーク t*=7.84 (DNS 8.98)・終値 −3.4% と ILES より悪化**
    (run_0042/0043)。→ **解像/遷移流の推奨は WALE off (LESorRANS:0) + matrix ES σ=0.02 jump=2**
    (−1.4%・ピーク 8.96)。WALE は未解像高 Re 乱流用オプションへ格下げ (適用検証は今後)。
  - 過去の「WALE 併用」表記 (case/09 32³/64³ L3・関連 plan・memory) を「ILES (WALE 不活性)」へ訂正。
    数値結果自体は ILES として全て有効。

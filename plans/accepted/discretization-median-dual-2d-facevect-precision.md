# 2D median-dual 面ベクトル・境界半割面の桁落ち除去 (ローカル原点 + double)

## メタ

- **area**: `discretization / architecture`
- **status**: `done`
- **related_docs**:
  - `methods/discretization.md` §2.5.6 (双対幾何の数値精度 — 本 plan で 2D 適用を追記)
- **related_plans**:
  - [`architecture-median-dual-3d-double-geometry.md`](../accepted/architecture-median-dual-3d-double-geometry.md) (同型修正の 3D 版・レシピの出典)
  - `boundary-node-nozzle-wall-outlet-stability.md` §2.9 (2D shoelace [体積・重心] の先行修正)
- **created**: `2026-08-31`
- **owner**: `sano + agent`

## 1. 目的

2D `buildMedianDual` (`mesh/gmshReader.hpp`) の**双対面ベクトルと境界半割面**が絶対座標
float32 のまま残っており、薄壁 y+~1 メッシュで median-dual 変換が閉性破綻する。
実害 (case/42 M6 NS, 2026-08-30): 第一セル ~2.4 μm (wf 4.5e-5 × r_t 53.75 mm, ni1700,
全長 6.3 m) で `dual faces not closed` (normalized 0.099, 壁 bcond 面積 3 % 欠損) となり
変換拒否 → y+~1.4 (第一セル ~3.5 μm) への退避を強いられた。§2.5.6 の規約
(ローカル原点 + double) を 2D 側にも適用して除去する。

## 2. 現状分析 (2026-08-31 コード検分)

- **修正済み**: 2D sub-CV shoelace (体積・重心, A 相対 + double)、3D 全幾何 (M/N 相対 Newell +
  double 蓄積)。
- **露出**: ① 2D 双対面ベクトル — `Mx=0.5(A+B)`, `sx=Gx−Mx`, rotate, **向き整合内積
  `nx·eABx+ny·eABy`** が全て geom_float (float32)・絶対座標。x~6.3 m で丸め ~0.5 μm。
  スリバー面では**向き内積が誤反転**し面 1 枚 (mm 級) が逆符号で閉性に入る — 観測された
  normalized 0.099 (=mm 級欠損) は丸め積算では説明できず、この誤反転が主犯。
  ② 2D 境界半割面 — `planes[].surfVect` (makeMesh の float 絶対座標差) の半分を
  geom_float で蓄積 (`bnodeAccum`/`halfByOwner`/`hcentByOwner`)。
- 閉性は「同一座標集合から厳密演算で作れば恒等的にゼロ」の組合せ恒等式なので、
  演算を double にすれば float32 座標のままでも閉性は面ベクトル格納丸め (~1e-7 相対) まで落ちる。
- **追加検分 (①② 修正後も ni1700 が 0.099 のまま再現) → 真犯人は makeMesh の 2D CW 判定**:
  `isCW` の shoelace が float **絶対座標の積和** (`x0*y1−x1*y0`, 積 ~3 の ulp ~2.4e-7) で、
  スリバーセルの符号付き面積 (~1e-9 m²) を 8 桁上回るノイズに埋まり判定がほぼランダム化。
  その結果**境界 surfVect が誤反転** (壁 ~50 本 = 面積 3 % 欠損・閉性 mm 級)。①②の双対側修正は
  「整向済み surfVect」を信頼する設計なので上流の誤反転が素通りしていた。

## 3. スコープ

- **やる**: ① 2D `buildMedianDual` の面ベクトル計算のエッジ中点 M 相対 + double 化
  (向き内積・面重心含む)、② 境界半割面の節点座標からの double 再構成 (向きは整向済み
  `surfVect` の符号に一致させる) + 蓄積 double 化、③ **makeMesh の 2D CW 判定 (isCW shoelace) の
  先頭ノード相対 + double 化** (スコープ拡張 2026-08-31: 真犯人がここと判明したため)。
- **やらない**: `geom_float` 自体の double 化 / 3D 側の makeMesh 整向内積 (`Db`/`D`) の
  double 化 (同種リスクは残るが 3D 薄壁の検証ケースが無いため別途 — 残件に記載) /
  ソルバ側の変更。

## 4. 設計方針

- 面ベクトル: `M`, `e_AB`, `G−M` を double (float 座標の promote は厳密) で計算し、
  rotate・向き判定・集約を double、最後に `geom_float` へ cast して格納。
  平行移動不変なので厳密演算では従来と同値 (健全メッシュはビット近傍で不変)。
- 境界半割面: `surfVect` を読む代わりにエッジ (A,B) 座標から double で法線を再構成し、
  `surfVect` との内積で外向き符号だけ合わせる (makeMesh の整向規約を保持)。
- 検証ゲート (`closTol 1e-3`) は据え置き — 修正後は正常メッシュで ~1e-7 台に落ちるため
  余裕が 4 桁できる。

## 5. 検証

1. **converter 診断 A/B**: M5 NS (800×97, 3.6 μm)・M6 NS (1250×97, 3.5 μm) — 修正前後で
   closure が桁で改善し、体積総和 relErr 不変であること。
2. **これまで不可能だったメッシュ**: M6 ni1700/wf4.5e-5 (2.4 μm) が閉性ゲート内で変換成功。
3. **CFD 回帰**: M5 Euler センター点 (run_0072 相当) を新変換メッシュで再実行し
   dM/ε_M が同水準 (run_0106)。
4. **payoff**: M6 y+~1 (ni1700/wf4.5e-5) の NS v1 を実走し (run_0107)、y+1.4 版 (run_0104)
   と出口コア M が整合すること。

## 変更ログ

- 2026-08-31: 起票・コード検分 (誤反転主犯説)・実装着手。
- 2026-08-31: 実装・検証完了 (①②③, `mesh/gmshReader.hpp` 内で完結)。検証結果:
  - converter 診断: M5 NS 800×97 / M6 NS 1250×97 とも閉性 normalized ~**1e-7** (4 桁改善)、
    壁 bcond 面積が primal と 1e-5 一致 (修正前 M6-1700 は 3 % 欠損)。体積総和 relErr 不変 (~e-10)。
  - **ni1700/wf4.5e-5 (第一セル 2.4 μm) が変換成功**: CW 誤検出 0 (修正前は ~50 本の
    境界 surfVect 誤反転)、閉性 9.5e-8。
  - CFD 回帰 (run_0106 = M5 Euler センター点を新変換で再実行): dM 0.000430 vs run_0096 0.000429、
    M_axis_exit 差 7e-6、ε_M 同一 — 健全メッシュの挙動は不変。
  - payoff (run_0107 = M6 y+~1 真値 NS v1): 出口コア M 5.9640 (−0.60 %) — y+1.4 版 run_0104 の
    5.9660 (−0.57 %) と 0.03 % 整合。**y+1.4 退避は不要になった**。
  - 残件: 3D makeMesh 整向内積 (`Db`/`D`) の double 化 (同種リスク、3D 薄壁ケース登場時に)。
- 2026-08-31 (残件消化): **3D makeMesh 整向内積 (`Db`/`D`) も double 化** — 面先頭ノードを
  ローカル原点に、面重心・セル重心 (ノード算術平均 = centCoords と同定義) を double で再構成して
  評価。A/B (case/35 tg.msh, node 3D): 整向 fix 数・閉性 (2.47e-6)・体積総和が不変、
  **出力 h5 の全 dataset ビット一致** (健全メッシュで挙動変化ゼロ)。3D 薄壁の正例検証は
  該当メッシュ登場時 (2D で同レシピの根治実績あり)。残件は解消。

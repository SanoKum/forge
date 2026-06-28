# cell 周期境界の対流保存バグ (seam 非周期化)

## メタ

- **area**: `boundary`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md)
  - 検証ケース: [`case/09.Taylor-Green`](../../case/09.Taylor-Green/README.md)
- **created**: `2026-06-28`
- **owner**: `CFD Dev`

## 1. 目的

cell (collocated primal) モードの triply-periodic 境界が**運動量を保存しない**バグを修正する。Taylor-Green (全周期箱) で全運動量 $\sum\rho\mathbf u V$ (理論=0) が線形に注入され、周期 seam に不連続が成長して「流れが周期的にならない」。node (median-dual) は同条件で厳密保存。

## 2. 切り分け済みの事実 (2026-06-28)

`case/09.Taylor-Green` の cell pure-KEEP 非粘性 (`run_0014`) + 一時診断 (main.cpp に `CONSV_DIAG` で Σres / ghost seam-flux を計測) で確認:

- **t=0 (解析 TGV, 厳密周期) では保存**: 対流流束直後の interior $\sum\mathrm{res}=O(10^{-8})$、周期ペア (bcond1+bcond2 等) の ghost seam-flux が**完全相殺**。
- **場が発達すると非保存**: RK ステージ2以降 interior $\sum\mathrm{res}_{\rho u}=+1.55\mathrm{e}{-2}$, $\rho v=-1.55\mathrm{e}{-2}$ (x↔y 反対称, TGV の対称性を保つ)。ghost seam-flux のペア相殺が崩れる (y-pair で normal 方向に +0.015 残る)。$\sum_{\mathrm{ALL}}\mathrm{res}=O(10^{-8})$ は保たれる (局所 atomicAdd は保存) ので、**「ghost に入れた運動量が partner 内部セルへ渡らず捨てられている」**。
- **場の非周期化を定量化**: y-seam を跨ぐ $\langle|\Delta u|\rangle$ は step0=1.6e-9 (連続) → step100=0.032 (内部の1.45倍) → step350=0.30 (内部の5.46倍)。seam に不連続が成長。
- **スキーム非依存**: pure KEEP (勾配・重心不使用) でも ROE (MUSCL+Roe) でも発生。KEEP は低散逸ゆえ顕在化が早い。
- **一様流は全方向で厳密保存** (`run_fs_x/y/z`, 整合した一様状態): 法線の反対称 $n_A=-n_B$・面積一致・ghost=partner コピーは健全。
- **partner 対応・並進ベクトルは幾何的に正しい** (python で setPeriodicPartner を再現し横方向誤差 0、ミスマッチ 0)。
- **node は厳密保存** (多重度補正後: 質量厳密・運動量 ~1e-7・seam 連続)。node は周期面を移流ループから除外し `periodicNodeGather` で残差合算する方式。

## 3. 矛盾点と次の計測

理屈上は「ghost=partner かつ KEEP/ROE 流束が $f(B,A,-n)=-f(A,B,n)$ で反対称」なら seam はステージに依らず厳密相殺し保存するはず。にもかかわらず場の発達とともに相殺が崩れる。残る容疑:

- RK4 多ステージでの ghost 更新の不整合 (中間状態で ghost プリミティブが partner と一致していない可能性)。`periodic_d` は毎ステージ呼ばれるが、`dependentVariables`→`applyBconds(periodic_d)`→`convectiveFlux` の順序と、ghost 保存量 (RK 中間 Q) の扱いを精査する。
- **次タスク**: 正しいローカル添字で flux 時点の `Ux/Uy/P[ghost]` vs `[partnerCell]` を直接ダンプ (前回診断は `iCells_ghst` と `partnerCellID` のローカル順序食い違いで 0.43 の偽差を出した。`bplane` 配列の構築順を揃えて再計測)。一致していれば流束式/幾何、不一致なら ghost 更新タイミングが原因と確定する。

## 4. 関連

- 別バグ (本件と独立): [`boundaryCond_d.cu`](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu) `periodic_d_wrapper` が `dtheta` を `bc.inputInts["dtheta"]` (int マップ・存在しないキー→0) から読む。回転周期 (type=1) が常に dtheta=0 で無効化される。Cartesian (type=0) では 0 が正しく無害。`bc.inputFloats["dtheta"]` から読むよう修正が必要。
- 暫定運用: 全周期の保存性が要るケース (TGV 等) は **node モード**を使う。

## 変更ログ

- `2026-06-28` — 調査着手。上記の切り分けを実施し、原因を cell 周期 seam の対流保存に局在化。修正は次タスク。

# 3D median-dual 幾何の桁落ち除去 (buildMedianDual3D のローカル原点化 + double 蓄積)

## メタ

- **area**: `architecture / discretization`
- **status**: `draft`
- **related_docs**:
  - `methods/discretization.md` (M4 3D median-dual の節に精度注意を追記予定)
- **related_plans**:
  - [`discretization-median-dual-3d.md`](discretization-median-dual-3d.md) (親: M4 3D 双対生成)
  - [`boundary-node-nozzle-wall-outlet-stability.md`](boundary-node-nozzle-wall-outlet-stability.md) §2.9 (2D 版の同種修正・実害の記録)
- **created**: `2026-08-11`
- **owner**: `sano`

## 1. 目的

2D `buildMedianDual` で実害が出た float32 幾何桁落ち (y+≈1 スリバー CV で双対重心が
~100 μm 狂い CV 外へ → node `wall_dist` 破壊 → SST 爆発; §2.9 で A 相対座標+double 化により
修正済み) と**同型のリスクが 3D 版 `buildMedianDual3D` に残っている**。3D で薄壁 (y+~1 級)
node メッシュを使う前に除去する。3D の RANS/DES 壁解像は WMLES・DDES 計画の前提であり、
未修正のまま 3D 細壁に進むと 2D と同じ「原因不明の初手発散」を踏む。

## 2. 現状分析 (2026-08-11 コード検分)

**既に安全な部分** (2D と違い 3D は最初から局所差分 double で書かれている):

- `tetVol` (双対体積): 4 点を先頭点相対の差分にしてから外積 — 桁落ちなし。
- 双対 CV 重心・双対面重心: 四面体/パッチ頂点の**算術平均の体積/面積加重和**であり、
  2D シューレースのような「相殺した和の比」を取らない — 桁落ちなし。

**露出している部分**:

1. **`newell()` (双対面/境界半割面の面積ベクトル)**: 項が `(a−b)·(a+b)` の形で、
   差分 (パッチ寸法 ~μm–100 μm) に**絶対座標の和 (~領域寸法)** が掛かる。座標は
   `nodes[].coords` = geom_float (float32) 格納なので丸め ~1e-9 m を持ち、これが絶対座標
   因子で増幅される。スリバーパッチ (例: 100 μm × 0.5 μm, |S|~5e-11 m²) では
   誤差 ~数 e-11 = **面ベクトルの数十 %** になり得る。向き整合の符号判定
   (`pv·eAB < 0` / `hv·S < 0`) が誤反転するリスクもある。演算自体は double なので、
   2D と同じく**パッチのローカル原点 (例: エッジ中点 M / 面重心 Fcen) へ平行移動してから
   Newell を評価**すれば増幅因子が領域寸法→パッチ寸法に落ち、~5 桁改善する。
2. **蓄積配列が geom_float**: `DFA` (双対面アキュムレータ) は double だが、
   `bnodeAccum` / `halfByOwner` / `hcentByOwner` は geom_float で加算 — スリバー面の
   小さい寄与同士の加算で相対誤差が積む。ローカル double 蓄積 → 最後に cast へ。
3. (確認事項) `planes[ip].centCoords` / `surfVect` (primal 面) は float32 で計算済みの値を
   使う — primal 面は通常寸法なので実害は小さい見込みだが、検証項目に含める。

## 3. スコープ

- **やる**: `buildMedianDual3D` 内の Newell 呼び出しのローカル原点化 (内部双対面パッチ
  (M,F1,G,F2) と境界半割面サブ四角 (N,Mn,Fcen,Mp) の両方)、境界半割面蓄積の double 化、
  閉性/体積診断の閾値確認。methods への精度注意の追記。
- **やらない**: 2D 側の再変更 (修正済み)。coords の double 格納化 (全体 FP64 化は
  [fp32-metric-freestream-fix] 系の別議論)。dcc/勾配など solver 側幾何の変更。

## 4. 設計方針

- `newell(p, out)` は現状シグネチャのまま、**呼び出し側で p を局所原点相対にする**
  (2D 修正と同じ流儀)。局所原点は当該パッチの第 1 頂点 (M / N) で十分。
  面ベクトルは平行移動不変なので結果は厳密演算で同値。
- 境界半割面: `hv` 計算を相対座標化し、`bnodeAccum`/`halfByOwner`/`hcentByOwner` を
  double 蓄積に変更、平坦化時に geom_float へ cast。
- 向き整合の符号判定は相対座標化後の `hv` で行う (精度改善後は誤反転リスク消滅)。
- 既存の自己診断 (`closure max|sum dS|`, `volume sum: dual= primal=`, 非正体積警告) を
  そのまま合否ゲートに使う。

## 5. 実装ステップ

1. methods/discretization.md の M4 節に「双対幾何は局所原点+double で評価する
   (float32 絶対座標の桁落ち対策, 2D §2.9 と同根)」を追記。
2. `buildMedianDual3D`: 内部双対面パッチの Newell を M 相対化。
3. 境界半割面: サブ四角の Newell を N 相対化 + 蓄積 double 化。
4. 検証 (§6) → 変更ログ更新 → accepted へ。

## 6. 検証

- [ ] **スリバー単体検査**: 押し出し 3D 薄壁構造格子 (例: 10 mm 級領域で壁第一間隔 1 μm,
  ヘックス) を生成・変換し、(a) `closure max|sum dS|` が修正前後で桁改善すること、
  (b) 双対体積和 = primal 体積和 (既存診断) の相対誤差 ≤1e-6、(c) 壁隣接ノードの
  `wall_dist` が解析値 Δy と ~1% で一致 (2D の実害と同じ観点)、(d) 双対重心が CV 内に
  収まること。
- [ ] **非退行**: 既存 3D M4 検証ケース (case/38 周期チャネル系 / 3D TG) のメッシュを
  再変換し、通常寸法 CV の幾何差が相対 ~1e-6 以下であること + 短尺 run (数百 step) の
  残差履歴が従来と表示桁一致。
- [ ] 3D 薄壁 node SST の段階起動スモーク (2D run_0034 のレシピ: 陰 cfl0.5 1次 →
  陰 cfl4 2次) が初手発散しないこと。

## 7. 変更ログ

- `2026-08-11` — 起票 (2D 版修正 commit 3cdddb9e のフォローアップ。現状分析まで実施、
  実装は未着手)。

# architecture-rans-sst — RANS SST 親計画

## メタ

- **area**: `architecture`
- **status**: `in_progress`
- **related_docs**:
  - [`docs/turbulence/theory.md`](../../docs/turbulence/theory.md)
  - [`docs/turbulence/implementation.md`](../../docs/turbulence/implementation.md)
- **related_plans**:
  - [`architecture-axisymmetric.md`](architecture-axisymmetric.md)
  - [`gpu-implicit-plan.md`](gpu-implicit-plan.md)
- **created**: `2026-06-02`
- **owner**: `Copilot`

## 1. 目的

forge の圧縮性 NS ソルバに、低 Re の **Menter SST** モデルを導入する。
本親計画では、まず既存の laminar / LES 経路を維持したまま、
`rho k`, `rho omega` を generic scalar transport 基盤の最初の 2 本として導入する。
初回マイルストーンは explicit の advection-only 輸送を成立させることとし、
その後に diffusion / source / 本格 SST 閉鎖へ広げる。

## 2. スコープ

### やる
- `docs/turbulence/` の theory / implementation を整備する。
- `solverConfig` に RANS SST 用の設定項目を追加する。
- generic scalar transport 基盤を追加する。
- 初回は explicit で `rho k`, `rho omega` の advection-only 輸送を追加する。
- `turbulent_viscosity_d_wrapper` を LES / RANS dispatcher に拡張する。
- turbulence scalar の境界条件・出力・残差監視を追加する。
- 既存の laminar / LES ケースで回帰確認を行う。
- 既存の axisymmetric nozzle 複製 run で scalar advection を確認する。

### やらない
- 軸対称 SST の diffusion / source / 幾何 source 対応 → 子 plan に分離する。
- SST の陰解法 Jacobian 連成 → 別 plan に分離する。
- 壁関数、遷移モデル、他の RANS モデル追加。

## 3. 関連 docs と前提

- 理論: [`docs/turbulence/theory.md`](../../docs/turbulence/theory.md)
- 実装: [`docs/turbulence/implementation.md`](../../docs/turbulence/implementation.md)
- 全体構成: [`docs/architecture/overview.md`](../../docs/architecture/overview.md)
- 軸対称の前提: [`docs/axisymmetric/theory.md`](../../docs/axisymmetric/theory.md)
  と [`architecture-axisymmetric.md`](architecture-axisymmetric.md)
- 検証運用: [`../forge-verification-cases.md`](../forge-verification-cases.md)

軸対称 plan 側で「LES / RANS の軸対称対応 → 別 plan」と明記されているため、
本計画で扱う axisymmetric は既存幾何の上で scalar advection plumbing を確認する
範囲に限定し、SST 特有の diffusion / source / 幾何項は完了条件に含めない。

## 4. 設計方針

- 保存変数は既存 5 変数に `rho k`, `rho omega` を追加するが、実装は
  「5 方程式 NS 本体」と「generic scalar transport 層」を分離して進める。
- 初期フェーズでは implicit 用 block system は拡張せず、explicit storage と
  residual 系のみを広げる。
- 渦粘性は既存 `vis_turb` を共用し、WALE と SST が同じ出力経路を使う。
- 乱流 scalar の対流・拡散は、既存 NS flux から分離した
  `scalarTransport_d.*` の共通輸送コアで扱う。
- SST 固有の boundary / source / closure は model-specific layer として分離し、
  `ransScalarBoundary_d.*`, `ransSource_d.*` のような専用 file 群に置く。
- 初回マイルストーンの scalar advection は 1 次 upwind を採用し、
  scalar gradient / limiter / diffusion / source は後段に送る。
- axisymmetric の SST 固有幾何 source は本親計画では扱わず、子 plan に切り出す。

## 5. 実装ステップ

1. `docs/turbulence/*` を追加し、`docs/index.md` に領域を追加する。
2. 本 plan を起票し、`.github/plans/README.md` に登録する。
3. `input/solverConfig.{hpp,cpp}` に `RANSmodel` と自由流乱流量設定を追加する。
4. `variables.{hpp,cpp}` に `rho k`, `rho omega`、関連勾配、残差、出力項を追加する。
5. `cuda_forge/turbulent_viscosity_d.cu` を SST dispatcher 対応に拡張する。
6. `scalarTransport_d.*` を追加し、まず advection-only の kernel / wrapper を実装する。
7. `ransScalarBoundary_d.*` を追加し、wall / inlet などの turbulence scalar 境界を既存 NS BC の後段で処理する。
8. `ransSource_d.*` を追加し、SST source を scalar transport core の後段に差し込む。
9. `boundaryCond.*`, `output/*`, `main.cpp` の残差監視を拡張する。
10. 既存 nozzle 複製 run を含む標準ケースで explicit scalar advection を検証する。
11. 完了後、軸対称 SST の子 plan を起票する。

## 6. 検証

- **単体 / ビルド**: `solver_density_cuda` がコンパイルできること。
- **検証ケース**:
  - 既定: `case/23.axi_nozzle/run_0033_slau_axisym_m4_publicrao_full_10k` の複製 run
  - 回帰: 既存 laminar / LES case
  - 後続: airfoil 系 case
- **判定基準**:
  - SST 無効時に既存 laminar / LES ケースが従来通り動く。
  - advection-only 段階で `k`, `omega` が downstream へ輸送され、負値化や発散を起こさない。
  - `residual_history.csv` と `residual_history.png` を新規 run に保存する。

## 7. 影響範囲

- 触る主要ファイル:
  - [`input/solverConfig.{hpp,cpp}`](../../solver_density_cuda/input/)
  - [`variables.{hpp,cpp}`](../../solver_density_cuda/)
  - [`main.cpp`](../../solver_density_cuda/main.cpp)
  - [`boundaryCond.{hpp,cpp}`](../../solver_density_cuda/)
  - `cuda_forge/turbulent_viscosity_d.*`
  - turbulence scalar 用の新規 kernel / wrapper
  - [`output/`](../../solver_density_cuda/output/)
- 既存ケースへの影響:
  - SST 無効時は完全互換を維持する。
- ドキュメント:
  - `docs/turbulence/*`
  - [`docs/index.md`](../../docs/index.md)
  - [`.github/plans/README.md`](README.md)

## 8. 完了条件

- [x] 関連 `docs/turbulence/theory.md` 追加済み
- [x] 関連 `docs/turbulence/implementation.md` 追加済み
- [ ] generic scalar transport 基盤の初回実装と advection-only 検証完了
- [ ] turbulence scalar boundary layer の分離完了
- [ ] SST source layer の分離完了
- [ ] `.github/plans/README.md` の状態を `done` に更新
- [ ] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載
- [ ] 軸対称 SST の子 plan を起票

## 9. 変更ログ

- `2026-06-02` — 初稿。親 plan と turbulence docs を作成。
- `2026-06-06` — 実装方針を更新。SST 専用 kernel 群ではなく generic scalar transport 基盤を採用し、初回マイルストーンを advection-only に再定義。
- `2026-06-06` — scalar transport に face ベース diffusion を追加し、`k` / `omega` の輸送を advection + diffusion の 2 段階で回す実装へ進めた。
- `2026-06-06` — 今後の凝縮・液滴・複数 scalar 拡張を見据え、`scalarTransport_d.*` は共通輸送コアに限定し、SST 固有の boundary / source / closure は model-specific layer に分離する方針へ更新した。
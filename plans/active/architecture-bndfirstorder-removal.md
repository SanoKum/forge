# `mesh.bndFirstOrder` の廃止 (使用禁止 → コード削除)

## メタ

- **area**: `architecture / boundary`
- **status**: `draft`  <!-- 使用禁止のルール化は完了、コード削除が残タスク -->
- **related_docs**:
  - [`procedures/solver-settings.md`](../../procedures/solver-settings.md) (「使用禁止」節が正)
  - [`methods/discretization.md`](../../methods/discretization.md) §境界弱形式 (旧 `bndFirstOrder` 診断の記述を要修正)
- **related_plans**:
  - [`boundary-node-nozzle-wall-outlet-stability.md`](boundary-node-nozzle-wall-outlet-stability.md) §2.7 (node ノズルレシピからの撤去は実施済み)
  - [`discretization-node-boundary-ghostless.md`](discretization-node-boundary-ghostless.md) (本フラグが代用していた「境界のロバスト化」の本筋)
- **created**: `2026-08-12`
- **owner**: `sano`

## 1. 目的

`mesh.bndFirstOrder` は「境界隣接 CV の 2 次 MUSCL 再構成だけを 1 次に落とす」意図の診断/安定化
フラグだったが、実装が意図と一致しておらず、**安定化した場合もその理由が交絡していて解釈できない**。
ユーザ指示 (2026-08-12) により**まず使用禁止をルール化**し (完了)、続いて**コードから削除**して
同種の誤用と誤った切り分けが再発しない状態を得る。

## 2. スコープ

- **やる**:
  - `bndFirstOrder` の config キー・`zeroBndNodeGradient_d` カーネル・呼び出しの削除。
  - 未知キーとして無視されず**config load でエラーにする**か、少なくとも警告する (既存 run の
    config に残っているため、黙って無視すると「効いているつもり」の再発を招く)。
  - `methods/discretization.md` の旧記述 (bump hiM が本フラグで直った、という記録) に
    「副作用による見かけの安定化」という訂正を入れる。
- **やらない**:
  - `bnode_flag_d` 自体の削除 (他用途があれば残す。ただし「いずれかの bcond の CV」という
    定義が疑似 2D で全域になる件は、利用箇所ごとに妥当性を確認する)。
  - 本フラグが代用していた近壁 2 次再構成のロバスト化そのもの
    (→ ghostless 境界弱形式 plan の担当)。

## 3. 廃止理由 (2026-08-12, case/26 の node 2 次切り分けで判明)

1. **粘性応力を破壊する**。`calcGradient_d.cu` の `zeroBndNodeGradient_d` は
   `drodx…dPdz` に加えて `dUxdx…dUzdz` を 0 にする。コメントは
   「divU・dT は粘性で使うため残す (本オプションは対流再構成のみを対象とする)」とあるが、
   **`viscousFlux_d.cu:114-118` は `dUxdx/dUxdy/dUydx/…` を面平均して Newton 応力を作る**ので、
   偏差応力が丸ごと消える。divU/dT を残しても救われない。
2. **適用範囲が意図より遥かに広い**。`mesh.cpp:876-884` の `bnode_flag` は
   *いずれかの* bcond の `iCells` を 1 にする。疑似 2D (1 セル厚) メッシュでは spanwise の
   `side1`/`side2` slip が**全ノードを覆う**ため flag が全域 1 になる
   (case/26 `flat_plate_yp30_node.h5`: 21510/21510 ノードが flag 1)。
   → 「境界近傍だけ 1 次化」ではなく **全域 1 次化 + 全域の粘性応力喪失**。
3. 結果として、本フラグを付けた A/B は**切り分けとして成立しない**。case/26 `run_0037` の
   「bndFirstOrder でも残差上昇 (改善せず)」は、2 次精度の不安定が直らなかったのではなく
   **粘性を失った別の計算**を見ていた可能性が高い。

過去に本フラグで「安定化した」と記録されたケース (bump hiM node / 旧 node ノズルレシピ) も、
上記の副作用 (実質 1 次化 + 粘性減) による見かけの安定化と解釈すべきで、正しい対策ではない。

## 4. 実装ステップ

1. `solver_density_cuda/input/solverConfig.{hpp,cpp}`: `bndFirstOrder` メンバとパースを削除し、
   キーが存在したら明示エラー (「廃止されたキー」) にする。
2. `solver_density_cuda/cuda_forge/calcGradient_d.cu`: `zeroBndNodeGradient_d` と呼び出し
   (L1015-1018 付近) を削除。
3. `solver_density_cuda/mesh/mesh.{cpp,hpp}`: `bnode_flag_d` の他用途を確認。用途が
   本フラグだけなら削除、残す場合は「いずれかの bcond」という定義の危険性をコメントに明記。
4. docs: `procedures/solver-settings.md` (完了)、`AGENTS.md` (完了)、
   `methods/discretization.md` の旧記述訂正、既存 `case/*/README.md` の該当行に注記。

## 5. 検証

- **ビルド**: full rebuild (config 構造体を触るため差分ビルド禁止 — [[stale-build-struct-layout-trap]])。
- **回帰**: `bndFirstOrder` を使っていない既定経路がビット不変であること
  (case/40 `run_0066` 生産構成の再走で `τ_w`/壁温/ṁ/推力が再現性帯内)。
- **エラー**: 旧キーを含む config を食わせて明示エラーになること。

## 6. 完了条件

- [x] 使用禁止のルール化 (`AGENTS.md` + `procedures/solver-settings.md`)
- [ ] コード削除 + 旧キーの明示エラー化
- [ ] `methods/discretization.md` の旧記述訂正
- [ ] 既定経路のビット不変回帰
- [ ] `status: done` にして `plans/accepted/` へ移動、[`plans/README.md`](../README.md) 同期

## 7. 変更ログ

- `2026-08-12` — 起票。ユーザ指示により使用禁止をルール化 (`AGENTS.md` / `procedures/solver-settings.md`)。
  コード削除は将来課題として本 plan に保持。廃止理由 §3 は case/26 の node 2 次切り分けで実測確認。

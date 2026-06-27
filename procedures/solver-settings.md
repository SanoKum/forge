# Forge ソルバー設定リファレンス

> **位置づけ**: 本文書は `solverConfig.yaml` の**運用上の設定リファレンス** (各設定の意味・正しい使い方・注意点) である。
> スキームの理論的背景や実装の詳細は `methods/<area>/` を参照すること (例: 対流スキームは [`../methods/convection/`](../methods/convection/)、
> 時間積分は [`../methods/time_integration/`](../methods/time_integration/))。ここは「どの値を選ぶか」、`methods/` は「なぜそう動くか」を担う。

この文書は `solverConfig.yaml` の主要な数値設定について、正しい使い方と注意点をまとめたものです。解析設定を変更するときは、**必ずこのファイルを参照してから**変更を行うこと。

## convMethod — 対流スキームの次数

| 値 | スキーム | 説明 |
|---|---|---|
| `0` | 1次風上 | 安定性が高い。収束困難なケースの初期化や発散回避に使う |
| `1` | 2次風上 | 通常の高精度計算に使う |
| `2` | 3次 MUSCL | より高精度だが計算コスト増。安定した解が得られてから使う |

## limiter — リミタ

| 値 | リミタ | 説明 |
|---|---|---|
| `0` | なし | **リミタ無しのため最も発散しやすい。** 基本的に単体では使わない |
| `1` | Barth | |
| `2` | Venkatakrishnan | 通常使用するリミタ。`convMethod: 1` または `2` と組み合わせる |

## 推奨組み合わせ

| 用途 | convMethod | limiter |
|---|---|---|
| 初期化・発散回避 | `0` | `2` |
| 通常の高精度計算 | `1` | `2` |
| 高精度（安定解から） | `2` | `2` |

**`limiter: 0` は原則使用しない。** 意図的に使う場合は理由を明記すること。

## CFL の定義 (`cfl` と `cfl_pseudo`) — 重要

`time.deltaT` には `cfl` と `cfl_pseudo` の 2 つがある。**定常計算 (`unsteady: 0`) と
二重時間刻み (`dualTime: 1`) では、実際の局所時間刻みを決めるのは `cfl` ではなく
`cfl_pseudo` である。** これを取り違えないこと。

根拠 ([`setDT_d.cu`](../solver_density_cuda/cuda_forge/setDT_d.cu) の `setDT_d_wrapper`):

- `unsteady == 0` または `dualTime == 1` のとき、局所刻みは `setDTlocal_pseudo_cell_d` で
  `dt_local = cfl_pseudo * dt / cfl_cell` と設定される。`cfl_cell` は同じ `dt` から計算されるため
  `dt` が相殺し、**実効局所 CFL = `cfl_pseudo`** になる。
- `cfg.cfl` は表示用の `cfg.dt` (および `max cfl` 表示) をスケールするだけで、
  定常計算の積分自体には効かない。したがって定常で `cfl` だけを上げても収束は速くならない。

運用ルール:

- **定常 (`unsteady: 0`) の収束を速めたいときは `cfl_pseudo` を上げる** (例: 陰解法 `blockDPLUR: 1`
  なら `cfl_pseudo` を 20〜50 程度まで)。`cfl` を上げても無意味。
- 非定常 (`unsteady: 1`, `dualTime: 0`) の物理時間刻みは `cfl` (または `dt`) で決まる。
- ログの `max cfl` 表示は `cfg.cfl` に追従する値であり、実効積分 CFL ではない点に注意する。

## detectNaN — NaN/Inf 検知診断モード

`time.deltaT.detectNaN` (任意, 既定 `0`)。発散の発生ステップと場を特定するための診断オプション。

- `0` (既定): 検査を一切行わない。通常実行はビット不変・性能影響なし。
- `1`: **毎ステップ終端**で保存量 `ro`,`roUx`,`roUy`,`roUz`,`roe` と圧力 `P` (RANS 時は `roK`,`roOmega` も) を内部セル (`nCells`) にわたって検査する。NaN/Inf が 1 つでもあれば、その時点の解を `res_nan_<step>.h5` / `.xmf` に強制ダンプし、最初に該当した変数名を `stderr` に出力して `exit(1)` で停止する。

```yaml
time:
  deltaT:
    detectNaN: 1   # 発散をその場で捕捉してダンプ・停止する
```

毎ステップ GPU リダクション (`thrust::any_of`) が走るため、原因調査時のみ有効化し、本番計算では `0` のままにする。ダンプされた `res_nan_<step>.h5` を ParaView で開けば、どの境界・どのセルで最初に NaN が出たかを直接確認できる。

## リスタート手順

前の run の結果から引き継ぐ場合は `valueFileName` に該当の `res_*.h5` を指定する。

```yaml
mesh:
  meshFileName: "axi_nozzle.h5"      # メッシュ (変わらない)
  valueFileName: "res_10000.h5"      # リスタート元の結果ファイル
```

`initial` キーの値はリスタート時には無視される（valueFileName が res_*.h5 を指していれば、そこから値を読む）。

## 段階的な計算戦略

均一初期値から超音速流れを計算する場合、いきなり2次以上で始めると全域超音速流出に発散することがある。次の2段階を踏む：

1. `convMethod: 0, limiter: 2` で安定した解を得る（例：10000ステップ）
2. その結果を `valueFileName` で引き継ぎ、`convMethod: 1, limiter: 2` で計算する

## discretization / bndFirstOrder — 離散化レイアウト (node-centered)

`mesh.discretization` (任意, 既定 `"cell"`)。`"node"` で node-centered (中点双対 median-dual) 化。
詳細は [`methods/discretization/`](../methods/discretization/)。

`mesh.bndFirstOrder` (任意, 既定 `0`)。`1` で境界隣接 CV の 2 次 MUSCL 再構成を 1 次に落とす。
node-centered の壁近傍高マッハ発散 (近壁 2 次再構成のロバスト性問題) の対策。explicit では
リミットサイクルが残ることがあり、implicit と併用すると完全収束する (bump Mach1.65 で確認)。
cell-centered では通常不要 (既定 0)。

```yaml
mesh:
  discretization: "node"   # cell | node
  bndFirstOrder: 1         # node-centered の壁近傍高マッハ安定化
```

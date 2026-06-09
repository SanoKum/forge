# Forge Calculation Workflow

この文書は、このリポジトリで計算ケースを準備して `forge` を実行する際の標準手順をまとめたものです。計算方法を案内するときは、まずこの文書の流れを基準にし、ケース側に専用スクリプトや README がある場合はそちらを優先します。

## 対象範囲

`case/` 配下の多くのケースは、概ね次の構成を取ります。

- `mesh/`: `.geo` や生成スクリプトなど、メッシュ作成用の入力
- `run_*`: `run_slau`、`run_roe`、`run_keep` など、ソルバー実行用ディレクトリ

質問が特定ケースに限定されていない場合も、基本的にはこの構成を前提に案内します。

## 標準手順

計算を実行するときは、既存の `run_*` ディレクトリをそのまま使い回さず、必ず新しい `run_*` ディレクトリを作ってから実行する。

理由:

- 既存の結果やログを上書きしないため
- 実行時の設定差分を追跡しやすくするため
- `residual_history.csv`、`res_*.h5`、実行ログなどを run ごとに保存するため
- 実行後に residual plot PNG を run ごとに残すため

## run ディレクトリ命名規則

`run_*` が増えたときに新旧や派生関係が分からなくならないよう、run 名は次の形式で統一する。

```text
run_<NNNN>_<scheme>_<short-slug>
```

例:

```text
run_0001_slau_baseline
run_0002_slau_lowbp_1st
run_0003_slau_lowbp_muscl_restart
run_0004_roe_baseline
```

ルール:

- `<NNNN>` は `0000` から `9999` までの 4 桁ゼロ埋め連番とする。
- `<scheme>` は `slau`, `roe`, `keep`, `ausm`, `ausmUP` などスキーム識別子を入れる。
- 連番はケースディレクトリごとに単調増加で使い、スキームを変えても振り直さない。
- 欠番は気にせず、既存最大番号の次を使う。
- `<short-slug>` は 2〜4 語程度の短い ASCII ラベルにし、変更意図を一目で分かるようにする。
- 日付は原則 run 名に入れず、時系列は連番で管理する。
- 既存 run を複製して派生 run を作るときも、新しい連番を必ず振り直す。

推奨運用:

- まず `ls` や `find` でケース内の既存 `run_` を確認する。
- その中の最大番号に `+1` した名前で新しい run を作る。
- `short-slug` は `baseline`, `lowbp`, `lowbp_1st`, `muscl_restart`, `axisym_m2` のように目的中心で付ける。

非推奨:

- `run_slau_test`, `run_new_slau`, `run_latest_slau` のような意味が時間で変わる名前
- `run_20260531_slau_xxx` のように日付だけで順序管理する名前
- `run_0007_slau_copy`, `run_0008_slau_final2` のような派生関係が読めない名前

1. 対象のケースディレクトリを `case/` 配下から特定する。
2. 基準にする既存の `run_*` を複製して、新しい `run_*` ディレクトリを作る。
3. 必要なら新しい `run_*` 内の `solverConfig.yaml` などを調整する。
4. `mesh/` に入り、Gmsh で `.geo` から `.msh` を生成する。
5. 生成した `.msh` を新しい `run_*` ディレクトリへコピーする。
6. その新しい `run_*` ディレクトリで `convertGmshToForge` を実行し、HDF5 に変換する。
7. 実行する `forge` バイナリが最新ソースで再ビルドされているか確認する (下記「実行前のバイナリ鮮度チェック」)。
8. 同じ新しい `run_*` ディレクトリで `forge` を実行する。
9. 実行後、その新しい `run_*` ディレクトリで `residual_history.csv` から `residual_history.png` を生成する。

重要なのは、`solverConfig.yaml` がある新しい `run_*` ディレクトリを起点に変換と実行を行うことです。

## 実行前のバイナリ鮮度チェック (重要)

`forge` を実行するコマンドは `solver_density_cuda/build/forge` を直接呼ぶだけで、**ソースの再ビルドを自動では行わない**。`cuda_forge/*.cu` などを編集した後に古い `build/forge` のまま計算すると、現行ソースと無関係な **stale バイナリの結果が silently に出る** (例: ソース側で実装した乱流モデルが反映されず、結果が壊れて見える)。

`forge` を実行する前に、必ずバイナリがソースより新しいことを確認すること。

```bash
cd /home/sano/work/forge/solver_density_cuda
# build/forge より新しいソースがあれば stale (= 再ビルドが必要)
find cuda_forge *.cpp *.cu *.hpp -newer build/forge 2>/dev/null
```

上のコマンドが何か出力したら、計算前に再ビルドする。Docker dev image 内で `tools/build.sh` を使う:

```bash
cd /home/sano/work/forge
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/workspace" \
  -e HDF5_INC=/usr/include/hdf5/serial \
  -e HDF5_LIBDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial \
  forge-solver:cuda-dev \
  bash -c "cd /workspace/solver_density_cuda && bash tools/build.sh"
```

設定変更 (例: `convMethod`) が解に与える影響を診断するときは、まず stale でないバイナリで「変更前設定」を再現する control run を取り、スキーム差とコード世代の差を取り違えないこと。

## 標準コマンド例

```bash
# 既存 run を複製して新しい run を作る
cp -r run_0003_slau_baseline run_0004_slau_implicit_check

# <case>/mesh で実行
gmsh -3 input.geo -o case.msh -format msh4

# 生成メッシュを run ディレクトリへコピー
cp case.msh ../run_0004_slau_implicit_check/

# <case>/run_0004_slau_implicit_check で実行
convertGmshToForge case.msh case.h5
# forge 実行前に build/forge が stale でないか確認 (「実行前のバイナリ鮮度チェック」参照)
forge
python3 /home/sano/work/forge/solver_density_cuda/tools/plot_implicit_residuals.py \
	--input residual_history.csv \
	--output residual_history.png
```

Docker 経由で上記を手動実行する場合は、必ず `--user "$(id -u):$(id -g)"` を付けること (詳細は `.github/forge-development-environment.md` のファイル所有者ルールを参照)。

読み替えポイント:

- `input.geo`: 実際のジオメトリ入力ファイル名
- `case.msh`: 生成するメッシュファイル名
- `run_0004_slau_implicit_check`: 実際に使う新しい `run_*` ディレクトリ名
- `case.h5`: `solverConfig.yaml` の `mesh.meshFileName` と `mesh.valueFileName` に一致する HDF5 名

## Docker を使う場合の流れ

`solver_density_cuda` の Docker 環境を使う場合も、計算の論理的な流れは同じです。

1. 必要なら Docker イメージをビルドする。
2. 基準にする既存の `run_*` を複製して、新しい `run_*` を作る。
3. ケースの `mesh/` で `.msh` を生成する。
4. `.msh` を新しい `run_*` にコピーする。
5. その新しい `run_*` ディレクトリで `convertGmshToForge` を実行する。
6. 同じ新しい `run_*` ディレクトリで `forge` を実行する。
7. 実行後、その新しい `run_*` ディレクトリで `residual_history.csv` から `residual_history.png` を生成する。

代表例:

```bash
cd /home/sano/work/forge/solver_density_cuda
docker build -f Dockerfile.cuda.dev -t forge-solver:cuda-dev .

cd /home/sano/work/forge/case/20.naca_ml/001.test
cp -r run_slau run_slau_20260510_implicit_check

cd /home/sano/work/forge/case/20.naca_ml/001.test/run_slau_20260510_implicit_check
./run_case.sh
```

ケースに `run_case.sh` のような専用ヘルパーがある場合は、手順を手で再構成するより、そのヘルパーを優先します。ただし、その場合も既存の `run_*` を直接使わず、複製した新しい `run_*` 側で実行する。

Residual plot は PNG で残せば十分とし、既定では `residual_history.csv` を入力に `residual_history.png` を生成する。

## Gmsh / ParaView の起動

Docker 環境経由で Gmsh や ParaView を開くときは、`solver_density_cuda/tools/` 配下のホスト側ラッパースクリプトを使います。

Gmsh:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_gmsh_gui.sh
```

ParaView:

```bash
cd /home/sano/work/forge/solver_density_cuda
./tools/run_paraview_gui.sh
```

ファイルを直接開く例:

```bash
./tools/run_gmsh_gui.sh /home/sano/work/forge/case/08.bump/mesh/bump.geo
./tools/run_paraview_gui.sh /home/sano/work/forge/case/08.bump/run_slau/res_0.xmf
```

補足:

- これらのスクリプトはコンテナ内ではなくホスト側から実行する。
- WSL2 では WSLg を前提に GUI 転送する。
- ParaView のアイコン表示が崩れる場合は `solver_density_cuda/Dockerfile.cuda.dev` から Docker イメージを再ビルドする。

## ケースごとの注意

- `run_case.sh` などのケース専用ヘルパーがあるなら、それを優先する。
- 実行前に、既存の `run_*` を複製して新しい `run_*` を作る。既存ディレクトリへそのまま再実行しない。
- HDF5 変換後のファイル名は、必ず `solverConfig.yaml` の `mesh.meshFileName` と `mesh.valueFileName` に合わせる。
- `forge` は必ず `solverConfig.yaml` がある新しい `run_*` ディレクトリで実行する。
- `forge` 実行後は、同じ新しい `run_*` ディレクトリで `residual_history.csv` から `residual_history.png` を生成する。
- `run_*` が複数ある場合は、ユーザー指定のもの、または使いたいスキームに対応するものを選ぶ。

## 結果ディレクトリの明示 (投入時・まとめ時)

検証計算の所在が後から分からなくならないよう、計算の **投入時** と **結果まとめ時** の両方で、出力先 `run_*` をリポジトリルートからの相対パスで明示する (正本ルールは `AGENTS.md` の「計算結果ディレクトリの明示」)。

- **投入時**: 起動する応答内で各 run の出力先パスを列挙し、run ごとの設定差分 (スキーム・`convMethod`・物性など) を 1 行で添える。
  - 例: `case/05.sod_shock_tube/run_0004_slau_tracer_n2n2/` — SLAU, `species:[N2,N2]` 識別子トレーサ。
- **まとめ時**: 結論・比較表とあわせて、根拠となった `run_*` パスと主要成果物 (`res_*.h5`, `residual_history.png`, 後処理図) を併記する。複数 run 比較では「どの数値/図がどの run か」を取り違えないよう run パスを明記する。
- 一時 run・破棄予定 run もパスと用途を明記し、所在不明にしない。

## 運用上の位置づけ

Copilot が計算手順を説明するときは、この文書を標準の参照先とします。個別ケースの検証方針や、Docker と native 実行の使い分けは、別の関連文書を優先して参照します。
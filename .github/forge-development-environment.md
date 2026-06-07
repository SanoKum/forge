# Forge Development Environment

この文書は、Forge の開発時に Docker と native 環境をどう使い分けるかを整理したものです。

## 基本方針

- 通常のコード修正、ビルド、ケース実行は Docker コンテナを基本とする。
- `solver_density_cuda` の開発環境は、`Dockerfile.cuda.dev` とそのコンテナ実行を既定の流れとして扱う。
- コンテナ内でビルドした成果物を使ってケースを実行する運用を標準とする。

## Docker を既定とする理由

- リポジトリのビルド手順や依存関係が Docker 前提で揃っている。
- 開発者ごとの差分を減らしやすい。
- `run_case.sh` のようなケース専用スクリプトも、Docker ベースの実行を前提にしているものがある。

## 例外: NVIDIA の計測ツールは native を優先する

NVIDIA 提供の GPU 速度計測・プロファイリングツールを使う場合は、Docker 内でうまくいかないことがあるため、WSL native または Linux native の環境を優先する。

対象として想定するツール:

- Nsight Compute (`ncu`)
- Nsight Systems (`nsys`)
- その他、NVIDIA 提供の GPU 計測ツール

通常の solver 実行は Docker を基準にしてよいが、これらのツールで計測したいときは、まず native 実行へ切り替える。

## native 実行時の目安

- WSL native / Linux native では `solver_density_cuda/tools/build_native_wsl.sh` をビルドの起点として扱う。
- Ubuntu 24.04 + CUDA toolkit 12.0 の実績では、CUDA host compiler に `gcc-12` / `g++-12` を使う構成が安定している。
- GPU アーキテクチャは環境に応じて `FORGE_CUDA_ARCHITECTURES` で上書きする。

代表例:

```bash
cd /home/sano/work/forge/solver_density_cuda
FORGE_BUILD_JOBS=1 FORGE_CUDA_ARCHITECTURES=86 ./tools/build_native_wsl.sh
```

この例は native 側のビルド確認用であり、通常開発の既定経路を Docker から native に置き換える意図ではない。

## Docker 実行時のファイル所有者ルール

`docker run` を含むすべてのスクリプト・コマンドには、必ず `--user "$(id -u):$(id -g)"` を付加すること。

```bash
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$repo_root":/workspace \
  ...
```

**理由:** Docker コンテナはデフォルトで root として実行される。`--user` を省略すると、コンテナ内で生成したファイル (`*.h5`, `*.msh`, `*.csv`, `*.png` など) がホスト側でも root 所有になり、通常ユーザーから読み書きできなくなる。

**適用範囲:**
- `case/` 配下の `run_case.sh`
- `solver_density_cuda/tools/` 配下のラッパースクリプト (既適用)
- エージェントや手動で書く `docker run` 呼び出しすべて

コンテナ内でツールが `$HOME` を参照する場合は、`-e HOME=/tmp/forge-home` も合わせて指定する (ParaView など)。

## 使い分けの実務ルール

- ビルドと通常のケース確認: まず Docker で行う。
- Gmsh / ParaView の GUI 起動: `solver_density_cuda/tools/` 配下のラッパースクリプトをホスト側から使う。
- NVIDIA の計測ツールを試す: Docker 内に固執せず、WSL native または Linux native に切り替える。
- Docker で計測が失敗した場合: 例外扱いではなく、native へ切り替えるのが既定対応。

## 関連文書

- 計算準備と `forge` 実行手順: `.github/forge-calculation-workflow.md`
- 標準検証ケース: `.github/forge-verification-cases.md`
- Docker 開発環境の補足: `solver_density_cuda/README_docker.md`
- native ビルドスクリプト: `solver_density_cuda/tools/build_native_wsl.sh`

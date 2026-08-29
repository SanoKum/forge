#!/usr/bin/env bash
# EC2 インスタンス初期セットアップ (Deep Learning Base GPU AMI, Ubuntu 22.04 想定)。
# SSH ログイン後に一度だけ実行する:
#   curl -fsSL https://raw.githubusercontent.com/SanoKum/forge/<branch>/solver_density_cuda/tools/cloud/setup_instance.sh | bash
# もしくは clone 済みなら: ./solver_density_cuda/tools/cloud/setup_instance.sh
#
# 環境変数:
#   FORGE_REPO_DIR   clone 先 (default: $HOME/forge)
#   FORGE_BRANCH     clone するブランチ (default: main)
#   FORGE_REPO_URL   リポジトリ URL (default: git@github.com:SanoKum/forge.git)
set -euo pipefail

repo_dir="${FORGE_REPO_DIR:-$HOME/forge}"
branch="${FORGE_BRANCH:-main}"
repo_url="${FORGE_REPO_URL:-git@github.com:SanoKum/forge.git}"

echo "== [1/5] GPU / Docker の確認 =="
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader
docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi -L

echo "== [2/5] リポジトリ取得 =="
if [[ ! -d "$repo_dir/.git" ]]; then
  if ! ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new -T git@github.com 2>&1 | grep -q "successfully authenticated"; then
    echo "ERROR: GitHub への SSH 認証に失敗。deploy key を ~/.ssh/ に配置し、" >&2
    echo "       ~/.ssh/config で Host github.com に IdentityFile を指定してから再実行する。" >&2
    exit 1
  fi
  # .git が ~900MB あるため shallow clone を既定とする
  git clone --depth 1 --single-branch -b "$branch" "$repo_url" "$repo_dir"
fi
git -C "$repo_dir" submodule update --init --depth 1 solver_density_cuda/third_party/HighFive

echo "== [3/5] 計算用 Docker イメージのビルド =="
docker build -f "$repo_dir/solver_density_cuda/Dockerfile.cuda.cloud" \
  -t forge-solver:cuda-cloud "$repo_dir/solver_density_cuda"

echo "== [4/5] コンテナ内ビルド (Release) =="
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$repo_dir/solver_density_cuda":/workspace -w /workspace \
  forge-solver:cuda-cloud ./tools/build.sh

echo "== [5/5] idle auto-stop の常駐化 (消し忘れ保険) =="
sudo install -m 0755 "$repo_dir/solver_density_cuda/tools/cloud/idle_autostop.sh" /usr/local/bin/forge-idle-autostop
echo "*/5 * * * * root /usr/local/bin/forge-idle-autostop >> /var/log/forge-idle-autostop.log 2>&1" \
  | sudo tee /etc/cron.d/forge-idle-autostop >/dev/null

echo "DONE. forge binary: $repo_dir/solver_density_cuda/build/forge"
echo "速度計測用の native ビルドは: cd $repo_dir/solver_density_cuda && ./tools/build_native_wsl.sh"

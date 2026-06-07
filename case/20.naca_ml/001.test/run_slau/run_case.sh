#!/usr/bin/env bash
set -euo pipefail

image_name="${FORGE_IMAGE_NAME:-forge-solver:cuda-dev}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case_dir="$(cd -- "$script_dir/.." && pwd)"
repo_root="$(cd -- "$case_dir/../../.." && pwd)"
mesh_dir="$case_dir/mesh"

if ! command -v docker >/dev/null 2>&1; then
  echo "docker command not found" >&2
  exit 1
fi

if ! docker image inspect "$image_name" >/dev/null 2>&1; then
  echo "docker image '$image_name' not found" >&2
  echo "Build it first from solver_density_cuda:" >&2
  echo "  docker build -f Dockerfile.cuda.dev -t $image_name ." >&2
  exit 1
fi

echo "[1/4] Generate naca.geo"
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$repo_root":/workspace \
  -w /workspace/case/20.naca_ml/001.test/mesh \
  "$image_name" \
  bash -lc 'python3 generate_naca_v3.py'

echo "[2/4] Generate naca.msh with gmsh"
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$repo_root":/workspace \
  -w /workspace/case/20.naca_ml/001.test/mesh \
  "$image_name" \
  bash -lc 'gmsh -3 naca.geo -o naca.msh -format msh4'

cp "$mesh_dir/naca.msh" "$script_dir/naca.msh"

echo "[3/4] Convert naca.msh to naca.h5"
docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$repo_root":/workspace \
  -w /workspace/case/20.naca_ml/001.test/run_slau \
  "$image_name" \
  bash -lc '/workspace/solver_density_cuda/build/convertGmshToForge naca.msh naca.h5'

echo "[4/4] Run forge"
exec docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$repo_root":/workspace \
  -w /workspace/case/20.naca_ml/001.test/run_slau \
  "$image_name" \
  bash -lc '/workspace/solver_density_cuda/build/forge'
#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
solver_dir="$(cd -- "$script_dir/.." && pwd)"

missing_packages=()

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

check_package() {
  local pkg="$1"
  if ! dpkg -s "$pkg" >/dev/null 2>&1; then
    missing_packages+=("$pkg")
  fi
}

echo "[check_wsl_profile_target] solver dir: $solver_dir"
echo "[check_wsl_profile_target] This checks whether WSL is ready for Windows-side Nsight Compute profiling."
echo

if have_cmd nvcc; then
  echo "[ok] nvcc: $(command -v nvcc)"
  nvcc --version | tail -n 1
else
  echo "[missing] nvcc not found"
  echo "          Install the CUDA toolkit for WSL-Ubuntu. Do not install cuda/cuda-drivers meta-packages in WSL."
fi

for tool in cmake ninja g++ ssh; do
  if have_cmd "$tool"; then
    echo "[ok] $tool: $(command -v "$tool")"
  else
    echo "[missing] $tool"
  fi
done

if have_cmd sshd; then
  echo "[ok] sshd: $(command -v sshd)"
else
  echo "[missing] sshd"
fi

check_package build-essential
check_package cmake
check_package ninja-build
check_package gfortran
check_package pkg-config
check_package git
check_package python3
check_package python3-pip
check_package python3-h5py
check_package libhdf5-dev
check_package libyaml-cpp-dev
check_package libmetis-dev
check_package openssh-server

echo
if (( ${#missing_packages[@]} == 0 )); then
  echo "[ok] apt packages: all expected packages are installed"
else
  echo "[missing] apt packages: ${missing_packages[*]}"
  echo "          Suggested command:"
  echo "          sudo apt update && sudo apt install -y ${missing_packages[*]}"
fi

echo
if [[ -e /usr/lib/wsl/lib/libcuda.so ]]; then
  echo "[ok] WSL CUDA driver stub: /usr/lib/wsl/lib/libcuda.so"
else
  echo "[warn] WSL CUDA driver stub not found at /usr/lib/wsl/lib/libcuda.so"
fi

if [[ -f "$solver_dir/third_party/HighFive/include/highfive/highfive.hpp" ]]; then
  echo "[ok] HighFive headers are present in third_party/HighFive"
else
  echo "[missing] HighFive headers not found in third_party/HighFive"
  echo "          Run: git submodule update --init --recursive"
fi

echo
echo "Next steps:"
echo "  1. Fix any [missing] items above."
echo "  2. Build natively with ./tools/build_native_wsl.sh"
echo "  3. Enable and start sshd in WSL so Windows Nsight Compute can use SSH remote launch."

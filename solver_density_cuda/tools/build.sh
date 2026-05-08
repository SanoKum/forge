#!/usr/bin/env bash
set -euo pipefail
# Simple helper invoked INSIDE the container.
# It passes the HDF5 include/lib to CMake to satisfy the current CMakeLists variables.

if [[ -f build/CMakeCache.txt ]] && ! grep -Fxq "CMAKE_HOME_DIRECTORY:INTERNAL=$PWD" build/CMakeCache.txt; then
	rm -f build/CMakeCache.txt
	rm -rf build/CMakeFiles
fi

mkdir -p build
cd build
cmake -G Ninja \
	-DCMAKE_BUILD_TYPE=Release \
	-Dhdf5_inc="$HDF5_INC" \
	-Dhdf5_libdir="$HDF5_LIBDIR" \
	..
cmake --build . -j

#!/usr/bin/env bash
set -euo pipefail
cmake -S . -B build-conda -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_PREFIX_PATH="$PREFIX" -DPHYLTER2_OPENMP=OFF
cmake --build build-conda --parallel "${CPU_COUNT:-2}"
ctest --test-dir build-conda --output-on-failure
cmake --install build-conda --component Runtime

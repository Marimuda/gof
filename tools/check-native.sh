#!/usr/bin/env bash
set -euo pipefail
gof_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
gof_standard=${1:-c++20}
case "$gof_standard" in c++20|c++23) ;; *) echo 'Use c++20 or c++23' >&2; exit 2 ;; esac
gof_cxx=${CXX:-c++}
gof_tmp=$(mktemp -d)
trap 'rm -rf -- "$gof_tmp"' EXIT
"$gof_cxx" -std="$gof_standard" -O1 -g -fsanitize=address,undefined \
    -fno-omit-frame-pointer -I "$gof_root/inst/include" \
    "$gof_root/tests/native/helpers.cpp" "$gof_root/src/utils.cpp" \
    "$gof_root/src/cumres.cpp" -larmadillo -o "$gof_tmp/native"
# LSan cannot inspect threads under the traced execution environment. ASan
# bounds/use-after-free checks and UBSan remain enabled, with fatal failures.
ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$gof_tmp/native"

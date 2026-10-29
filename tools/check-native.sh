#!/usr/bin/env bash
set -euo pipefail
gof_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
gof_standard=${1:-c++20}
case "$gof_standard" in c++14|c++17|c++20|c++23) ;; *) echo 'Use c++14, c++17, c++20 or c++23' >&2; exit 2 ;; esac
gof_cxx=${CXX:-c++}
gof_leaks=${GOF_CHECK_LEAKS:-0}
case "$gof_leaks" in 0|1) ;; *) echo 'GOF_CHECK_LEAKS must be 0 or 1' >&2; exit 2 ;; esac
gof_tmp=$(mktemp -d)
trap 'rm -rf -- "$gof_tmp"' EXIT
"$gof_cxx" -std="$gof_standard" -O1 -g -fsanitize=address,undefined \
    -fno-omit-frame-pointer -I "$gof_root/inst/include" \
    "$gof_root/tests/native/helpers.cpp" "$gof_root/src/utils.cpp" "$gof_root/src/glm.cpp" \
    "$gof_root/src/cumres.cpp" -larmadillo -o "$gof_tmp/native"
# LSan cannot inspect threads in the traced local environment; untraced CI
# enables it with GOF_CHECK_LEAKS=1. ASan/UBSan errors remain fatal in both modes.
ASAN_OPTIONS="${ASAN_OPTIONS:+$ASAN_OPTIONS:}detect_leaks=$gof_leaks:halt_on_error=1" \
UBSAN_OPTIONS="${UBSAN_OPTIONS:+$UBSAN_OPTIONS:}halt_on_error=1:print_stacktrace=1" "$gof_tmp/native"

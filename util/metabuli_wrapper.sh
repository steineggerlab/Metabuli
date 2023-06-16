#!/bin/sh
FLAGS="$(grep -m 1 '^flags' /proc/cpuinfo)"
case "${FLAGS}" in
  *avx2*)
    exec /usr/local/bin/metabuli_avx2 "$@"
    ;;
  *)
    exec /usr/local/bin/metabuli_sse2 "$@"
    ;;
esac

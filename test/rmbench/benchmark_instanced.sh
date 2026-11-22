#!/bin/bash
# Benchmark instanced (materialised and native compressed) against the
# conventional redrawn machine solve, reporting wall time and memory.
#
# Usage:
#   benchmark_instanced.sh [--repeats N] [--pin CORE] [--full] [--quick]
#                          [--with-conventional] [--tiled FILE] [--redraw FILE]
#
# Each method runs in its own process so the reported peak RSS is not
# contaminated by an earlier method. --quick (the default) shrinks the checked-in
# radial-machine tiled fixture to a 60-degree sector by dividing every tile's
# repeat count; --full uses the whole 360-degree machine.
#
# The materialised and native rows use the same instanced topology and differ
# only in whether the mesh is expanded (materialised) or consumed through the
# compact LogicalMeshView (native). The conventional row is the redrawn .fem
# reference and is included only with --with-conventional.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$(dirname "$HERE")")"
BIN="$ROOT/cfemm/bin/instanced_solver_benchmark"
TILED="$ROOT/mfemm/testing/radial_machine/data/radial_machine_tiled.json"
REDRAW="$ROOT/mfemm/testing/radial_machine/data/radial_machine_redraw_01.fem"

NREPS=2
DIVISOR=6
PIN=""
WITH_CONVENTIONAL=0

while [ $# -gt 0 ]; do
    case "$1" in
        --repeats) NREPS="$2"; shift 2 ;;
        --pin) PIN="$2"; shift 2 ;;
        --full) DIVISOR=1; shift ;;
        --quick) DIVISOR=6; shift ;;
        --with-conventional) WITH_CONVENTIONAL=1; shift ;;
        --tiled) TILED="$2"; shift 2 ;;
        --redraw) REDRAW="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [ ! -x "$BIN" ]; then
    echo "benchmark binary not found: $BIN" >&2
    echo "build it with: cmake --build <build-dir> --target instanced_solver_benchmark" >&2
    exit 1
fi

run_method() { # run_method <method> <extra args...>
    local method="$1"; shift
    local -a cmd=("$BIN" --method "$method" --repeats "$NREPS" "$@")
    if [ -n "$PIN" ]; then
        cmd=(taskset -c "$PIN" "${cmd[@]}")
    fi
    # The solver prints progress to stdout; the data row is the last line.
    "${cmd[@]}" 2>/dev/null | tail -n 1
}

printf '%-14s %9s %9s %9s %9s %8s %8s %8s %8s %12s %8s\n' \
    "method" "nodes" "elements" "store_MB" "expand_MB" "mesh_s" "view_s" "solve_s" \
    "band" "sumA2" "rss_MB"

if [ "$DIVISOR" -gt 1 ]; then
    echo "# tiled fixture reduced by instance divisor $DIVISOR (60-degree sector)" >&2
else
    echo "# full 360-degree tiled fixture" >&2
fi

run_method materialized --tiled "$TILED" --instance-divisor "$DIVISOR"
run_method native --tiled "$TILED" --instance-divisor "$DIVISOR"

if [ "$WITH_CONVENTIONAL" -eq 1 ]; then
    run_method conventional --redraw "$REDRAW"
fi

#!/bin/bash
# Benchmark the legacy vs PETSc linear solver backends on the radial-flux PM
# machine models.
#
# Usage: benchmark_radial_machine.sh [model.fem ...] [--iterations N] [--pin CORE]
# Defaults to the sliding model plus all redraw fixtures, 3 repetitions each.
# Both backends are pinned to the same CPU core (--pin) for a fair comparison
# on loaded machines.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$(dirname "$HERE")")"
BIN="$ROOT/cfemm/bin"
DATA="$ROOT/mfemm/testing/radial_machine/data"
WORK="${TMPDIR:-/tmp}/rmbench_work"
NREPS=3
PIN=""

MODELS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --iterations) NREPS="$2"; shift 2 ;;
        --pin) PIN="$2"; shift 2 ;;
        --) shift; MODELS+=("$@"); break ;;
        *) MODELS+=("$1"); shift ;;
    esac
done

if [ ${#MODELS[@]} -eq 0 ]; then
    MODELS=("$DATA/radial_machine_sliding.fem")
    for i in $(seq -w 1 10); do MODELS+=("$DATA/radial_machine_redraw_$i.fem"); done
fi

# prefix the solver with taskset when a core was requested
PIN_CMD=()
if [ -n "$PIN" ]; then
    PIN_CMD=(taskset -c "$PIN")
fi

mkdir -p "$WORK"

measure() { # measure <backend> <base> <meshdir>
    local backend="$1" base="$2" meshdir="$3"
    local total=0.0
    for rep in $(seq 1 "$NREPS"); do
        cp "$meshdir"/"$base".node "$meshdir"/"$base".ele "$meshdir"/"$base".pbc "$meshdir"/"$base".edge "$WORK"/ 2>/dev/null
        cd "$WORK"
        local start end
        start=$(date +%s.%N)
        XFEMM_SOLVER_BACKEND="$backend" "${PIN_CMD[@]}" "$BIN/fsolver" "$base" > "out_${backend}_${base}.log" 2>&1
        local rc=$?
        end=$(date +%s.%N)
        if [ $rc -ne 0 ]; then
            echo "  ERROR: $backend $base exited $rc"
            return 1
        fi
        total=$(python3 -c "print($total + ($end - $start))")
        cp "$WORK/$base.ans" "$WORK/$base.$backend.$rep.ans" 2>/dev/null
    done
    local avg
    avg=$(python3 -c "print($total / $NREPS)")
    local iters
    iters=$(grep -c "Newton Iteration" "$WORK/out_${backend}_${base}.log" 2>/dev/null || echo 0)
    echo "  $backend: avg $avg s over $NREPS runs, $iters nonlinear iterations"
}

printf "%-32s %-6s %s\n" "model" "nodes" "result"
for fem in "${MODELS[@]}"; do
    # resolve the model to an absolute path (we change directory below)
    case "$fem" in
        /*) : ;;
        *) fem="$(pwd)/$fem" ;;
    esac
    base="$(basename "${fem%.fem}")"
    cd "$WORK"
    cp "$fem" .
    "$BIN/fmesher" "$(basename "$fem")" >/dev/null 2>&1
    mkdir -p "mesh_$base"
    cp "$base".node "$base".ele "$base".pbc "$base".edge "mesh_$base"/ 2>/dev/null
    nodes=$(head -1 "$base".node | awk '{print $1}')
    printf "\n== %s (nodes=%s) ==\n" "$base" "$nodes"
    measure legacy "$base" "mesh_$base"
    measure petsc "$base" "mesh_$base"
done

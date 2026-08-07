#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
xfemm_root=$(CDPATH= cd -- "$script_dir/../../.." && pwd)
reference_root=${REFERENCE_XFEMM_ROOT:-}
reference_tag=${REFERENCE_XFEMM_TAG:-v4.0}
n_positions=${N_POSITIONS:-10}

if [[ ! "$n_positions" =~ ^[0-9]+$ ]] || (( n_positions < 1 || n_positions > 10 )); then
    echo "N_POSITIONS must be an integer from 1 to 10." >&2
    exit 2
fi

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/xfemm-radial-machine.XXXXXXXX")
cleanup() {
    if [[ ${KEEP_TEST_OUTPUTS:-0} == 1 ]]; then
        echo "Keeping test outputs in $work_dir"
    elif [[ -n "$work_dir" && -d "$work_dir" ]]; then
        rm -rf -- "$work_dir"
    fi
}
trap cleanup EXIT

if [[ -z "$reference_root" ]]; then
    reference_root="$work_dir/xfemm-$reference_tag"
    mkdir -p "$reference_root"
    if git -C "$xfemm_root" rev-parse --verify --quiet "refs/tags/$reference_tag" >/dev/null; then
        git -C "$xfemm_root" archive "$reference_tag" | tar -xf - -C "$reference_root"
    else
        rmdir "$reference_root"
        git clone --depth 1 --branch "$reference_tag" \
            https://github.com/crobarcro/xfemm.git "$reference_root"
    fi
fi

if [[ ! -f "$reference_root/mfemm/mfemm_setup.m" ]]; then
    echo "Reference xfemm root is invalid: $reference_root" >&2
    exit 2
fi

current_output="$work_dir/current.mat"
reference_output="$work_dir/reference.mat"

run_case() {
    local case_root=$1
    local output_file=$2
    XFEMM_CASE_ROOT="$case_root" \
    XFEMM_TEST_SUPPORT="$script_dir" \
    XFEMM_RESULT_FILE="$output_file" \
    XFEMM_N_POSITIONS="$n_positions" \
    octave --no-gui --quiet --no-init-file --eval \
      "xr=getenv('XFEMM_CASE_ROOT'); addpath(fullfile(xr,'mfemm')); mfemm_setup('RunTests',false); addpath(getenv('XFEMM_TEST_SUPPORT')); n=str2double(getenv('XFEMM_N_POSITIONS')); inds=unique(round(linspace(1,10,n))); radial_machine_fixture_case('redraw','PositionIndices',inds,'OutputFile',getenv('XFEMM_RESULT_FILE'));"
}

# Configuring generates cfemm/libfemm/femmversion.h in release source trees.
cmake -S "$reference_root/cfemm" -B "$work_dir/reference-build" \
      -DCMAKE_BUILD_TYPE=Release >/dev/null

echo "Running current xfemm redraw case..."
run_case "$xfemm_root" "$current_output"
echo "Running xfemm $reference_tag redraw case..."
run_case "$reference_root" "$reference_output"

XFEMM_TEST_SUPPORT="$script_dir" \
XFEMM_CURRENT_RESULT="$current_output" \
XFEMM_REFERENCE_RESULT="$reference_output" \
octave --no-gui --quiet --no-init-file --eval \
  "addpath(getenv('XFEMM_TEST_SUPPORT')); compare_radial_machine_result_files(getenv('XFEMM_REFERENCE_RESULT'),getenv('XFEMM_CURRENT_RESULT'));"

echo "Current xfemm versus $reference_tag static machine test passed."

#!/usr/bin/env bash
# Checks typed command-line params against the schema check; starts no processes.
# Usage: NEXTFLOW=/path/to/nextflow test/test_cli_params.sh

set -uo pipefail

NEXTFLOW=${NEXTFLOW:-nextflow}
REPO=$(cd "$(dirname "$0")/.." && pwd)
work_root=$(mktemp -d "${TMPDIR:-/tmp}/vimop_cli_test.XXXXXX") || { echo "cannot create a temporary directory"; exit 1; }
trap 'rm -rf "${work_root:?}"' EXIT

failures=0

run_case() {
    local name="$1" expected="$2"
    shift 2
    local launch_dir="$work_root/$name"
    mkdir -p "$launch_dir"
    (cd "$launch_dir" && "$NEXTFLOW" run "$REPO/main.nf" "$@" > out.txt 2>&1)
    local out="$launch_dir/out.txt"

    local actual
    if grep -q "Validation of pipeline parameters failed" "$out"; then
        actual=reject
    elif grep -q "Typical pipeline command" "$out"; then
        actual=help
    elif grep -q "No fastqs provided" "$out"; then
        actual=pass
    else
        actual=other
    fi

    if [[ "$actual" == "$expected" ]]; then
        echo "ok      $name ($*)"
    else
        echo "FAILED  $name ($*): expected $expected, got $actual"
        grep -E "expected type|ERROR|Error" "$out" | head -3 | sed 's/^/          /'
        failures=$((failures + 1))
    fi
}

run_case should_pass_number          pass   --min_cpus 1
run_case should_pass_integer         pass   --trim_length 30
run_case should_pass_fraction        pass   --centrifuge_plot_min_frac 0.001
run_case should_pass_bare_flag       pass   --output_krona_plot
run_case should_pass_false_flag      pass   --output_krona_plot false
run_case should_pass_string          pass   --consensus_method simple
run_case should_reject_text_number   reject --min_cpus abc
run_case should_reject_text_flag     reject --output_krona_plot maybe
run_case should_reject_fraction_int  reject --trim_length 2.5
run_case should_show_help            help   --help

if [[ $failures -gt 0 ]]; then
    echo "$failures case(s) failed"
    exit 1
fi
echo "all cases passed"

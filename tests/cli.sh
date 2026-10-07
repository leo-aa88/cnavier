#!/bin/sh
# Command-line checks: every invalid invocation must be refused the way the
# program refuses input (exit status 1 with an error or usage message) and
# write nothing; a crash or a missing binary is a failure, not a refusal. The
# valid edge cases must run.
# Usage: tests/cli.sh path/to/cnavier
set -u
bin=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
[ -x "$bin" ] || { echo "no executable at $bin"; exit 1; }
failed=0
count=0

# expect <0|nonzero> <arguments...>
expect() {
    want=$1
    shift
    count=$((count + 1))
    dir=$(mktemp -d)
    mkdir "$dir/output"
    out=$(cd "$dir" && "$bin" "$@" 2>&1)
    status=$?
    files=$(ls "$dir/output" | wc -l)
    rm -rf "$dir"
    if [ "$want" = 0 ]; then
        ok=$([ $status -eq 0 ] && echo yes || echo no)
    else
        refused=$(echo "$out" | grep -cE '\*\* Error|^Usage:')
        ok=$([ $status -eq 1 ] && [ "$refused" -gt 0 ] && [ "$files" -eq 0 ] && echo yes || echo no)
    fi
    if [ "$ok" = yes ]; then
        echo "  [ ok ] exit $status, $files files: cnavier $*"
    else
        echo "  [FAIL] exit $status, $files files: cnavier $*"
        echo "$out" | tail -3 | sed 's/^/         /'
        failed=$((failed + 1))
    fi
}

echo "Command line: invalid input is refused before anything is written"
expect nonzero --dt nan
expect nonzero --dt inf
expect nonzero --dt -1
expect nonzero --dt 0
expect nonzero --tf 1e9 --dt 1e-9
expect nonzero --tf 0.001
expect nonzero --n 64abc
expect nonzero --n ""
expect nonzero --n 7
expect nonzero --nx 16385
expect nonzero --ny 4
expect nonzero --n 16384
expect nonzero --output-interval -1
expect nonzero --output-interval 2x
expect nonzero --dt 0.02
expect nonzero --dt 0.0059
expect nonzero --re 0
expect nonzero --re -100
expect nonzero --re nan
expect nonzero --case nonsense
expect nonzero --case ""
expect nonzero --case taylor-green --dt 0.02
expect nonzero --bogus
expect nonzero stray-argument

echo "Command line: valid edge cases run"
expect 0 --help
expect 0 --tf 0.005 --output-interval 0
expect 0 --nx 33 --ny 17 --tf 0.01 --output-interval 0
expect 0 --case taylor-green --nx 32 --ny 24 --tf 0.01 --output-interval 0
expect 0 --case shear-layer --n 33 --re 1000 --tf 0.01 --output-interval 0

echo
echo "$count checks, $failed failed"
[ $failed -eq 0 ]

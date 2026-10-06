#!/bin/sh
# End-to-end checks of the default lid-driven cavity case (Re=100, 64x64).
# Usage: tests/regression.sh /absolute/path/to/cnavier
set -eu
bin=$1
here=$(cd "$(dirname "$0")" && pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
mkdir "$dir/output"
cd "$dir"

# 100 steps, compared with results stored in tests/reference. They are
# deterministic, so they move only if the numerics change: then update the
# stored files deliberately (see the README).
echo "Regression: 100 steps of the default case against stored results"
"$bin" --tf 0.5 --output-interval 0 > run_short.txt
python3 "$here/check_centerline.py" reference output "$here/reference" 1e-5

# The full default run (t = 30, steady) against the reference data. The
# current code is within 0.0024 (u) and 0.0074 (v).
echo "Regression: steady default case against Ghia et al. (1982)"
"$bin" --output-interval 0 > run_full.txt
python3 "$here/check_centerline.py" ghia output 0.004 0.010

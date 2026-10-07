#!/bin/sh
# End-to-end checks of the default lid-driven cavity case (Re=100, 64x64).
# Usage: tests/regression.sh path/to/cnavier
set -eu
bin=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
here=$(cd "$(dirname "$0")" && pwd)
dir=$(mktemp -d)
trap 'rm -rf "$dir"' EXIT
mkdir "$dir/output" "$dir/nan"
cd "$dir"

# The checker must fail on values that are not finite numbers: a solution
# that blew up would otherwise compare as a perfect match.
echo "Regression: the checker rejects a NaN solution"
printf 'y,u\n0,0\n0.5,nan\n1,1\n' > nan/centerline_u_sim.csv
printf 'x,v\n0,0\n0.5,nan\n1,0\n' > nan/centerline_v_sim.csv
if python3 "$here/check_centerline.py" reference nan "$here/reference" 1e-5 > /dev/null; then
    echo "  [FAIL] check_centerline.py accepted NaN values"
    exit 1
fi
echo "  [ ok ] check_centerline.py rejects NaN values"

# 100 steps, compared with results stored in tests/reference. They are
# deterministic, so they move only if the numerics change: then update the
# stored files deliberately (see the README).
echo "Regression: 100 steps of the default case against stored results"
"$bin" --tf 0.5 --output-interval 0 > run_short.txt
python3 "$here/check_centerline.py" reference output "$here/reference" 1e-5

# The full default run (t = 30, steady) against the reference data of Ghia
# et al. (1982), stored in tests/reference. The current code is within
# 0.0023 (u) and 0.0073 (v).
echo "Regression: steady default case against Ghia et al. (1982)"
"$bin" --output-interval 0 > run_full.txt
python3 "$here/check_centerline.py" ghia output "$here/reference" 0.004 0.010

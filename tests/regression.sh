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

# The periodic Taylor-Green vortex decays like the exact solution: on 32^2
# and 64^2 (sixth-order operators) the error at t = 0.5 is about 4e-8 and
# 6e-10 of max |w|, a factor of 64 = 2^6 apart.
echo "Regression: periodic Taylor-Green vortex against the exact solution"
"$bin" --case taylor-green --n 32 --dt 0.001 --tf 0.5 --output-interval 0 > tg32.txt
"$bin" --case taylor-green --n 64 --dt 0.001 --tf 0.5 --output-interval 0 > tg64.txt
e32=$(sed -n 's/.*relative to max |w| \([^ ]*\)$/\1/p' tg32.txt)
e64=$(sed -n 's/.*relative to max |w| \([^ ]*\)$/\1/p' tg64.txt)
if awk -v a="$e32" -v b="$e64" 'BEGIN { ok = (a + 0 > 0 && a < 1e-7 && b > 0 && a / b > 48); exit !ok }'; then
    echo "  [ ok ] relative error $e32 (32^2), $e64 (64^2): order $(awk -v a="$e32" -v b="$e64" 'BEGIN { printf "%.2f", log(a / b) / log(2) }')"
else
    echo "  [FAIL] relative error '$e32' (32^2), '$e64' (64^2): need < 1e-7 at 32^2 and a ratio above 48"
    exit 1
fi

# The energy the Taylor-Green run writes to integrals.csv decays as
# E = 1/4 exp(-4 nu k^2 t), k = 2 pi, nu = 1/Re = 0.01
echo "Regression: Taylor-Green energy in output/integrals.csv"
last=$(tail -1 output/integrals.csv)
if echo "$last" | awk -F, '{ e = 0.25 * exp(-4 * 0.01 * (2 * 3.14159265358979) ^ 2 * $2); r = ($3 - e) / e; if (r < 0) r = -r; exit !(NF == 5 && r < 1e-6) }'; then
    echo "  [ ok ] E = $(echo "$last" | cut -d, -f3) at t = $(echo "$last" | cut -d, -f2), within 1e-6 of the exact decay"
else
    echo "  [FAIL] last line of integrals.csv: '$last'"
    exit 1
fi

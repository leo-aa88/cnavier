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
if echo "$last" | awk -F, '{ e = 0.25 * exp(-4 * 0.01 * (2 * 3.14159265358979) ^ 2 * $2); r = ($3 - e) / e; if (r < 0) r = -r; exit !(NF == 7 && r < 1e-6) }'; then
    echo "  [ ok ] E = $(echo "$last" | cut -d, -f3) at t = $(echo "$last" | cut -d, -f2), within 1e-6 of the exact decay"
else
    echo "  [FAIL] last line of integrals.csv: '$last'"
    exit 1
fi

# Kolmogorov flow at Re = 50 settles to the laminar state: with u_L =
# A / (nu k^2), A = 1, k = 8 pi, E = u_L^2 / 4 and the energy input I = u_L / 2
echo "Regression: laminar Kolmogorov flow (energy and energy input)"
"$bin" --case kolmogorov --n 64 --re 50 --dt 0.002 --tf 2 --output-interval 0 > kolmogorov.txt
last=$(tail -1 output/integrals.csv)
if echo "$last" | awk -F, '{ ul = 1 / (0.02 * (8 * 3.14159265358979) ^ 2); e = ul * ul / 4; i = ul / 2;
                             re = ($3 - e) / e; ri = ($6 - i) / i; if (re < 0) re = -re; if (ri < 0) ri = -ri;
                             exit !(NF == 7 && re < 1e-4 && ri < 1e-4) }'; then
    echo "  [ ok ] E = $(echo "$last" | cut -d, -f3), I = $(echo "$last" | cut -d, -f6): the laminar values"
else
    echo "  [FAIL] last line of integrals.csv: '$last'"
    exit 1
fi

# Forced 2-D turbulence has a dual cascade: forced at |k|/2 pi = 8 on 64^2,
# with hyperviscosity and hypodrag confining the dissipation to the ends of
# the spectrum, the time-mean energy flux is negative (to large scales) below
# the forcing and the enstrophy flux positive (to small scales) above it, and
# not merely positive: at least 50, a fifth of the enstrophy input (seeds 1-3
# give 75-85). With the skew-symmetric nonlinear term the net enstrophy
# transfer is round-off.
echo "Regression: dual cascade of forced 2-D turbulence (flux signs)"
rm -f output/spectrum-1-*.csv
"$bin" --case forced --n 64 --re 1e7 --drag 0 --hypodrag 400 --hyperviscosity 1.6e-16 --hyper-order 4 \
    --forcing-rate 0.1 --forcing-k 8 --dt 2e-3 --tf 8 --output-interval 0 --spectrum-interval 50 \
    --integrals-interval 0 --advection skew > cascade.txt
python3 "$here/../tools/cascade.py" check output --kf 8 --from 3 --below 1,6 --above 10,20 --conserves-enstrophy \
    --min-pi-z 50

# The flow past a cylinder at Re = 40 on a coarse polar grid (--geometry
# polar; 96 x 57, r_inf = 20, about 20 s). It is deterministic, so its final
# C_D, L_w/D, theta_s and the area means E, Z, P (the last line of
# output/cylinder_history.csv) are compared with the stored values in
# tests/reference/cylinder_regression.csv to 1e-4 (relative): the run stops at
# max |dw/dt| < 1e-6, so this gate is tight against the numerics yet not
# against the step at which the run happens to stop. Regenerate the file
# deliberately if the numerics change (the README says how). Separately, the
# run must be steady, symmetric and divergence-free to round-off, and close to
# the reference solution of Gautier et al. (2013): a sanity check only, the
# coarse grid and the small domain being 0.2% (C_D), 0.3% (L_w/D) and 0.06
# deg away.
echo "Regression: cylinder at Re = 40 (polar grid) against stored results and Gautier et al. (2013)"
"$bin" --geometry polar --ntheta 96 --r-inf 20 --diag-interval 500 > cylinder.txt
python3 - "$here/reference/cylinder_regression.csv" output/cylinder_history.csv <<'PY'
import csv, math, sys
ref = {r["quantity"]: float(r["value"]) for r in csv.DictReader(open(sys.argv[1]))}
last = list(csv.DictReader(open(sys.argv[2])))[-1]
cols = {"cd": "cd", "wake_length": "wake_length", "separation_angle": "separation_angle", "E": "E", "Z": "Z", "P": "P"}
worst, bad = 0.0, []
for q, c in cols.items():
    if c not in last:
        bad.append("%s missing from the history" % c)
        continue
    v = float(last[c])
    e = abs(v - ref[q]) / abs(ref[q]) if math.isfinite(v) else math.inf
    worst = max(worst, e)
    if not e <= 1e-4:
        bad.append("%s %.10g (stored %.10g)" % (q, v, ref[q]))
if bad:
    print("  [FAIL] against tests/reference/cylinder_regression.csv (1e-4): " + "; ".join(bad))
    sys.exit(1)
print("  [ ok ] C_D, L_w/D, theta_s, E, Z, P vs stored results: max relative diff %.1e (limit 1e-4)" % worst)
PY
cd_=$(sed -n 's/^  C_D *\([^ ]*\) .*/\1/p' cylinder.txt)
lw=$(sed -n 's/^  L_w\/D *\([^ ]*\) .*/\1/p' cylinder.txt)
ts=$(sed -n 's/^  theta_s *\([^ ]*\) .*/\1/p' cylinder.txt)
asym=$(sed -n 's/.*asymmetry max |w(theta) + w(-theta)| = \([^ ]*\)$/\1/p' cylinder.txt)
cont=$(sed -n 's/^Continuity max: \([^ ]*\) | min: \([^ ]*\) .*/\1 \2/p' cylinder.txt)
if grep -q "^Steady at" cylinder.txt && awk -v c="$cd_" -v l="$lw" -v t="$ts" -v a="$asym" -v k="$cont" 'BEGIN {
    split(k, m, " "); e = m[1] + 0; if (-m[2] > e) e = -m[2]
    ok = (c / 1.4931 - 1) ^ 2 < 0.01 ^ 2 && (l / 2.2360 - 1) ^ 2 < 0.01 ^ 2 && (t - 126.3945) ^ 2 < 0.2 ^ 2 &&
         a + 0 < 1e-12 && e < 1e-12
    exit !ok }'; then
    echo "  [ ok ] steady, asymmetry $asym, continuity $cont; sanity vs Gautier et al.: C_D $cd_, L_w/D $lw, theta_s $ts (within 1%, 1%, 0.2 deg)"
else
    echo "  [FAIL] C_D '$cd_', L_w/D '$lw', theta_s '$ts', asymmetry '$asym', continuity '$cont'"
    tail -12 cylinder.txt
    exit 1
fi

# The steady-state check measures the rate of change over the steps actually
# taken since the previous diagnostic. Here the last interval is 1 step of a
# --diag-interval 250 (251 steps), far from steady (max |dw/dt| about 0.27 at
# t = 10.8): the run must stop unsteady, not count that step as 250 and pass
# --steady-tol 2e-3.
echo "Regression: cylinder steady-state check on a short last interval"
"$bin" --geometry polar --ntheta 32 --r-inf 5 --tf 10.83 --diag-interval 250 --steady-tol 2e-3 > unsteady.txt
if grep -q "^Stopped at t = 10.8181 after 251 steps" unsteady.txt; then
    echo "  [ ok ] $(grep '^Stopped' unsteady.txt)"
else
    echo "  [FAIL] expected 'Stopped at t = 10.8181 after 251 steps', got: $(grep -E '^(Steady|Stopped)' unsteady.txt)"
    exit 1
fi

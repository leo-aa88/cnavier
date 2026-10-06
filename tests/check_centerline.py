#!/usr/bin/env python3
"""Compare centerline CSVs written by cnavier.

check_centerline.py reference <output dir> <reference dir> <tolerance>
    The simulated profiles must match stored ones to within <tolerance>.

check_centerline.py ghia <output dir> <reference dir> <u tolerance> <v tolerance>
    Interpolated to the interior points of Ghia et al. (1982), stored in
    <reference dir>/ghia_{u,v}.csv, the simulated profiles must be within the
    tolerances of the reference values.

Any value that is not a finite number, or a missing or malformed file, is a
failure.
"""
import csv
import math
import sys


def load(path):
    with open(path) as f:
        rows = list(csv.reader(f))[1:]
    xs, ys = [float(r[0]) for r in rows], [float(r[1]) for r in rows]
    if not rows or not all(math.isfinite(v) for v in xs + ys):
        raise ValueError(f"{path}: empty, or contains a value that is not a finite number")
    return xs, ys


def interp(xs, ys, x):
    for k in range(len(xs) - 1):
        if xs[k] <= x <= xs[k + 1]:
            t = (x - xs[k]) / (xs[k + 1] - xs[k])
            return (1 - t) * ys[k] + t * ys[k + 1]
    raise ValueError(f"{x} is outside the sampled range")


def compare(mode, out, ref, tol, c):
    sx, sy = load(f"{out}/centerline_{c}_sim.csv")
    if mode == "reference":
        rx, ry = load(f"{ref}/centerline_{c}_short.csv")
        if len(rx) != len(sx):
            raise ValueError(f"{len(sx)} points, the reference has {len(rx)}")
        diffs = [abs(a - b) for a, b in zip(sx + sy, rx + ry)]
    else:
        gx, gy = load(f"{ref}/ghia_{c}.csv")
        # the end points lie on the walls, where the values are imposed
        diffs = [abs(interp(sx, sy, x) - y) for x, y in zip(gx[1:-1], gy[1:-1])]
    return max(diffs)


def main():
    mode, out, ref = sys.argv[1], sys.argv[2], sys.argv[3]
    what = "vs stored reference" if mode == "reference" else "vs Ghia et al. (1982)"
    failed = 0
    for i, c in enumerate(("u", "v")):
        tol = float(sys.argv[4] if mode == "reference" else sys.argv[4 + i])
        try:
            err = compare(mode, out, ref, tol, c)
            ok = math.isfinite(err) and err <= tol
            msg = f"max |diff| = {err:.2e} (limit {tol:.0e})"
        except (OSError, ValueError, IndexError) as e:
            ok, msg = False, str(e)
        failed += not ok
        print(f"  [{' ok ' if ok else 'FAIL'}] centerline {c} {what}: {msg}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

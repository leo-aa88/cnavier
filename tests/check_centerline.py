#!/usr/bin/env python3
"""Compare centerline CSVs written by cnavier.

check_centerline.py reference <output dir> <reference dir> <tolerance>
    The simulated profiles must match stored ones to within <tolerance>.

check_centerline.py ghia <output dir> <u tolerance> <v tolerance>
    Interpolated to the interior points of Ghia et al. (1982), the simulated
    profiles must be within the tolerances of the reference values.
"""
import csv
import sys


def load(path):
    with open(path) as f:
        rows = list(csv.reader(f))[1:]
    return [float(r[0]) for r in rows], [float(r[1]) for r in rows]


def interp(xs, ys, x):
    for k in range(len(xs) - 1):
        if xs[k] <= x <= xs[k + 1]:
            t = (x - xs[k]) / (xs[k + 1] - xs[k])
            return (1 - t) * ys[k] + t * ys[k + 1]
    raise ValueError(f"{x} is outside the sampled range")


def main():
    mode, out = sys.argv[1], sys.argv[2]
    failed = 0
    for c in ("u", "v"):
        sx, sy = load(f"{out}/centerline_{c}_sim.csv")
        if mode == "reference":
            rx, ry = load(f"{sys.argv[3]}/centerline_{c}_short.csv")
            tol = float(sys.argv[4])
            if len(rx) != len(sx):
                err = float("inf")
            else:
                err = max(max(abs(a - b) for a, b in zip(sx, rx)), max(abs(a - b) for a, b in zip(sy, ry)))
            what = "vs stored reference"
        else:
            gx, gy = load(f"{out}/centerline_{c}_ghia.csv")
            tol = float(sys.argv[3] if c == "u" else sys.argv[4])
            # the end points lie on the walls, where the values are imposed
            err = max(abs(interp(sx, sy, x) - y) for x, y in zip(gx[1:-1], gy[1:-1]))
            what = "vs Ghia et al. (1982)"
        ok = err <= tol
        failed += not ok
        print(f"  [{' ok ' if ok else 'FAIL'}] centerline {c} {what}: max |diff| = {err:.2e} (limit {tol:.0e})")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

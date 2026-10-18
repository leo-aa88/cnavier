#!/usr/bin/env python3
"""Validate the polar solver on the flow past a cylinder at Re = 40 against
the reference solution of Gautier, Biau and Lamballais (2013, arXiv:1310.6641;
tests/reference/cylinder_re40.csv).

    tools/cylinder_validate.py run BIN OUT_DIR [--grids 128,256,512]
                               [--r-inf 30,40,50] [-- EXTRA_ARGS...]
        runs BIN --geometry polar to steady state for every grid (ntheta,
        square cells) and outer radius, each in OUT_DIR/n<N>_r<R>/ with its
        log.txt and output/; EXTRA_ARGS go to every run (e.g.
        --outer-bc freestream)
    tools/cylinder_validate.py table RUN_DIR [...]
        the summary of each run against the reference: C_D, L_w/D, theta_s,
        the recirculation centre and the velocity at (3, 3), with the
        observed order between successive grids of one outer radius and the
        Richardson extrapolation of C_D and L_w/D from the two finest
    tools/cylinder_validate.py plot OUT.png RUN_DIR [...]
        the centreline velocity u_x on the wake axis (Gautier et al., Fig. 1:
        it must rise monotonically to 1) and the wall vorticity of each run
"""
import csv
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF_FILE = os.path.join(HERE, "..", "tests", "reference", "cylinder_re40.csv")


def reference():
    with open(REF_FILE) as f:
        return {row["quantity"]: float(row["value"]) for row in csv.DictReader(f)}


def run(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    binary, out = os.path.abspath(argv[0]), argv[1]
    grids, radii, extra = [128, 256], [40.0], []
    rest = argv[2:]
    if "--" in rest:
        extra = rest[rest.index("--") + 1:]
        rest = rest[: rest.index("--")]
    for k in range(0, len(rest), 2):
        if rest[k] == "--grids":
            grids = [int(x) for x in rest[k + 1].split(",")]
        elif rest[k] == "--r-inf":
            radii = [float(x) for x in rest[k + 1].split(",")]
        else:
            sys.exit(__doc__)
    for r in radii:
        for n in grids:
            d = os.path.join(out, "n%d_r%g" % (n, r))
            os.makedirs(os.path.join(d, "output"), exist_ok=True)
            args = [binary, "--geometry", "polar", "--ntheta", str(n), "--r-inf", str(r), "--tf", "600",
                    "--diag-interval", str(max(200, n * 4))] + extra
            print("running", " ".join(args[1:]), "in", d, flush=True)
            with open(os.path.join(d, "log.txt"), "w") as log:
                subprocess.run(args, cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)


NUM = r"([-+0-9.eE]+)"


def summary(run_dir):
    """The final summary of a run's log.txt"""
    text = open(os.path.join(run_dir, "log.txt")).read()

    def grab(pattern, n=1):
        m = re.search(pattern, text)
        if not m:
            return None if n == 1 else (None,) * n
        return float(m.group(1)) if n == 1 else tuple(float(m.group(k + 1)) for k in range(n))

    s = {"dir": run_dir}
    m = re.search(r"Grid: (\d+) \(theta\) x (\d+) \(xi\), r from [-+0-9.eE]+ to " + NUM, text)
    s["ntheta"], s["nr"], s["r_inf"] = int(m.group(1)), int(m.group(2)), float(m.group(3))
    s["outer"] = "free stream" if "outer boundary: free stream" in text else "far field"
    m = re.search(r"Poisson: order (\d)", text)
    s["poisson"] = int(m.group(1)) if m else 2
    s["cd"] = grab(r"C_D\s+" + NUM + r"\s+reference")
    s["lw"] = grab(r"L_w/D\s+" + NUM + r"\s+reference")
    s["ts"] = grab(r"theta_s\s+" + NUM + r"\s+reference")
    s["a"], s["b"] = grab(r"a/D, b/D\s+" + NUM + r"\s+" + NUM, 2)
    s["ux"], s["uy"] = grab(r"at \(3, 3\) u_x " + NUM + r" \(reference [^)]*\), u_y " + NUM, 2)
    s["steady"] = "Steady at" in text
    s["ms"] = grab(r"\| " + NUM + r" ms per step")
    s["t"] = grab(r"(?:Steady|Stopped) at t = " + NUM)
    return s


def table(dirs):
    ref = reference()
    runs = sorted((summary(d) for d in dirs), key=lambda s: (s["r_inf"], s["ntheta"]))
    print("| grid | r_inf | outer | C_D | err | L_w/D | err | theta_s | err (deg) | a/D | b/D | u_x(3,3) | u_y(3,3) | t | ms/step |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in runs:
        print("| %dx%d | %g | %s | %.5f | %+.2f%% | %.5f | %+.2f%% | %.4f | %+.3f | %.4f | %.4f | %.6f | %.6f | %.0f%s | %.2f |" % (
            s["ntheta"], s["nr"], s["r_inf"], s["outer"], s["cd"], 100 * (s["cd"] / ref["cd"] - 1), s["lw"],
            100 * (s["lw"] / ref["wake_length"] - 1), s["ts"], s["ts"] - ref["separation_angle"], s["a"], s["b"],
            s["ux"], s["uy"], s["t"], "" if s["steady"] else " (not steady)", s["ms"]))
    print("| reference | 40 | | %.4f | | %.4f | | %.4f | | %.2f | %.2f | %.6f | %.6f | | |" % (
        ref["cd"], ref["wake_length"], ref["separation_angle"], ref["vortex_a"], ref["vortex_b"], ref["ux_3_3"],
        ref["uy_3_3"]))
    # Observed order from three successive grids (refinement ratio 2) of one
    # outer radius, and Richardson extrapolation from the two finest
    # (one series per outer radius, outer boundary and Poisson order)
    for key in sorted({(s["r_inf"], s["outer"], s["poisson"]) for s in runs}):
        g = [s for s in runs if (s["r_inf"], s["outer"], s["poisson"]) == key]
        r = "%g (%s, Poisson order %d)" % key
        for q in ("cd", "lw", "ts"):
            if len(g) >= 3:
                d1, d2 = g[-2][q] - g[-3][q], g[-1][q] - g[-2][q]
                if d1 != 0 and d2 != 0 and d1 / d2 > 0:
                    p = math.log2(d1 / d2)
                    ext = g[-1][q] + d2 / (2 ** p - 1)
                    print("r_inf %s, %s: observed order %.2f, extrapolated %.5f" % (r, q, p, ext))
            elif len(g) == 2:
                print("r_inf %s, %s: change %.2e from %d to %d" % (r, q, g[1][q] - g[0][q], g[0]["ntheta"],
                                                                   g[1]["ntheta"]))


def plot(out, dirs):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for d in dirs:
        s = summary(d)
        label = "%dx%d, r_inf %g, %s" % (s["ntheta"], s["nr"], s["r_inf"], s["outer"])
        with open(os.path.join(d, "output", "cylinder_centreline.csv")) as f:
            rows = [(float(r["x"]), float(r["u_x"])) for r in csv.DictReader(f)]
        wake = [(x, u) for x, u in rows if x >= 0.5]
        ax[0].plot([x for x, _ in wake], [u for _, u in wake], label=label)
        with open(os.path.join(d, "output", "cylinder_surface.csv")) as f:
            rows = [(float(r["theta_deg"]), float(r["w_wall"])) for r in csv.DictReader(f)]
        ax[1].plot([180 - t for t, _ in rows if t <= 180], [w for t, w in rows if t <= 180], label=label)
    ax[0].axhline(1.0, color="0.6", lw=0.8)
    ax[0].set_xlabel("x / D (wake axis)")
    ax[0].set_ylabel("u_x / U")
    ax[0].set_title("Centreline velocity")
    ax[1].axhline(0.0, color="0.6", lw=0.8)
    ax[1].set_xlabel("angle from the front stagnation point (deg)")
    ax[1].set_ylabel("wall vorticity (upper surface)")
    ax[1].set_title("Wall vorticity")
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print("wrote", out)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    cmd, argv = sys.argv[1], sys.argv[2:]
    if cmd == "run":
        run(argv)
    elif cmd == "table":
        table(argv)
    elif cmd == "plot" and len(argv) >= 2:
        plot(argv[0], argv[1:])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compare the periodic schemes of one case against a reference run.

    tools/compare_schemes.py CASE_DIR REF_DIR --at T [--grids 256,512]
        CASE_DIR holds runs named n<N>_<scheme> (scheme: o6, compact6,
        spectral, spectral23), each with output/spectrum-1-*.csv and a
        log.txt with the solver's "ms per step" line; REF_DIR is the
        reference run. For each grid and scheme, at the frame nearest to T:
        the shells within 10 % and 2 % of the reference, the largest ratio
        E / E_ref (pile-up at the cutoff), the reference's E at the grid's
        Nyquist shell relative to its peak, the time per step, and the cost
        at equal resolved range relative to explicit order 6: at the
        advective limit the step goes as 1 / (k* h)max, and explicit order 6
        resolves a range proportional to n at a cost per unit time going as
        n^3, so matching a range K costs (K / K_o6)^3 times its own.
    tools/compare_schemes.py crossover OUT.png CASE:REF:T [CASE:REF:T ...]
        for every case and grid: the resolved range (to 10 %) of each scheme
        relative to explicit order 6, against the reference's energy at the
        grid's Nyquist shell relative to its peak; does the ordering of the
        schemes follow that one parameter?
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cascade  # noqa: E402

KMAX_H = {"o6": 1.586, "compact6": 1.989, "spectral": 3.14159265, "spectral23": 2.0943951}
SCHEMES = ["o6", "compact6", "spectral", "spectral23"]


def ms_per_step(run_dir):
    try:
        with open(os.path.join(run_dir, "log.txt")) as f:
            m = re.search(r"([0-9.]+) ms per step", f.read())
            return float(m.group(1)) if m else None
    except OSError:
        return None


def measure(ref, run_dir, t):
    _, rows = cascade.frame_at(os.path.join(run_dir, "output"), t)
    nyq = cascade.nyquist(rows)
    ratio = [rows[k]["E"] / ref[k]["E"] for k in range(1, nyq + 1)]
    first = {}
    for th in (0.02, 0.1):
        bad = next((k for k, q in enumerate(ratio, start=1) if abs(q - 1) > th), None)
        first[th] = (bad or nyq + 1) - 1
    large = max(abs(q - 1) for q in ratio[:20])
    return first[0.1], first[0.02], max(ratio), nyq, large


def crossover(png, specs):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    points = {s: [] for s in SCHEMES}
    print(f"{'case':<10} {'grid':>5} {'E(K_N)/E_max':>12}  K10 relative to order 6: compact6, spectral, spectral23")
    for spec in specs:
        case, ref_dir, t = spec.split(":")
        _, ref = cascade.frame_at(os.path.join(ref_dir, "output"), float(t))
        peak = max(r["E"] for r in ref[1:])
        for n in (256, 512):
            res = {}
            for s in SCHEMES:
                d = os.path.join(case, f"n{n}_{s}")
                if os.path.isdir(d):
                    res[s] = measure(ref, d, float(t))
            if "o6" not in res:
                continue
            occ = ref[res["o6"][3]]["E"] / peak
            rel = {s: res[s][0] / res["o6"][0] for s in res}
            for s in rel:
                points[s].append((occ, rel[s]))
            print(f"{os.path.basename(case):<10} {n:>5} {occ:12.1e}  " +
                  "  ".join(f"{rel.get(s, float('nan')):.2f}" for s in SCHEMES[1:]))
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    style = {"o6": "k.", "compact6": "C1o", "spectral": "C0s", "spectral23": "C3^"}
    for s in SCHEMES:
        if points[s]:
            x, y = zip(*sorted(points[s]))
            ax.semilogx(x, y, style[s], label=s)
    ax.axhline(1, color="gray", lw=0.8)
    ax.set_xlabel("reference E(K_Nyquist) / E_peak")
    ax.set_ylabel("shells within 10 %, relative to order 6")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=90)
    print(f"wrote {png}")


def main(argv):
    if len(argv) >= 4 and argv[1] == "crossover":
        crossover(argv[2], argv[3:])
        return 0
    if len(argv) < 5 or "--at" not in argv:
        print(__doc__)
        return 1
    case, ref_dir = argv[1], argv[2]
    t = float(argv[argv.index("--at") + 1])
    grids = [int(g) for g in argv[argv.index("--grids") + 1].split(",")] if "--grids" in argv else [256, 512]
    tr, ref = cascade.frame_at(os.path.join(ref_dir, "output"), t)
    peak = max(r["E"] for r in ref[1:])
    print(f"{case}: reference {ref_dir} at t = {tr:g}")
    print(f"{'grid':>5} {'scheme':<11} {'K10':>4} {'K2':>4} {'max E/Eref':>10} {'K<=20 dev':>9} {'ms':>6} "
          f"{'cost vs o6':>10}")
    for n in grids:
        res = {}
        for s in SCHEMES:
            d = os.path.join(case, f"n{n}_{s}")
            if os.path.isdir(d):
                res[s] = measure(ref, d, t) + (ms_per_step(d),)
        if "o6" not in res:
            continue
        k_o6, ms_o6 = res["o6"][0], res["o6"][5]
        nyq = res["o6"][3]
        print(f"{n:>5} reference E(K_Nyquist = {nyq}) / E_peak = {ref[nyq]['E'] / peak:.1e}")
        for s, (k10, k2, pile, _, large, ms) in res.items():
            cost = ""
            if ms and ms_o6 and k10 > 0:
                own = ms * KMAX_H[s]
                explicit = ms_o6 * KMAX_H["o6"] * (k10 / k_o6) ** 3
                cost = f"{own / explicit:10.2f}"
            print(f"{'':>5} {s:<11} {k10:4d} {k2:4d} {pile:10.2f} {large:9.1e} {ms or float('nan'):6.1f} {cost}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/env python3
"""Time averages, plots and checks of the spectra cnavier writes on periodic grids.

    tools/cascade.py average OUTPUT_DIR [--from T]
        time-average output/spectrum-1-*.csv over the frames with t >= T and
        write OUTPUT_DIR/spectrum-mean.csv (same columns, plus the number of
        frames in a comment)
    tools/cascade.py plot OUTPUT_DIR [--from T] [--kf K] [--png FILE] [--title S] [--compare DIR2 --label S2 --name S1]
        plot the averaged energy spectrum and the energy and enstrophy fluxes,
        up to the Nyquist shell (needs matplotlib); --compare overlays the
        fluxes of a second run (dashed), averaged over the same times, with
        legend labels S2 and S1
    tools/cascade.py check OUTPUT_DIR --kf K [--from T] [--below a,b] [--above c,d] [--conserves-enstrophy]
                     [--min-pi-e X] [--min-pi-z Y] [--spread-e S] [--spread-z S]
        print the mean budget, and check the dual cascade: the energy flux
        negative for K in [a, b] and the enstrophy flux positive for K in
        [c, d], with K = |k| / (2 pi / L); with --conserves-enstrophy (runs
        with --advection skew) also that the net nonlinear enstrophy transfer
        is round-off. Plateaus, optionally: -Pi_E >= X on every shell of
        [a, b], Pi_Z >= Y on every shell of [c, d], and (max - min) / mean of
        the flux over its range at most S. Exit status 1 if any check fails
    tools/cascade.py compare REF_DIR DIR [DIR ...] --at T [--names a,b,...] [--png FILE]
        the energy spectra of runs of the same flow at the frame nearest to
        time T, and their ratio to the reference run's, against K / K_Nyquist
        of each run; prints where each departs from the reference by 2 % and
        10 % (needs matplotlib for the plot)
    tools/cascade.py evolution OUTPUT_DIR --times a,b,... [--png FILE]
        E(t) and Z(t) from integrals.csv, relative to their initial values,
        and the energy spectra at the frames nearest to the given times
        (needs matplotlib)
    tools/cascade.py budget OUTPUT_DIR --from T --eps EPS --kf K [--dk DK]
        the energy budget of a randomly forced run over t >= T: the drift
        dE/dt from integrals.csv, the time-mean dissipation at small and large
        scales and the net nonlinear transfer from the spectra (with standard
        errors from block averages), the input they imply, and how far that
        may differ from EPS because the kicks inject EPS only on average: over
        a window T_w the realized input has a standard deviation of about
        sqrt(2 EPS E_f / (M T_w)), E_f the energy of the M forced modes

Wavenumbers are printed and compared in units of 2 pi / L (shell index).
Needs only the Python standard library, except `plot`.
"""
import csv
import glob
import math
import os
import re
import sys

COLUMNS = ["k", "E", "Z", "Pi_E", "Pi_Z", "D_E", "D_Z", "F_E", "F_Z"]


def read_spectrum(path):
    t = None
    rows = []
    with open(path) as f:
        first = f.readline()
        m = re.match(r"#\s*t\s*=\s*(\S+)", first)
        if m:
            t = float(m.group(1))
        reader = csv.reader(f)
        header = next(reader)
        for r in reader:
            vals = [float(x) for x in r]
            if any(not math.isfinite(v) for v in vals):
                raise ValueError(f"{path}: value that is not a finite number")
            rows.append(dict(zip(header, vals)))
    return t, rows


def frames(directory, t_from):
    paths = glob.glob(os.path.join(directory, "spectrum-1-*.csv"))
    paths.sort(key=lambda p: int(re.search(r"spectrum-1-(\d+)\.csv$", p).group(1)))
    out = []
    for p in paths:
        t, rows = read_spectrum(p)
        if t is not None and t >= t_from:
            out.append((t, rows))
    if not out:
        raise SystemExit(f"no spectrum-1-*.csv with t >= {t_from} in {directory}")
    return out


def average(directory, t_from):
    fr = frames(directory, t_from)
    n = len(fr[0][1])
    cols = [c for c in COLUMNS if c in fr[0][1][0]]
    mean = [{c: 0.0 for c in cols} for _ in range(n)]
    for _, rows in fr:
        for b in range(n):
            for c in cols:
                mean[b][c] += rows[b][c] / len(fr)
    dk = mean[1]["k"] if n > 1 else 1.0
    return mean, cols, dk, len(fr), fr[0][0], fr[-1][0]


def cmd_average(directory, t_from):
    mean, cols, dk, count, t0, t1 = average(directory, t_from)
    path = os.path.join(directory, "spectrum-mean.csv")
    with open(path, "w") as f:
        f.write(f"# mean over {count} frames, t = {t0:g} .. {t1:g}\n")
        f.write(",".join(cols) + "\n")
        for row in mean:
            f.write(",".join(f"{row[c]:.17g}" for c in cols) + "\n")
    print(f"wrote {path} ({count} frames, t = {t0:g} .. {t1:g})")


def nyquist(mean):
    """Shell index of the Nyquist wavenumber of a square grid: the bins run
    out to the corners of the spectrum, sqrt(2) times further"""
    return int((len(mean) - 1) / math.sqrt(2))


def budget(mean):
    """Time-mean sums: the dissipation at small scales (D) and at large scales
    (F), and the net transfer of the nonlinear term, -Pi through the last
    shell, which a conservative discretisation would make zero"""
    out = {}
    for q in ("E", "Z"):
        if f"D_{q}" in mean[0]:
            out[f"D_{q}"] = sum(r[f"D_{q}"] for r in mean)
        if f"F_{q}" in mean[0]:
            out[f"F_{q}"] = sum(r[f"F_{q}"] for r in mean)
        out[f"T_{q}"] = -mean[-1][f"Pi_{q}"]
    return out


def cmd_check(directory, t_from, kf, below, above, conserves, plateau):
    mean, _, dk, count, t0, t1 = average(directory, t_from)
    ok = True
    print(f"Mean over {count} frames, t = {t0:g} .. {t1:g}; K = |k| dk^-1, forcing at K = {kf:g}")
    b = budget(mean)
    for q, name in (("E", "energy"), ("Z", "enstrophy")):
        parts = [f"{label} {b[key]:.4g}" for key, label in ((f"D_{q}", "small scales"), (f"F_{q}", "large scales"))
                 if key in b]
        print(f"  {name}: removed at {', '.join(parts)}; net nonlinear transfer {b[f'T_{q}']:.3g}")
    peak_z = max(r["Pi_Z"] for r in mean)
    print(f"  peak fluxes: Pi_E {min(r['Pi_E'] for r in mean):.4g}, Pi_Z {peak_z:.4g}")
    if conserves:
        if abs(b["T_Z"]) <= 1e-10 * abs(peak_z):
            print(f"  [ ok ] net nonlinear enstrophy transfer {b['T_Z']:.1e}: round-off")
        else:
            print(f"  [FAIL] net nonlinear enstrophy transfer {b['T_Z']:.3e}, expected round-off")
            ok = False
    for b, row in enumerate(mean):
        if below[0] <= b <= below[1] and not row["Pi_E"] < 0.0:
            print(f"  [FAIL] Pi_E(K = {b}) = {row['Pi_E']:.3e}, expected < 0 (inverse energy cascade)")
            ok = False
        if above[0] <= b <= above[1] and not row["Pi_Z"] > 0.0:
            print(f"  [FAIL] Pi_Z(K = {b}) = {row['Pi_Z']:.3e}, expected > 0 (direct enstrophy cascade)")
            ok = False
    e = [mean[b]["Pi_E"] for b in range(below[0], below[1] + 1)]
    z = [mean[b]["Pi_Z"] for b in range(above[0], above[1] + 1)]
    for name, flux, sign, rng, low, spread in (("Pi_E", e, -1.0, below, plateau["--min-pi-e"], plateau["--spread-e"]),
                                               ("Pi_Z", z, 1.0, above, plateau["--min-pi-z"], plateau["--spread-z"])):
        mag = [sign * x for x in flux]
        if low is not None:
            good = min(mag) >= low
            print(f"  [{' ok ' if good else 'FAIL'}] |{name}| >= {low:g} for K = {rng[0]}..{rng[1]} "
                  f"(smallest {min(mag):.4g})")
            ok = ok and good
        if spread is not None:
            m = sum(mag) / len(mag)
            sp = (max(mag) - min(mag)) / abs(m) if m else math.inf
            good = sp <= spread
            print(f"  [{' ok ' if good else 'FAIL'}] {name} over K = {rng[0]}..{rng[1]}: mean {sign * m:.4g}, "
                  f"spread (max - min) / mean {sp:.3f} <= {spread:g}")
            ok = ok and good
    if ok:
        print(f"  [ ok ] Pi_E < 0 for K = {below[0]}..{below[1]} ({min(e):.3e} .. {max(e):.3e})")
        print(f"  [ ok ] Pi_Z > 0 for K = {above[0]}..{above[1]} ({min(z):.3e} .. {max(z):.3e})")
    return ok


def cmd_plot(directory, t_from, kf, png, title, compare, label, name):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mean, _, dk, count, t0, t1 = average(directory, t_from)
    K = [b for b in range(1, nyquist(mean) + 1)]
    E = [mean[b]["E"] / dk for b in K]
    PE = [mean[b]["Pi_E"] for b in K]
    PZ = [mean[b]["Pi_Z"] for b in K]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].loglog(K, E, "k-", label="E(k)")
    if kf:
        # Reference slopes, each through the spectrum at a point inside its
        # range (and a factor 3 above it, to keep the curve visible): the
        # middle of 1..kf for k^-5/3, of kf..K_max/2 for k^-3
        def ref(k0, slope, ks, style, label):
            i0 = min(range(len(K)), key=lambda i: abs(K[i] - k0))
            if ks:
                ax[0].loglog(ks, [3 * E[i0] * (k / K[i0]) ** slope for k in ks], style, label=label)

        ref(math.sqrt(kf), -5 / 3, [k for k in K if k <= kf], "b--", "k^-5/3")
        ref(math.sqrt(kf * K[-1] / 2), -3, [k for k in K if kf <= k <= K[-1] / 2], "r--", "k^-3")
        for a in ax:
            a.axvline(kf, color="gray", lw=0.8)
    ax[0].set_xlabel("K = |k| / (2 pi / L)")
    ax[0].set_title("energy spectrum (density)")
    ax[0].legend()
    if compare:
        other = average(compare, t_from)[0]
        Ko = [b for b in range(1, nyquist(other) + 1)]
        ax[1].semilogx(Ko, [other[b]["Pi_E"] for b in Ko], "k--", lw=1, label=label or compare)
        ax[2].semilogx(Ko, [other[b]["Pi_Z"] for b in Ko], "k--", lw=1, label=label or compare)
    ax[1].semilogx(K, PE, "b-", label=name or directory)
    ax[1].axhline(0, color="gray", lw=0.8)
    ax[1].set_title("energy flux Pi_E")
    ax[1].set_xlabel("K")
    ax[2].semilogx(K, PZ, "r-", label=name or directory)
    if compare:
        ax[1].legend()
        ax[2].legend()
    ax[2].axhline(0, color="gray", lw=0.8)
    ax[2].set_title("enstrophy flux Pi_Z")
    ax[2].set_xlabel("K")
    fig.suptitle(f"{title or directory}: mean over {count} frames, t = {t0:g} .. {t1:g}")
    fig.tight_layout()
    png = png or os.path.join(directory, "cascade.png")
    fig.savefig(png, dpi=90)
    print(f"wrote {png}")


def frame_at(directory, t):
    """The spectrum frame of directory nearest to time t"""
    best = None
    for p in glob.glob(os.path.join(directory, "spectrum-1-*.csv")):
        tp, rows = read_spectrum(p)
        if tp is not None and (best is None or abs(tp - t) < abs(best[0] - t)):
            best = (tp, rows)
    if best is None:
        raise SystemExit(f"no spectrum-1-*.csv in {directory}")
    return best


def cmd_compare(ref_dir, dirs, t, names, png):
    tr, ref = frame_at(ref_dir, t)
    names = names or dirs
    curves = []
    print(f"Reference {ref_dir} at t = {tr:g}")
    for d, name in zip(dirs, names):
        td, rows = frame_at(d, t)
        nyq = nyquist(rows)
        K = list(range(1, nyq + 1))
        ratio = [rows[k]["E"] / ref[k]["E"] for k in K]
        first = {}
        for th in (0.02, 0.1):
            first[th] = next((k for k, q in zip(K, ratio) if abs(q - 1) > th), None)
        peak = max(r["E"] for r in ref[1:])
        print(f"  {name} (t = {td:g}, K_Nyquist {nyq}): within 2 % up to K = "
              f"{(first[0.02] or nyq + 1) - 1}, within 10 % up to K = {(first[0.1] or nyq + 1) - 1}; "
              f"max E/E_ref {max(ratio):.2f}; reference E(K_Nyquist)/E_peak {ref[nyq]['E'] / peak:.1e}")
        curves.append((name, K, nyq, [rows[k]["E"] for k in K], ratio))
    if png is None:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    nyq_ref = nyquist(ref)
    ax[0].loglog(range(1, nyq_ref + 1), [ref[k]["E"] for k in range(1, nyq_ref + 1)], "k-", lw=2, label="reference")
    for name, K, nyq, E, ratio in curves:
        ax[0].loglog(K, E, label=name)
        ax[1].plot([k / nyq for k in K], ratio, label=name)
    ax[0].set_xlabel("K = |k| / (2 pi / L)")
    ax[0].set_title(f"energy spectrum (shell sums), t = {tr:g}")
    ax[0].legend(fontsize=8)
    ax[1].axhline(1, color="gray", lw=0.8)
    ax[1].set_ylim(0, 1.5)
    ax[1].set_xlabel("K / K_Nyquist")
    ax[1].set_title("ratio to the reference")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=90)
    print(f"wrote {png}")


def cmd_evolution(directory, times, png):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(os.path.join(directory, "integrals.csv")) as f:
        rows = [dict(zip(["step", "t", "E", "Z", "P", "I", "I_disc"], map(float, line.split(","))))
                for line in f.read().split("\n")[1:] if line]
    t = [r["t"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].semilogy(t, [r["E"] / rows[0]["E"] for r in rows], "b-", label="E / E(0)")
    ax[0].semilogy(t, [r["Z"] / rows[0]["Z"] for r in rows], "r-", label="Z / Z(0)")
    ax[0].set_xlabel("t")
    ax[0].set_title("energy and enstrophy")
    ax[0].legend()
    for tt in times:
        tf, spec = frame_at(directory, tt)
        nyq = nyquist(spec)
        dk = spec[1]["k"]
        ax[1].loglog(range(1, nyq + 1), [spec[k]["E"] / dk for k in range(1, nyq + 1)], label=f"t = {tf:.3g}")
    K = [k for k in range(10, nyquist(spec) // 2)]
    if K:
        ref = spec[K[0]]["E"] / dk * 3
        ax[1].loglog(K, [ref * (k / K[0]) ** -3 for k in K], "k--", lw=0.8, label="k^-3")
    ax[1].set_xlabel("K = |k| / (2 pi / L)")
    ax[1].set_title("energy spectrum (density)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    png = png or os.path.join(directory, "evolution.png")
    fig.savefig(png, dpi=90)
    print(f"wrote {png}")
    for r in rows[:: max(1, len(rows) // 10)] + [rows[-1]]:
        print(f"  t = {r['t']:8.3f}  E = {r['E']:.6g}  Z = {r['Z']:.6g}")


def block_error(x):
    """Mean and its standard error, from the variance of block means: blocks
    of increasing size until the estimate stops growing (correlated samples)"""
    n = len(x)
    m = sum(x) / n
    best = 0.0
    size = 1
    while n // size >= 8:
        blocks = [sum(x[i * size:(i + 1) * size]) / size for i in range(n // size)]
        nb = len(blocks)
        var = sum((b - m) ** 2 for b in blocks) / (nb - 1)
        best = max(best, math.sqrt(var / nb))
        size *= 2
    return m, best


def cmd_budget(directory, t_from, eps, kf, dk):
    with open(os.path.join(directory, "integrals.csv")) as f:
        rows = [list(map(float, line.split(","))) for line in f.read().split("\n")[1:] if line]
    rows = [r for r in rows if r[1] >= t_from]
    fr = frames(directory, t_from)
    t0, t1 = fr[0][0], fr[-1][0]
    e0 = min(rows, key=lambda r: abs(r[1] - t0))
    e1 = min(rows, key=lambda r: abs(r[1] - t1))
    drift = (e1[2] - e0[2]) / (e1[1] - e0[1])
    series = {q: [sum(r[q] for r in rows_) for _, rows_ in fr] for q in ("D_E", "F_E")}
    series["T_E"] = [-rows_[-1]["Pi_E"] for _, rows_ in fr]
    stats = {q: block_error(v) for q, v in series.items()}
    implied = drift + stats["D_E"][0] + stats["F_E"][0] - stats["T_E"][0]
    # The forced modes: the half plane of the square lattice within kf +- dk
    nyq = nyquist(fr[0][1])
    modes = sum(1 for n in range(0, nyq + 1) for m in range(-nyq + 1, nyq)
                if not (n == 0 and m <= 0) and abs(math.hypot(m, n) - kf) <= dk)
    mean = [sum(rows_[b]["E"] for _, rows_ in fr) / len(fr) for b in range(len(fr[0][1]))]
    e_forced = sum(mean[b] for b in range(len(mean)) if abs(b - kf) <= dk)
    sd_input = math.sqrt(2.0 * eps * e_forced / (modes * (t1 - t0)))
    se = math.sqrt(stats["D_E"][1] ** 2 + stats["F_E"][1] ** 2 + stats["T_E"][1] ** 2 + sd_input ** 2)
    print(f"Energy budget over t = {t0:g} .. {t1:g} ({len(fr)} spectrum frames)")
    print(f"  drift dE/dt                   {drift:+.4g}   (E = {e0[2]:.4g} -> {e1[2]:.4g})")
    print(f"  small-scale dissipation D_E   {stats['D_E'][0]:.4g} +- {stats['D_E'][1]:.2g}")
    print(f"  large-scale dissipation F_E   {stats['F_E'][0]:.4g} +- {stats['F_E'][1]:.2g}")
    print(f"  net nonlinear transfer T_E    {stats['T_E'][0]:+.3g} +- {stats['T_E'][1]:.2g}")
    print(f"  implied input dE/dt + D + F - T = {implied:.4g}, against eps = {eps:g}: "
          f"difference {implied - eps:+.3g} ({100 * (implied - eps) / eps:+.1f} %)")
    print(f"  realized input of the kicks: std about {sd_input:.2g} ({modes} forced modes, "
          f"E_f = {e_forced:.3g} in shells {kf - dk:g}..{kf + dk:g})")
    print(f"  difference / combined standard error = {(implied - eps) / se:+.2f}")


def main(argv):
    if len(argv) < 3 or argv[1] not in ("average", "plot", "check", "compare", "evolution", "budget"):
        print(__doc__)
        return 1
    if argv[1] == "evolution":
        times, png = None, None
        for i in range(3, len(argv) - 1):
            if argv[i] == "--times":
                times = [float(x) for x in argv[i + 1].split(",")]
            if argv[i] == "--png":
                png = argv[i + 1]
        if not times:
            print(__doc__)
            return 1
        cmd_evolution(argv[2], times, png)
        return 0
    if argv[1] == "compare":
        dirs, at, names, png, i = [], None, None, None, 3
        while i < len(argv):
            if argv[i] in ("--at", "--names", "--png") and i + 1 < len(argv):
                at, names, png = (float(argv[i + 1]) if argv[i] == "--at" else at,
                                  argv[i + 1].split(",") if argv[i] == "--names" else names,
                                  argv[i + 1] if argv[i] == "--png" else png)
                i += 2
            else:
                dirs.append(argv[i])
                i += 1
        if at is None or not dirs:
            print(__doc__)
            return 1
        cmd_compare(argv[2], dirs, at, names, png)
        return 0
    directory = argv[2]
    opts = {"--from": "0", "--kf": "0", "--below": None, "--above": None, "--png": None, "--title": None,
            "--compare": None, "--label": None, "--name": None}
    plateau = {"--min-pi-e": None, "--min-pi-z": None, "--spread-e": None, "--spread-z": None}
    opts.update({"--eps": None, "--dk": "1"})
    opts.update(plateau)
    conserves = "--conserves-enstrophy" in argv
    argv = [a for a in argv if a != "--conserves-enstrophy"]
    i = 3
    while i < len(argv):
        if argv[i] not in opts or i + 1 >= len(argv):
            print(__doc__)
            return 1
        opts[argv[i]] = argv[i + 1]
        i += 2
    t_from, kf = float(opts["--from"]), float(opts["--kf"])
    if argv[1] == "budget":
        if not kf or opts["--eps"] is None:
            print("budget needs --kf and --eps")
            return 1
        cmd_budget(directory, t_from, float(opts["--eps"]), kf, float(opts["--dk"]))
        return 0
    if argv[1] == "average":
        cmd_average(directory, t_from)
    elif argv[1] == "plot":
        cmd_plot(directory, t_from, kf, opts["--png"], opts["--title"], opts["--compare"], opts["--label"],
                 opts["--name"])
    else:
        if not kf:
            print("check needs --kf")
            return 1
        below = [int(x) for x in (opts["--below"] or f"1,{max(1, int(kf) // 2)}").split(",")]
        above = [int(x) for x in (opts["--above"] or f"{int(kf) + 2},{2 * int(kf)}").split(",")]
        plateau = {k: (float(opts[k]) if opts[k] is not None else None) for k in plateau}
        return 0 if cmd_check(directory, t_from, kf, below, above, conserves, plateau) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

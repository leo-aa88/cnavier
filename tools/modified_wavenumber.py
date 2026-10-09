#!/usr/bin/env python3
"""Modified wavenumbers of the solver's periodic difference schemes.

    tools/modified_wavenumber.py [--png FILE]
        print the fraction of the wavenumbers up to Nyquist (theta = k h in
        [0, pi]) that each scheme differentiates within 1 % and 10 %, and plot
        k* h against k h for the first and the second derivative (needs
        matplotlib for the plot)

A scheme acting on exp(i k x) gives i k* exp(i k x) for the first derivative
and -(k*)^2 exp(i k x) for the second; the exact values are k* = k. Explicit
centred stencils of order 2, 4 and 6 (finitediff.c), and Lele's (1992)
sixth-order tridiagonal compact schemes (FD_COMPACT6), whose symbols are
compact6_symbol() in finitediff.c.
"""
import math
import sys

SCHEMES = ("order 2", "order 4", "order 6", "compact 6")


def first(scheme, t):
    """k* h of the first derivative at theta = k h"""
    s1, s2, s3 = math.sin(t), math.sin(2 * t), math.sin(3 * t)
    if scheme == "order 2":
        return s1
    if scheme == "order 4":
        return (8 * s1 - s2) / 6
    if scheme == "order 6":
        return (45 * s1 - 9 * s2 + s3) / 30
    return (14 / 9 * s1 + 1 / 18 * s2) / (1 + 2 / 3 * math.cos(t))


def second(scheme, t):
    """(k* h)^2 of the second derivative at theta = k h"""
    c1, c2, c3 = math.cos(t), math.cos(2 * t), math.cos(3 * t)
    if scheme == "order 2":
        return 2 - 2 * c1
    if scheme == "order 4":
        return (30 - 32 * c1 + 2 * c2) / 12
    if scheme == "order 6":
        return (490 - 540 * c1 + 54 * c2 - 4 * c3) / 180
    return -(24 / 11 * (c1 - 1) + 3 / 22 * (c2 - 1)) / (1 + 4 / 11 * c1)


def resolved(f, exact, tol, n=20000):
    """Largest theta / pi up to which |f / exact - 1| <= tol"""
    for i in range(1, n + 1):
        t = math.pi * i / n
        if abs(f(t) / exact(t) - 1) > tol:
            return (i - 1) / n
    return 1.0


def main(argv):
    png = argv[argv.index("--png") + 1] if "--png" in argv else None
    print("Fraction of [0, pi] (k h up to Nyquist) differentiated within 1 % and 10 %")
    print(f"{'scheme':<11} {'1st, 1 %':>9} {'1st, 10 %':>10} {'2nd, 1 %':>9} {'2nd, 10 %':>10}  max k*h  max (k*h)^2")
    for s in SCHEMES:
        r = [resolved(lambda t: first(s, t), lambda t: t, tol) for tol in (0.01, 0.1)]
        r += [resolved(lambda t: second(s, t), lambda t: t * t, tol) for tol in (0.01, 0.1)]
        kmax = max(first(s, math.pi * i / 2000) for i in range(2001))
        print(f"{s:<11} {r[0]:9.3f} {r[1]:10.3f} {r[2]:9.3f} {r[3]:10.3f}  {kmax:7.3f}  {second(s, math.pi):11.3f}")
    if not png:
        return 0
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = [math.pi * i / 400 for i in range(401)]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(t, t, "k--", lw=0.8, label="exact")
    ax[1].plot(t, [x * x for x in t], "k--", lw=0.8, label="exact")
    for s in SCHEMES:
        ax[0].plot(t, [first(s, x) for x in t], label=s)
        ax[1].plot(t, [second(s, x) for x in t], label=s)
    ax[0].set_title("first derivative: k* h")
    ax[1].set_title("second derivative: (k* h)^2")
    for a in ax:
        a.set_xlabel("k h")
        a.set_xlim(0, math.pi)
        a.legend(fontsize=8)
    ax[0].set_ylim(0, math.pi)
    ax[1].set_ylim(0, math.pi ** 2)
    fig.tight_layout()
    fig.savefig(png, dpi=90)
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

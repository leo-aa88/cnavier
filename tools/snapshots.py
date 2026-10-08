#!/usr/bin/env python3
"""Panels of vorticity frames that cnavier wrote (needs numpy and matplotlib).

    tools/snapshots.py OUTPUT_DIR N1 [N2 ...] [--titles "a,b,..."] [--png FILE]
        plot output/vorticity-1-<N>.vtk side by side, each with its own colour
        scale, symmetric about zero; --titles labels the panels (for example
        with their times, which the VTK files do not hold)
"""
import os
import sys


def read_vtk(path):
    import numpy as np

    with open(path) as f:
        lines = f.read().split("\n")
    nx = ny = None
    for i, line in enumerate(lines):
        if line.startswith("DIMENSIONS"):
            nx, ny = int(line.split()[1]), int(line.split()[2])
        if line.startswith("LOOKUP_TABLE"):
            values = np.array(" ".join(lines[i + 1:]).split(), dtype=float)
            return values.reshape(ny, nx)
    raise SystemExit(f"{path}: no data")


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 1
    directory, frames, titles, png = argv[1], [], None, None
    i = 2
    while i < len(argv):
        if argv[i] == "--titles" and i + 1 < len(argv):
            titles = argv[i + 1].split(",")
            i += 2
        elif argv[i] == "--png" and i + 1 < len(argv):
            png = argv[i + 1]
            i += 2
        else:
            frames.append(int(argv[i]))
            i += 1
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, len(frames), figsize=(4.2 * len(frames), 4.2), squeeze=False)
    for k, n in enumerate(frames):
        w = read_vtk(os.path.join(directory, f"vorticity-1-{n}.vtk"))
        lim = abs(w).max()
        a = ax[0][k]
        a.imshow(w, origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim, extent=(0, 1, 0, 1))
        a.set_xticks([])
        a.set_yticks([])
        a.set_title(titles[k] if titles and k < len(titles) else f"frame {n}")
    fig.tight_layout()
    png = png or os.path.join(directory, "snapshots.png")
    fig.savefig(png, dpi=60)
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

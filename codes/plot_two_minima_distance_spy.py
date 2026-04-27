"""Spy plot for the distance graph in the two-minima MDGP instance."""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("codes/.matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from simple_mdgp_two_minima_section import make_atoms, make_bounds


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = make_atoms()
    edges, lower, upper = make_bounds(atoms)
    n_atoms = len(atoms)

    known = np.eye(n_atoms, dtype=bool)
    distances = np.zeros((n_atoms, n_atoms), dtype=float)
    for (i, j), lo, up in zip(edges, lower, upper, strict=True):
        known[i, j] = True
        known[j, i] = True
        distances[i, j] = distances[j, i] = 0.5 * (lo + up)

    fig, ax = plt.subplots(figsize=(6.8, 6.4))
    ax.spy(known, markersize=13, color="#111827")
    ax.set_title("Spy Plot of the 10-Atom Distance Matrix")
    ax.set_xlabel("atom index")
    ax.set_ylabel("atom index")
    ax.set_xticks(np.arange(n_atoms))
    ax.set_yticks(np.arange(n_atoms))
    ax.set_xticklabels(np.arange(1, n_atoms + 1))
    ax.set_yticklabels(np.arange(1, n_atoms + 1))
    ax.grid(True, color="#d1d5db", linewidth=0.8)

    for i, j in edges:
        if 7 in (i, j):
            ax.scatter([j, i], [i, j], marker="s", s=145, facecolors="none", edgecolors="#ef4444", linewidths=2.1)

    ax.text(
        0.5,
        -0.12,
        "Black marks: known distances. Red boxes: distances involving the varied atom 8.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
    )
    fig.tight_layout()

    fig.savefig(output_dir / "two_minima_distance_spy.png", dpi=220)
    fig.savefig(output_dir / "two_minima_distance_spy.svg")

    print(f"atoms: {n_atoms}")
    print(f"known off-diagonal distances: {len(edges)}")
    print("edges, 1-based:")
    for i, j in edges:
        print(f"  ({i + 1}, {j + 1})")
    print("Saved images/two_minima_distance_spy.png")
    print("Saved images/two_minima_distance_spy.svg")


if __name__ == "__main__":
    main()

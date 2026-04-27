"""Simple 2D section of a small MDGP instance under hyperbolic smoothing."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("codes/.matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np


def make_atoms() -> np.ndarray:
    """Return a small deterministic 10-atom configuration in R^3."""

    t = np.linspace(0.0, 2.8 * np.pi, 10)
    return np.column_stack(
        (
            1.15 * np.cos(t),
            1.15 * np.sin(t),
            0.32 * np.arange(10),
        )
    )


def make_bounds(atoms: np.ndarray, cutoff: float = 2.15, width: float = 0.10) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a sparse bounded-distance graph from a cutoff."""

    edges: list[tuple[int, int]] = []
    lower: list[float] = []
    upper: list[float] = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            distance = float(np.linalg.norm(atoms[i] - atoms[j]))
            if distance <= cutoff or j == i + 1:
                edges.append((i, j))
                lower.append(max(0.05, distance - width))
                upper.append(distance + width)
    return np.array(edges, dtype=int), np.array(lower), np.array(upper)


def phi(y: np.ndarray, lam: float, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of max(y, 0)."""

    return lam * y + np.sqrt(np.square(lam * y) + tau**2)


def theta(diff: np.ndarray, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of Euclidean distances."""

    return np.sqrt(np.sum(np.square(diff), axis=-1) + tau**2)


def mdgp_objective(
    atoms: np.ndarray,
    edges: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
) -> float:
    """Original nonsmooth bounded-distance penalty."""

    diffs = atoms[edges[:, 0]] - atoms[edges[:, 1]]
    distances = np.linalg.norm(diffs, axis=1)
    return float(np.sum(np.maximum(lower - distances, 0.0) + np.maximum(distances - upper, 0.0)))


def smoothed_mdgp_objective(
    atoms: np.ndarray,
    edges: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    lam: float,
    tau: float,
) -> float:
    """Hyperbolically smoothed bounded-distance penalty."""

    diffs = atoms[edges[:, 0]] - atoms[edges[:, 1]]
    distances = theta(diffs, tau=tau)
    return float(np.sum(phi(lower - distances, lam=lam, tau=tau) + phi(distances - upper, lam=lam, tau=tau)))


def evaluate_section(
    base_atoms: np.ndarray,
    edges: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    atom_index: int,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    tau: float | None,
    lam: float,
) -> np.ndarray:
    """Evaluate the objective varying only x,y of one atom."""

    values = np.empty_like(x_grid)
    atoms = base_atoms.copy()
    for row in range(x_grid.shape[0]):
        for col in range(x_grid.shape[1]):
            atoms[atom_index, 0] = x_grid[row, col]
            atoms[atom_index, 1] = y_grid[row, col]
            if tau is None:
                values[row, col] = mdgp_objective(atoms, edges, lower, upper)
            else:
                values[row, col] = smoothed_mdgp_objective(atoms, edges, lower, upper, lam=lam, tau=tau)
    return values


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = make_atoms()
    edges, lower, upper = make_bounds(atoms)

    atom_index = 7
    lam = 0.5
    true_xy = atoms[atom_index, :2].copy()
    span = 2.25
    axis = np.linspace(-span, span, 180)
    xx, yy = np.meshgrid(true_xy[0] + axis, true_xy[1] + axis)

    panels: list[tuple[str, float | None]] = [
        ("original", None),
        (r"$\tau=0.05$", 0.05),
        (r"$\tau=0.35$", 0.35),
        (r"$\tau=1.40$", 1.40),
    ]

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.4), constrained_layout=True)

    for ax, (title, tau) in zip(axes.ravel(), panels, strict=True):
        zz = evaluate_section(
            atoms,
            edges,
            lower,
            upper,
            atom_index=atom_index,
            x_grid=xx,
            y_grid=yy,
            tau=tau,
            lam=lam,
        )
        levels = np.linspace(float(np.min(zz)), float(np.percentile(zz, 96)), 16)
        ax.contourf(xx, yy, zz, levels=levels, cmap="viridis", alpha=0.78)
        lines = ax.contour(xx, yy, zz, levels=levels[::2], colors="black", linewidths=0.7)
        ax.clabel(lines, inline=True, fontsize=7, fmt="%.2f")
        ax.scatter([true_xy[0]], [true_xy[1]], s=55, color="crimson", zorder=5, label="true atom position")
        ax.set_title(title)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"$x_{8,1}$")
        ax.set_ylabel(r"$x_{8,2}$")
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle("2D section of a 10-atom MDGP objective", fontsize=16)
    fig.savefig(output_dir / "simple_mdgp_hyperbolic_section.png", dpi=220)
    fig.savefig(output_dir / "simple_mdgp_hyperbolic_section.svg")

    print(f"atoms: {len(atoms)}")
    print(f"edges: {len(edges)}")
    print(f"varied atom: {atom_index + 1}")
    print(f"max upper bound: {np.max(upper):.3f}")
    print("Saved images/simple_mdgp_hyperbolic_section.png")
    print("Saved images/simple_mdgp_hyperbolic_section.svg")


if __name__ == "__main__":
    main()

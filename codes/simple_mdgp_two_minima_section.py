"""Small MDGP section with two local minima that merge under smoothing."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("codes/.matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


def make_atoms() -> np.ndarray:
    """Return a deterministic 10-atom configuration in R^3."""

    atoms = np.array(
        [
            [-1.00, 0.00, 0.00],
            [1.00, 0.00, 0.00],
            [-1.35, -1.15, 0.25],
            [1.35, -1.15, 0.25],
            [-1.35, 1.15, -0.25],
            [1.35, 1.15, -0.25],
            [0.00, -1.75, 0.45],
            [0.00, 0.72, 0.00],
            [0.00, 1.75, -0.45],
            [0.00, 0.00, 1.10],
        ],
        dtype=float,
    )
    return atoms


def make_bounds(atoms: np.ndarray, width: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build bounds that make atom 8 ambiguous in the plotted section."""

    varied = 7
    distance_edges = [(0, varied), (1, varied)]
    fixed_edges = [(0, 2), (2, 6), (6, 3), (3, 1), (0, 4), (4, 8), (8, 5), (5, 1), (0, 9), (1, 9)]
    edges = np.array(distance_edges + fixed_edges, dtype=int)

    lower = []
    upper = []
    for i, j in edges:
        distance = float(np.linalg.norm(atoms[i] - atoms[j]))
        lower.append(max(0.05, distance - width))
        upper.append(distance + width)
    return edges, np.array(lower), np.array(upper)


def phi(y: np.ndarray, lam: float, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of max(y, 0)."""

    return lam * y + np.sqrt(np.square(lam * y) + tau**2)


def theta(diff: np.ndarray, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of Euclidean distances."""

    return np.sqrt(np.sum(np.square(diff), axis=-1) + tau**2)


def original_objective(
    atoms: np.ndarray,
    edges: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
) -> float:
    diffs = atoms[edges[:, 0]] - atoms[edges[:, 1]]
    distances = np.linalg.norm(diffs, axis=1)
    return float(np.sum(np.maximum(lower - distances, 0.0) + np.maximum(distances - upper, 0.0)))


def smoothed_objective(
    atoms: np.ndarray,
    edges: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    lam: float,
    tau: float,
) -> float:
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
    values = np.empty_like(x_grid)
    atoms = base_atoms.copy()
    for row in range(x_grid.shape[0]):
        for col in range(x_grid.shape[1]):
            atoms[atom_index, 0] = x_grid[row, col]
            atoms[atom_index, 1] = y_grid[row, col]
            values[row, col] = (
                original_objective(atoms, edges, lower, upper)
                if tau is None
                else smoothed_objective(atoms, edges, lower, upper, lam=lam, tau=tau)
            )
    return values


def count_low_basin_components(values: np.ndarray, fraction: float = 0.02) -> int:
    """Count connected components in a near-minimum sublevel set."""

    threshold = float(np.min(values) + fraction * (np.percentile(values, 90) - np.min(values)))
    _, count = ndimage.label(values <= threshold)
    return int(count)


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = make_atoms()
    edges, lower, upper = make_bounds(atoms)

    atom_index = 7
    lam = 0.5
    axis = np.linspace(-1.35, 1.35, 220)
    xx, yy = np.meshgrid(axis, axis)

    panels: list[tuple[str, float | None]] = [
        ("original", None),
        (r"$\tau=0.05$", 0.05),
        (r"$\tau=0.45$", 0.45),
        (r"$\tau=1.35$", 1.35),
    ]

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.6), constrained_layout=True)

    basin_counts: list[tuple[str, int]] = []
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
        basin_count = count_low_basin_components(zz)
        basin_counts.append((title, basin_count))

        levels = np.linspace(float(np.min(zz)), float(np.percentile(zz, 94)), 18)
        ax.contourf(xx, yy, zz, levels=levels, cmap="viridis", alpha=0.80)
        lines = ax.contour(xx, yy, zz, levels=levels[::3], colors="black", linewidths=0.75)
        ax.clabel(lines, inline=True, fontsize=7, fmt="%.2f")
        ax.scatter([-1.0, 1.0], [0.0, 0.0], marker="s", s=65, color="#111827", label="fixed anchors")
        ax.scatter([atoms[atom_index, 0], 0.0], [atoms[atom_index, 1], -atoms[atom_index, 1]], s=55, color="crimson", label="two feasible positions")
        ax.set_title(f"{title}; low basins = {basin_count}")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"$x_{8,1}$")
        ax.set_ylabel(r"$x_{8,2}$")
        ax.legend(loc="upper right", fontsize=7)

    fig.suptitle("Two wells in a 10-atom MDGP section, then one smoothed basin", fontsize=15)
    fig.savefig(output_dir / "simple_mdgp_two_minima_section.png", dpi=220)
    fig.savefig(output_dir / "simple_mdgp_two_minima_section.svg")

    print(f"atoms: {len(atoms)}")
    print(f"edges: {len(edges)}")
    print(f"varied atom: {atom_index + 1}")
    print(f"anchor edges for varied atom: (1, 8), (2, 8)")
    print(f"max upper bound: {np.max(upper):.3f}")
    for title, count in basin_counts:
        print(f"{title}: low basins = {count}")
    print("Saved images/simple_mdgp_two_minima_section.png")
    print("Saved images/simple_mdgp_two_minima_section.svg")


if __name__ == "__main__":
    main()

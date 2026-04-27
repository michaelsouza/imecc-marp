"""Demo for Gaussian smoothing on a one-dimensional nonconvex function."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("codes/.matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np

from gaussian_transform import gaussian_transform_1d


def rugged_function(x: np.ndarray) -> np.ndarray:
    """A smooth function with many local minima."""

    return (
        0.08 * np.square(x)
        + 0.55 * np.sin(6.0 * x)
        + 0.28 * np.sin(13.0 * x + 0.4)
        + 0.12 * np.cos(21.0 * x - 0.7)
    )


def count_local_minima(y: np.ndarray) -> int:
    """Count strict grid-local minima."""

    return int(np.sum((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:])))


def main() -> None:
    output_dir = Path("codes/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    x = np.linspace(-4.0, 4.0, 2401)
    lambdas = [0.0, 0.10, 0.20, 0.35, 0.55, 0.80]
    curves = {
        lam: gaussian_transform_1d(rugged_function, x, lam=lam, q=80)
        for lam in lambdas
    }

    print("Gaussian smoothing demo for a rugged one-dimensional function")
    print("lambda    local minima    min value")
    for lam, y in curves.items():
        print(f"{lam:5.2f} {count_local_minima(y):14d} {np.min(y):12.6f}")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(10.5, 6.0))

    colors = plt.cm.viridis(np.linspace(0.06, 0.92, len(lambdas)))
    for color, lam in zip(colors, lambdas, strict=True):
        label = rf"$\lambda={lam:.2f}$, minima={count_local_minima(curves[lam])}"
        linewidth = 2.8 if lam == 0 else 2.1
        ax.plot(x, curves[lam], color=color, linewidth=linewidth, label=label)

    ax.set_title("Gaussian transform reduces local minima")
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\langle f\rangle_\lambda(x)$")
    ax.legend(loc="upper center", ncols=2, frameon=True)
    ax.margins(x=0)
    fig.tight_layout()

    output_path = output_dir / "gaussian_transform_1d.svg"
    fig.savefig(output_path)
    print(f"\nSaved figure: {output_path}")


if __name__ == "__main__":
    main()

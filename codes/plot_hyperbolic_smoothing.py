"""Hand-drawn illustration of hyperbolic smoothing for the MDGP."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("codes/.matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np


def positive_part(y: np.ndarray) -> np.ndarray:
    """Return max(y, 0)."""

    return np.maximum(y, 0.0)


def phi(y: np.ndarray, lam: float, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of the positive-part penalty."""

    return lam * y + np.sqrt(np.square(lam * y) + tau**2)


def theta(x: np.ndarray, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of the Euclidean norm in one dimension."""

    return np.sqrt(np.square(x) + tau**2)


def configure_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#d6d0c4", linewidth=0.9, alpha=0.45)
    ax.tick_params(labelsize=11)


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.xkcd(scale=0.85, length=110, randomness=2.0)
    plt.rcParams.update(
        {
            "font.family": ["Comic Neue", "DejaVu Sans"],
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
            "svg.fonttype": "none",
        }
    )

    ink = "#111827"
    red = "#ef4444"
    teal = "#0f8f7a"
    blue = "#2563eb"
    gray = "#6b7280"
    tau_values = [1.20, 0.60, 0.25]
    tau_colors = [blue, teal, red]

    fig = plt.figure(figsize=(15.5, 8.7), facecolor="#fbf6ed")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.52], hspace=0.42, wspace=0.24)

    ax_penalty = fig.add_subplot(gs[0, 0])
    ax_norm = fig.add_subplot(gs[0, 1])
    ax_flow = fig.add_subplot(gs[1, :])

    y = np.linspace(-3.2, 3.2, 600)
    ax_penalty.plot(y, positive_part(y), color=ink, linewidth=3.0, label=r"$\max\{y,0\}$")
    for tau, color in zip(tau_values, tau_colors, strict=True):
        ax_penalty.plot(
            y,
            phi(y, lam=0.5, tau=tau),
            color=color,
            linewidth=2.2,
            label=rf"$\phi_{{1/2,{tau:g}}}(y)$",
        )
    ax_penalty.set_title("Smooth the penalty kink", fontsize=22, pad=14)
    ax_penalty.set_xlabel(r"constraint violation $y$", fontsize=15)
    ax_penalty.set_ylabel("penalty", fontsize=15)
    ax_penalty.set_xlim(-3.2, 3.2)
    ax_penalty.set_ylim(-0.15, 3.65)
    ax_penalty.legend(frameon=True, loc="upper left", fontsize=11)
    configure_axis(ax_penalty)
    ax_penalty.annotate(
        r"$\tau$ controls smoothness",
        xy=(0.0, phi(np.array([0.0]), 0.5, 1.2)[0]),
        xytext=(-1.45, 2.55),
        color=blue,
        fontsize=14,
        arrowprops={"arrowstyle": "->", "color": blue, "linewidth": 2.0},
    )

    x = np.linspace(-3.0, 3.0, 600)
    ax_norm.plot(x, np.abs(x), color=ink, linewidth=3.0, label=r"$\|x\|$")
    for tau, color in zip(tau_values, tau_colors, strict=True):
        ax_norm.plot(
            x,
            theta(x, tau),
            color=color,
            linewidth=2.2,
            label=rf"$\theta_{{{tau:g}}}(x)$",
        )
    ax_norm.set_title("Smooth the Euclidean norm", fontsize=22, pad=14)
    ax_norm.set_xlabel(r"coordinate difference $x_i-x_j$", fontsize=15)
    ax_norm.set_ylabel("distance", fontsize=15)
    ax_norm.set_xlim(-3.0, 3.0)
    ax_norm.set_ylim(-0.15, 3.55)
    ax_norm.legend(frameon=True, loc="upper center", fontsize=11)
    configure_axis(ax_norm)
    ax_norm.vlines(0, 0, 1.2, colors=blue, linestyles="dashed", linewidth=2.0)
    ax_norm.annotate(
        r"max gap = $\tau$",
        xy=(0.05, 1.2),
        xytext=(0.75, 1.95),
        color=blue,
        fontsize=14,
        arrowprops={"arrowstyle": "->", "color": blue, "linewidth": 2.0},
    )

    ax_flow.set_axis_off()
    ax_flow.set_xlim(0, 1)
    ax_flow.set_ylim(0, 1)
    xs = [0.12, 0.36, 0.61, 0.86]
    labels = [
        r"nonsmooth MDGP",
        r"large $\tau$",
        r"smaller $\tau$",
        r"$\tau \to 0$",
    ]
    subtitles = [
        r"$\max\{\cdot,0\}$ and $\|\cdot\|$",
        "easy smooth problem",
        "continuation path",
        "original objective",
    ]
    colors = [gray, blue, teal, red]
    for idx, (xpos, label, subtitle, color) in enumerate(zip(xs, labels, subtitles, colors, strict=True)):
        ax_flow.scatter([xpos], [0.56], s=860, color="#fff7ed", edgecolor=color, linewidth=3.0)
        ax_flow.text(xpos, 0.82, label, ha="center", va="center", fontsize=17, color=ink, fontweight="bold")
        ax_flow.text(xpos, 0.30, subtitle, ha="center", va="center", fontsize=13, color=color)
        ax_flow.text(xpos, 0.56, str(idx + 1), ha="center", va="center", fontsize=18, color=color, fontweight="bold")
        if idx < len(xs) - 1:
            ax_flow.annotate(
                "",
                xy=(xs[idx + 1] - 0.085, 0.56),
                xytext=(xpos + 0.085, 0.56),
                arrowprops={"arrowstyle": "->", "color": ink, "linewidth": 2.6},
            )
    ax_flow.text(
        0.5,
        0.08,
        r"Replace each pair term by "
        r"$\phi_{\lambda,\tau}(l_{ij}-\theta_\tau^{ij}(x))"
        r" + \phi_{\lambda,\tau}(\theta_\tau^{ij}(x)-u_{ij})$",
        ha="center",
        va="center",
        fontsize=15,
        color=ink,
    )

    fig.suptitle("Hyperbolic Smoothing for Distance Geometry", fontsize=30, fontweight="bold", y=0.985)
    fig.savefig(output_dir / "hyperbolic_smoothing.svg", bbox_inches="tight", pad_inches=0.18)
    fig.savefig(output_dir / "hyperbolic_smoothing.png", dpi=220, bbox_inches="tight", pad_inches=0.18)
    print(f"Saved {output_dir / 'hyperbolic_smoothing.svg'}")
    print(f"Saved {output_dir / 'hyperbolic_smoothing.png'}")


if __name__ == "__main__":
    main()

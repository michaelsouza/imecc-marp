"""One-dimensional hyperbolic smoothing demo for a piecewise-linear function."""

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
    """Hyperbolic smoothing of max(y, 0)."""

    return lam * y + np.sqrt(np.square(lam * y) + tau**2)


def piecewise_linear_function(x: np.ndarray) -> np.ndarray:
    """Piecewise-linear function with several nonsmooth kinks."""

    return (
        0.35
        + 0.18 * x
        + 0.90 * positive_part(-x - 1.65)
        + 0.65 * positive_part(x + 0.95)
        + 1.05 * positive_part(0.25 - x)
        + 0.75 * positive_part(x - 1.15)
        + 0.55 * positive_part(2.05 - x)
    )


def hyperbolic_smoothed_function(x: np.ndarray, lam: float, tau: float) -> np.ndarray:
    """Smooth the positive-part terms in ``piecewise_linear_function``."""

    return (
        0.35
        + 0.18 * x
        + 0.90 * phi(-x - 1.65, lam=lam, tau=tau)
        + 0.65 * phi(x + 0.95, lam=lam, tau=tau)
        + 1.05 * phi(0.25 - x, lam=lam, tau=tau)
        + 0.75 * phi(x - 1.15, lam=lam, tau=tau)
        + 0.55 * phi(2.05 - x, lam=lam, tau=tau)
    )


def configure_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="#d6d0c4", linewidth=0.9, alpha=0.48)
    ax.tick_params(labelsize=12)


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

    x = np.linspace(-3.0, 3.0, 1200)
    lam = 0.5
    tau_values = [0.85, 0.45, 0.20, 0.07]
    colors = ["#2563eb", "#0f8f7a", "#f59e0b", "#ef4444"]
    ink = "#111827"
    original = piecewise_linear_function(x)
    smoothed = {
        tau: hyperbolic_smoothed_function(x, lam=lam, tau=tau)
        for tau in tau_values
    }

    fig, ax = plt.subplots(figsize=(12.8, 7.2), facecolor="#fbf6ed")
    ax.set_facecolor("#fffdf8")

    ax.plot(
        x,
        original,
        color=ink,
        linewidth=3.2,
        label=r"piecewise-linear $f$",
        zorder=4,
    )
    for tau, color in zip(tau_values, colors, strict=True):
        ax.plot(
            x,
            smoothed[tau],
            color=color,
            linewidth=2.35,
            label=rf"$f_{{\lambda,\tau}}$, $\lambda=1/2$, $\tau={tau:g}$",
        )

    kink_positions = [-1.65, -0.95, 0.25, 1.15, 2.05]
    for xpos in kink_positions:
        ax.axvline(xpos, color="#9ca3af", linewidth=1.0, linestyle="--", alpha=0.55)

    ax.annotate(
        "large tau rounds the corners",
        xy=(0.25, hyperbolic_smoothed_function(np.array([0.25]), lam=lam, tau=0.85)[0]),
        xytext=(-2.55, 3.2),
        color="#2563eb",
        fontsize=15,
        arrowprops={"arrowstyle": "->", "color": "#2563eb", "linewidth": 2.2},
    )
    ax.annotate(
        r"as $\tau \to 0$, $f_{\lambda,\tau}$ follows $f$",
        xy=(1.15, hyperbolic_smoothed_function(np.array([1.15]), lam=lam, tau=0.07)[0]),
        xytext=(0.55, 1.45),
        color="#ef4444",
        fontsize=15,
        arrowprops={"arrowstyle": "->", "color": "#ef4444", "linewidth": 2.2},
    )

    ax.set_title("Hyperbolic Smoothing of a Piecewise-Linear Function", fontsize=24, pad=16)
    ax.set_xlabel(r"$x$", fontsize=17)
    ax.set_ylabel(r"$f(x)$ and $f_{\lambda,\tau}(x)$", fontsize=17)
    ax.set_xlim(-3.0, 3.0)
    all_values = np.concatenate([original, *smoothed.values()])
    ax.set_ylim(float(np.min(all_values) - 0.25), float(np.max(all_values) + 0.28))
    ax.legend(loc="upper right", frameon=True, fontsize=12)
    configure_axis(ax)

    formula = r"$f(x)=c+0.18x+\sum_r w_r[a_r x+b_r]_+$"
    smoothing_rule = (
        r"$f_{\lambda,\tau}(x)=c+0.18x+\sum_r w_r"
        r"\phi_{\lambda,\tau}(a_r x+b_r)$, "
        r"$\phi_{\lambda,\tau}(y)=\lambda y+\sqrt{\lambda^2y^2+\tau^2}$"
    )
    fig.text(0.5, 0.055, formula, ha="center", va="center", fontsize=14, color=ink)
    fig.text(0.5, 0.022, smoothing_rule, ha="center", va="center", fontsize=13, color=ink)
    fig.tight_layout(rect=(0.03, 0.10, 0.98, 0.98))

    svg_path = output_dir / "hyperbolic_smoothing_1d.svg"
    png_path = output_dir / "hyperbolic_smoothing_1d.png"
    fig.savefig(svg_path, bbox_inches="tight", pad_inches=0.18)
    fig.savefig(png_path, dpi=220, bbox_inches="tight", pad_inches=0.18)
    print(f"Saved {svg_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()

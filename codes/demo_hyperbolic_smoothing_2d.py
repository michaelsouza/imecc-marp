"""Two-dimensional contour demo for MDGP hyperbolic smoothing."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("codes/.matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np


def theta_tau(x: np.ndarray, y: np.ndarray, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of the Euclidean norm in two dimensions."""

    return np.sqrt(tau**2 + x**2 + y**2)


def phi(y: np.ndarray, lam: float, tau: float) -> np.ndarray:
    """Hyperbolic smoothing of the positive-part penalty."""

    return lam * y + np.sqrt(np.square(lam * y) + tau**2)


def mdgp_edge_term(
    x: np.ndarray,
    y: np.ndarray,
    lower: float,
    upper: float,
    lam: float,
    tau: float,
) -> np.ndarray:
    """Smoothed penalty for one bounded distance constraint."""

    r_tau = theta_tau(x, y, tau)
    return phi(lower - r_tau, lam=lam, tau=tau) + phi(r_tau - upper, lam=lam, tau=tau)


def radial_profile(
    radius: np.ndarray,
    lower: float,
    upper: float,
    lam: float,
    tau: float,
) -> np.ndarray:
    """Evaluate the same edge term as a radial function."""

    zeros = np.zeros_like(radius)
    return mdgp_edge_term(radius, zeros, lower=lower, upper=upper, lam=lam, tau=tau)


def sampled_min_second_radial_derivative(
    lower: float,
    upper: float,
    lam: float,
    tau: float,
    max_radius: float,
) -> float:
    """Return a simple sampled convexity diagnostic for the radial profile."""

    radius = np.linspace(0.0, max_radius, 2001)
    values = radial_profile(radius, lower=lower, upper=upper, lam=lam, tau=tau)
    step = radius[1] - radius[0]
    second = (values[2:] - 2.0 * values[1:-1] + values[:-2]) / step**2
    return float(np.min(second))


def configure_axis(ax: plt.Axes) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.set_xlabel(r"$x_1$", fontsize=12)
    ax.set_ylabel(r"$x_2$", fontsize=12)


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.xkcd(scale=0.78, length=95, randomness=1.7)
    plt.rcParams.update(
        {
            "font.family": ["Comic Neue", "DejaVu Sans"],
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
            "svg.fonttype": "none",
        }
    )

    lower = 0.65
    upper = 1.00
    lam = 0.5
    tau_values = [0.04, 0.25, 0.70, 1.15]
    grid = np.linspace(-1.75, 1.75, 380)
    x, y = np.meshgrid(grid, grid)

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2), facecolor="#fbf6ed")
    axes_flat = axes.ravel()

    for ax, tau in zip(axes_flat, tau_values, strict=True):
        z = mdgp_edge_term(x, y, lower=lower, upper=upper, lam=lam, tau=tau)
        min_second = sampled_min_second_radial_derivative(
            lower=lower,
            upper=upper,
            lam=lam,
            tau=tau,
            max_radius=float(np.max(grid)),
        )
        filled_levels = np.linspace(float(np.min(z)), float(np.max(z)), 18)
        line_levels = np.linspace(float(np.min(z)), float(np.max(z)), 8)[1:-1]
        ax.contourf(x, y, z, levels=filled_levels, cmap="YlGnBu_r", alpha=0.78)
        contour = ax.contour(x, y, z, levels=line_levels, colors="#111827", linewidths=1.15)
        ax.clabel(contour, inline=True, fontsize=7, fmt="%.1f")
        ax.contour(x, y, np.sqrt(x**2 + y**2), levels=[lower, upper], colors=["#ef4444", "#2563eb"], linewidths=2.2)
        ax.scatter([0.0], [0.0], s=40, color="#111827", zorder=4)
        status = "convex by corollary" if tau > upper else "not convex"
        sign = "> u" if tau > upper else "< u"
        ax.set_title(
            rf"$\tau={tau:g}$ ({sign}), {status}",
            fontsize=15,
            pad=10,
        )
        ax.text(
            0.03,
            0.04,
            rf"sampled min $g''(r)$ = {min_second:.2f}",
            transform=ax.transAxes,
            fontsize=10,
            color="#111827",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "#fffdf8", "edgecolor": "#d1d5db", "alpha": 0.88},
        )
        configure_axis(ax)

    fig.suptitle("Hyperbolic Smoothing Can Convexify an MDGP Edge Term", fontsize=24, fontweight="bold", y=0.975)
    fig.text(
        0.5,
        0.045,
        rf"One bounded distance constraint: "
        rf"$f_{{\lambda,\tau}}^{{ij}}(x)=\phi_{{\lambda,\tau}}(l-\theta_\tau(x))"
        rf"+\phi_{{\lambda,\tau}}(\theta_\tau(x)-u)$, "
        rf"$l={lower:g}$, $u={upper:g}$, $\lambda=1/2$.",
        ha="center",
        va="center",
        fontsize=13,
        color="#111827",
    )
    fig.text(
        0.5,
        0.017,
        r"Dashed intuition: small $\tau$ preserves the annular nonconvex valley; once $\tau>u$, the corollary guarantees convexity.",
        ha="center",
        va="center",
        fontsize=12,
        color="#374151",
    )
    fig.subplots_adjust(left=0.07, right=0.98, top=0.88, bottom=0.15, wspace=0.18, hspace=0.34)

    svg_path = output_dir / "hyperbolic_smoothing_2d_contours.svg"
    png_path = output_dir / "hyperbolic_smoothing_2d_contours.png"
    fig.savefig(svg_path, bbox_inches="tight", pad_inches=0.18)
    fig.savefig(png_path, dpi=220, bbox_inches="tight", pad_inches=0.18)
    print(f"Saved {svg_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()

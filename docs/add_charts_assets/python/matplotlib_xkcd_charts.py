"""Generate hand-drawn charts with Matplotlib's xkcd style.

Run from the repository root:

    python docs/add_charts_assets/python/matplotlib_xkcd_charts.py

The script writes SVG files into:

    docs/add_charts_assets/images/
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT_DIR = Path(__file__).resolve().parents[1] / "images"


def save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / name, bbox_inches="tight", format="svg")
    plt.close(fig)


def bar_chart() -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(
        ["BP", "SBP", "DGPw"],
        [30, 18, 22],
        color=["#ffb4a2", "#99f6e4", "#bfdbfe"],
        edgecolor="#111827",
        linewidth=2,
        hatch="////",
    )
    ax.set_title("Runtime Comparison")
    ax.set_ylabel("time (s)")
    save(fig, "bar_chart.svg")


def line_chart() -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(1, 8)
    ax.plot(x, [42, 31, 25, 19, 15, 12, 10], marker="o", linewidth=2.5, label="BP")
    ax.plot(x, [42, 28, 20, 14, 10, 8, 7], marker="s", linewidth=2.5, label="SBP")
    ax.set_title("Pruning Over Iterations")
    ax.set_xlabel("iteration")
    ax.set_ylabel("active branches")
    ax.legend(frameon=False)
    save(fig, "line_chart.svg")


def scatter_plot() -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.array([8, 12, 14, 17, 18, 22, 24, 29, 32])
    y = np.array([1.5, 2.4, 2.1, 3.2, 3.7, 4.1, 4.8, 6.2, 6.8])
    ax.scatter(x, y, s=110, color="#fb7185", edgecolor="#111827", linewidth=1.5)
    ax.set_title("Instance Size vs Runtime")
    ax.set_xlabel("number of vertices")
    ax.set_ylabel("runtime (s)")
    save(fig, "scatter_plot.svg")


def pie_chart() -> None:
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pie(
        [48, 32, 20],
        labels=["pruned", "feasible", "symmetric"],
        colors=["#99f6e4", "#bfdbfe", "#fecdd3"],
        wedgeprops={"edgecolor": "#111827", "linewidth": 2},
        autopct="%1.0f%%",
    )
    ax.set_title("Branch Outcomes")
    save(fig, "pie_chart.svg")


def boxplot() -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    data = [
        [25, 29, 31, 34, 38, 41],
        [14, 16, 17, 19, 22, 23],
        [18, 20, 21, 24, 27, 29],
    ]
    box = ax.boxplot(data, patch_artist=True, labels=["BP", "SBP", "DGPw"])
    for patch, color in zip(box["boxes"], ["#ffb4a2", "#99f6e4", "#bfdbfe"]):
        patch.set(facecolor=color, edgecolor="#111827", linewidth=2)
    for item in box["whiskers"] + box["caps"] + box["medians"]:
        item.set(color="#111827", linewidth=2)
    ax.set_title("Runtime Distribution")
    ax.set_ylabel("time (s)")
    save(fig, "boxplot.svg")


def main() -> None:
    with plt.xkcd(scale=1.0, length=90, randomness=2):
        plt.rcParams["figure.facecolor"] = "#f8faf7"
        plt.rcParams["axes.facecolor"] = "#f8faf7"
        bar_chart()
        line_chart()
        scatter_plot()
        pie_chart()
        boxplot()


if __name__ == "__main__":
    main()


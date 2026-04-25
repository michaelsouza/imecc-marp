"""Generate a hand-drawn DGP graph with NetworkX and Matplotlib.

Run from the repository root:

    python docs/add_images_assets/python/dgp_graph_networkx.py

The script writes:

    docs/add_images_assets/images/dgp_graph_networkx.svg
"""

from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx


OUT = Path(__file__).resolve().parents[1] / "images" / "dgp_graph_networkx.svg"


def main() -> None:
    graph = nx.Graph()
    graph.add_edges_from(
        [
            (1, 2),
            (2, 3),
            (1, 3),
            (1, 4),
            (2, 4),
            (3, 4),
            (2, 5),
            (3, 5),
            (4, 5),
        ]
    )
    pos = {
        1: (0.0, 0.0),
        2: (1.6, 0.0),
        3: (0.8, 1.25),
        4: (2.35, 1.05),
        5: (3.3, 0.15),
    }

    with plt.xkcd(scale=1.0, length=90, randomness=2):
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        fig.patch.set_facecolor("#f8faf7")
        ax.set_facecolor("#f8faf7")

        nx.draw_networkx_edges(graph, pos, ax=ax, width=2.5, edge_color="#111827")
        nx.draw_networkx_nodes(
            graph,
            pos,
            ax=ax,
            node_size=900,
            node_color=["#99f6e4", "#99f6e4", "#99f6e4", "#fecdd3", "#fecdd3"],
            edgecolors="#111827",
            linewidths=2.3,
        )
        nx.draw_networkx_labels(graph, pos, ax=ax, font_size=15, font_weight="bold")

        ax.text(0.0, 1.75, "DGP graph", fontsize=22, weight="bold")
        ax.text(2.15, 1.45, "new vertices", color="#be123c", fontsize=13)
        ax.text(0.05, -0.45, "known clique", color="#0f766e", fontsize=13)
        ax.set_axis_off()

        OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUT, bbox_inches="tight", format="svg")


if __name__ == "__main__":
    main()


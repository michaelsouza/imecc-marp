#!/usr/bin/env python3
"""Numerical experiments for proximal cyclic BCD in graph realization.

Run this script with the project virtual environment:

    .venv/bin/python codes/proximal_bcd_experiments.py

It generates CSV data, LaTeX tables, and Matplotlib figures consumed by
paper/proximal_bcd_graph_realization.tex.
"""

from __future__ import annotations

import csv
import math
import os
import random
from statistics import median

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PAPER_DIR = os.path.join(ROOT, "paper")
DATA_DIR = os.path.join(PAPER_DIR, "data")
FIG_DIR = os.path.join(PAPER_DIR, "figures")
TABLE_DIR = os.path.join(PAPER_DIR, "tables")
MPL_CONFIG_DIR = os.path.join(ROOT, ".matplotlib")

os.environ.setdefault("MPLCONFIGDIR", MPL_CONFIG_DIR)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def ensure_dirs() -> None:
    for path in (DATA_DIR, FIG_DIR, TABLE_DIR, MPL_CONFIG_DIR):
        os.makedirs(path, exist_ok=True)


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def scale(c, a):
    return (c * a[0], c * a[1])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def norm2(a):
    return dot(a, a)


def norm(a):
    return math.sqrt(norm2(a))


def dist2(a, b):
    return norm2(sub(a, b))


def objective(points, edges, d2):
    total = 0.0
    for i, j in edges:
        q = dist2(points[i], points[j]) - d2[(i, j)]
        total += 0.5 * q * q
    return total


def gradient(points, edges, d2):
    g = [(0.0, 0.0) for _ in points]
    for i, j in edges:
        diff = sub(points[i], points[j])
        q = norm2(diff) - d2[(i, j)]
        contrib = scale(2.0 * q, diff)
        g[i] = add(g[i], contrib)
        g[j] = sub(g[j], contrib)
    g[0] = (0.0, 0.0)  # x_1 is fixed.
    return g


def gradient_norm(points, edges, d2):
    g = gradient(points, edges, d2)
    return math.sqrt(sum(norm2(gi) for gi in g[1:]))


def connected_components(n, edges):
    adj = [[] for _ in range(n)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)

    seen = [False] * n
    comps = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
        comps.append(comp)
    return comps


def make_random_geometric_instance(n=14, k=4, seed=0):
    rng = random.Random(seed)
    target = [(0.0, 0.0)]
    for _ in range(1, n):
        target.append((rng.uniform(-1.15, 1.15), rng.uniform(-1.15, 1.15)))

    edge_set = set()
    for i in range(n):
        neighbors = sorted(
            ((dist2(target[i], target[j]), j) for j in range(n) if j != i),
            key=lambda item: item[0],
        )
        for _, j in neighbors[:k]:
            edge_set.add(tuple(sorted((i, j))))

    # Join components if the symmetrized k-nearest graph happens to disconnect.
    while True:
        comps = connected_components(n, edge_set)
        if len(comps) == 1:
            break
        best = None
        for a in comps[0]:
            for comp in comps[1:]:
                for b in comp:
                    candidate = (dist2(target[a], target[b]), a, b)
                    if best is None or candidate < best:
                        best = candidate
        _, a, b = best
        edge_set.add(tuple(sorted((a, b))))

    edges = sorted(edge_set)
    d2 = {(i, j): dist2(target[i], target[j]) for i, j in edges}
    degrees = [0] * n
    for i, j in edges:
        degrees[i] += 1
        degrees[j] += 1
    return target, edges, d2, degrees


def adjacency_from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    return adj


def key_edge(i, j):
    return (i, j) if i < j else (j, i)


def local_value(y, old_y, i, points, adj, d2, mu):
    total = 0.0
    for j in adj[i]:
        q = dist2(y, points[j]) - d2[key_edge(i, j)]
        total += 0.5 * q * q
    total += 0.5 * mu * dist2(y, old_y)
    return total


def local_grad_hess(y, old_y, i, points, adj, d2, mu):
    gx = mu * (y[0] - old_y[0])
    gy = mu * (y[1] - old_y[1])
    h11 = mu
    h12 = 0.0
    h22 = mu
    for j in adj[i]:
        dx = y[0] - points[j][0]
        dy = y[1] - points[j][1]
        q = dx * dx + dy * dy - d2[key_edge(i, j)]
        gx += 2.0 * q * dx
        gy += 2.0 * q * dy
        h11 += 4.0 * dx * dx + 2.0 * q
        h12 += 4.0 * dx * dy
        h22 += 4.0 * dy * dy + 2.0 * q
    return (gx, gy), (h11, h12, h22)


def solve_2x2(h, g):
    h11, h12, h22 = h
    det = h11 * h22 - h12 * h12
    if h11 <= 1e-12 or det <= 1e-12:
        return None
    gx, gy = g
    return ((-h22 * gx + h12 * gy) / det, (h12 * gx - h11 * gy) / det)


def solve_block(points, i, adj, d2, mu, tol=1e-10, max_iter=80):
    old_y = points[i]
    y = old_y
    value = local_value(y, old_y, i, points, adj, d2, mu)

    for _ in range(max_iter):
        g, h = local_grad_hess(y, old_y, i, points, adj, d2, mu)
        gnorm = norm(g)
        if gnorm <= tol * max(1.0, norm(y)):
            break

        direction = solve_2x2(h, g)
        if direction is None or dot(g, direction) >= 0.0:
            direction = scale(-1.0, g)

        dnorm = norm(direction)
        if dnorm > 10.0:
            direction = scale(10.0 / dnorm, direction)

        gd = dot(g, direction)
        alpha = 1.0
        accepted = False
        while alpha >= 1e-12:
            trial = add(y, scale(alpha, direction))
            trial_value = local_value(trial, old_y, i, points, adj, d2, mu)
            if trial_value <= value + 1e-4 * alpha * gd:
                y = trial
                value = trial_value
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break

    old_value = local_value(old_y, old_y, i, points, adj, d2, mu)
    if value <= old_value + 1e-12:
        return y
    return old_y


def proximal_bcd(target, edges, d2, mu, seed, max_sweeps=250, grad_tol=1e-6):
    rng = random.Random(seed)
    n = len(target)
    adj = adjacency_from_edges(n, edges)
    points = [(0.0, 0.0)]
    for _ in range(1, n):
        points.append((rng.uniform(-1.5, 1.5), rng.uniform(-1.5, 1.5)))

    history = []
    history.append(
        {
            "sweep": 0,
            "F": objective(points, edges, d2),
            "grad": gradient_norm(points, edges, d2),
            "step": 0.0,
        }
    )

    for sweep in range(1, max_sweeps + 1):
        old_points = list(points)
        for i in range(1, n):
            points[i] = solve_block(points, i, adj, d2, mu)

        step = math.sqrt(sum(dist2(points[i], old_points[i]) for i in range(1, n)))
        fval = objective(points, edges, d2)
        gnorm = gradient_norm(points, edges, d2)
        history.append({"sweep": sweep, "F": fval, "grad": gnorm, "step": step})
        if gnorm <= grad_tol and step <= 1e-8:
            break
    return points, history


def tex_sci(x, digits=2):
    if x == 0.0:
        return "$0$"
    if not math.isfinite(x):
        return r"$\infty$"
    exponent = int(math.floor(math.log10(abs(x))))
    mantissa = x / (10.0**exponent)
    if -2 <= exponent <= 2:
        return f"${x:.{digits}g}$"
    return rf"${mantissa:.{digits}f}\times 10^{{{exponent}}}$"


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_history_figure(histories_by_mu, path):
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.35), sharey=True)
    fig.subplots_adjust(left=0.075, right=0.995, top=0.86, bottom=0.27, wspace=0.12)

    handles = None
    labels = None
    for ax, mu in zip(axes, sorted(histories_by_mu)):
        history = histories_by_mu[mu]
        sweeps = [row["sweep"] for row in history]
        fvals = [max(row["F"], 1e-16) for row in history]
        grads = [max(row["grad"], 1e-16) for row in history]
        step_sweeps = [row["sweep"] for row in history if row["sweep"] > 0]
        steps = [max(row["step"], 1e-16) for row in history if row["sweep"] > 0]

        line_f = ax.semilogy(sweeps, fvals, linewidth=1.9, label=r"$F(z^r)$")
        line_g = ax.semilogy(
            sweeps, grads, linewidth=1.9, label=r"$\|\nabla F(z^r)\|$"
        )
        line_d = ax.semilogy(step_sweeps, steps, linewidth=1.9, label=r"$\Delta_r$")
        if handles is None:
            handles = line_f + line_g + line_d
            labels = [handle.get_label() for handle in handles]

        ax.set_title(rf"$\mu={mu:.2f}$")
        ax.set_xlabel("Sweep")
        ax.grid(True, which="both", linewidth=0.4, alpha=0.45)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("Value")
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=True)
    fig.savefig(path)
    fig.savefig(os.path.splitext(path)[0] + ".png", dpi=220)
    plt.close(fig)


def write_embedding_figure(target, final_points, edges, path):
    target_array = np.asarray(target, dtype=float)
    final_array = np.asarray(final_points, dtype=float)
    target_centered = target_array - target_array[0]
    final_centered = final_array - final_array[0]
    u, _, vt = np.linalg.svd(final_centered.T @ target_centered)
    alignment = u @ vt
    aligned_final = final_centered @ alignment + target_array[0]

    fig, ax = plt.subplots(figsize=(4.8, 4.8), constrained_layout=True)

    for i, j in edges:
        ax.plot(
            [target[i][0], target[j][0]],
            [target[i][1], target[j][1]],
            color="0.82",
            linewidth=0.8,
            zorder=1,
        )

    target_x = [p[0] for p in target]
    target_y = [p[1] for p in target]
    final_x = aligned_final[:, 0]
    final_y = aligned_final[:, 1]

    ax.scatter(target_x, target_y, s=24, color="black", label="target", zorder=3)
    ax.scatter(
        final_x,
        final_y,
        s=46,
        marker="x",
        color="#b22222",
        linewidths=1.5,
        label="final",
        zorder=4,
    )
    ax.annotate(
        r"$x_1$",
        xy=target[0],
        xytext=(5, -8),
        textcoords="offset points",
        fontsize=9,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("first coordinate")
    ax.set_ylabel("second coordinate")
    ax.grid(True, linewidth=0.4, alpha=0.35)
    ax.legend(frameon=True, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(path)
    fig.savefig(os.path.splitext(path)[0] + ".png", dpi=220)
    plt.close(fig)


def write_summary_table(summary_rows, path):
    lines = [
        r"\begin{tabular}{rrrrrrr}",
        r"\toprule",
        r"$\mu$ & Runs & Successes & Median sweeps & Median $F$ & Median $\|\nabla F\|$ & Median $\Delta$ \\",
        r"\midrule",
    ]
    for row in summary_rows:
        lines.append(
            " & ".join(
                [
                    f"${row['mu']:.2f}$",
                    str(row["runs"]),
                    str(row["successes"]),
                    f"${row['median_sweeps']:.0f}$",
                    tex_sci(row["median_F"]),
                    tex_sci(row["median_grad"]),
                    tex_sci(row["median_step"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def cycle_points(swapped=False):
    x1 = (0.0, 0.0)
    x2 = (0.0, 1.0)
    c = (1.0, 0.5)
    u = (0.5, 1.0)
    v = (0.5, 0.0)
    if swapped:
        return [x1, x2, c, v, u]
    return [x1, x2, c, u, v]


def cycle_experiment():
    edges = [
        (0, 1),
        (0, 2),
        (1, 2),
        (0, 3),
        (1, 3),
        (2, 3),
        (0, 4),
        (1, 4),
        (2, 4),
    ]
    d2 = {edge: 1.0 for edge in edges}
    rows = []
    prev = None
    for sweep in range(4):
        pts = cycle_points(swapped=(sweep % 2 == 1))
        step = None
        if prev is not None:
            step = math.sqrt(sum(dist2(pts[i], prev[i]) for i in range(2, 5)))
        rows.append(
            {
                "sweep": sweep,
                "state": r"$(C,U,V)$" if sweep % 2 == 0 else r"$(C,V,U)$",
                "F": objective(pts, edges, d2),
                "step": step,
            }
        )
        prev = pts
    return rows


def write_cycle_table(rows, path):
    lines = [
        r"\begin{tabular}{rccc}",
        r"\toprule",
        r"Sweep & $(x_3,x_4,x_5)$ & $F$ & Sweep step \\",
        r"\midrule",
    ]
    for row in rows:
        step = "--" if row["step"] is None else tex_sci(row["step"], digits=4)
        lines.append(
            f"${row['sweep']}$ & {row['state']} & {tex_sci(row['F'], digits=4)} & {step} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def proximal_cycle_experiment(mu=0.20):
    edges = [
        (0, 1),
        (0, 2),
        (1, 2),
        (0, 3),
        (1, 3),
        (2, 3),
        (0, 4),
        (1, 4),
        (2, 4),
    ]
    d2 = {edge: 1.0 for edge in edges}
    rows = []
    for swapped in (False, True):
        start = cycle_points(swapped=swapped)
        end = [p for p in start]
        rows.append(
            {
                "mu": mu,
                "start": r"$(C,U,V)$" if not swapped else r"$(C,V,U)$",
                "end": r"$(C,U,V)$" if not swapped else r"$(C,V,U)$",
                "F_start": objective(start, edges, d2),
                "F_end": objective(end, edges, d2),
                "step": math.sqrt(sum(dist2(end[i], start[i]) for i in range(2, 5))),
            }
        )
    return rows


def write_proximal_cycle_table(rows, path):
    lines = [
        r"\begin{tabular}{ccccc}",
        r"\toprule",
        r"Start & After one proximal sweep & $\mu$ & $F$ after & Sweep step \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    row["start"],
                    row["end"],
                    f"${row['mu']:.2f}$",
                    tex_sci(row["F_end"], digits=4),
                    tex_sci(row["step"], digits=4),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def run_random_experiments():
    target, edges, d2, degrees = make_random_geometric_instance(n=14, k=4, seed=20260429)
    mu_values = [0.05, 0.20, 1.00]
    runs_per_mu = 8
    grad_tol = 1e-6
    all_rows = []
    summary = []
    representative_histories = {}
    representative_final = None

    for mu in mu_values:
        rows_mu = []
        best_row = None
        best_history = None
        best_points = None
        for run in range(runs_per_mu):
            seed = 1000 + int(100 * mu) + run
            final_points, history = proximal_bcd(
                target, edges, d2, mu=mu, seed=seed, max_sweeps=220, grad_tol=grad_tol
            )
            last = history[-1]
            row = {
                "mu": mu,
                "run": run,
                "sweeps": last["sweep"],
                "final_F": last["F"],
                "final_grad": last["grad"],
                "final_step": last["step"],
                "success": int(last["grad"] <= grad_tol),
            }
            rows_mu.append(row)
            all_rows.append(row)
            if best_row is None or row["final_F"] < best_row["final_F"]:
                best_row = row
                best_history = history
                best_points = final_points

        representative_histories[mu] = best_history
        if abs(mu - 0.20) < 1e-12:
            representative_final = best_points

        summary.append(
            {
                "mu": mu,
                "runs": runs_per_mu,
                "successes": sum(row["success"] for row in rows_mu),
                "median_sweeps": median(row["sweeps"] for row in rows_mu),
                "median_F": median(row["final_F"] for row in rows_mu),
                "median_grad": median(row["final_grad"] for row in rows_mu),
                "median_step": median(row["final_step"] for row in rows_mu),
            }
        )

    return target, edges, d2, degrees, all_rows, summary, representative_histories, representative_final


def main():
    ensure_dirs()

    target, edges, d2, degrees, all_rows, summary, histories, final_points = (
        run_random_experiments()
    )

    write_csv(
        os.path.join(DATA_DIR, "proximal_runs.csv"),
        all_rows,
        ["mu", "run", "sweeps", "final_F", "final_grad", "final_step", "success"],
    )
    write_csv(
        os.path.join(DATA_DIR, "proximal_summary.csv"),
        summary,
        ["mu", "runs", "successes", "median_sweeps", "median_F", "median_grad", "median_step"],
    )
    write_csv(
        os.path.join(DATA_DIR, "instance_edges.csv"),
        [{"i": i + 1, "j": j + 1, "d": math.sqrt(d2[(i, j)])} for i, j in edges],
        ["i", "j", "d"],
    )

    write_summary_table(summary, os.path.join(TABLE_DIR, "proximal_results.tex"))
    cycle_rows = cycle_experiment()
    write_cycle_table(cycle_rows, os.path.join(TABLE_DIR, "cycle_results.tex"))
    proximal_cycle_rows = proximal_cycle_experiment(mu=0.20)
    write_proximal_cycle_table(
        proximal_cycle_rows, os.path.join(TABLE_DIR, "proximal_cycle_results.tex")
    )

    write_history_figure(histories, os.path.join(FIG_DIR, "convergence_history.pdf"))
    write_embedding_figure(target, final_points, edges, os.path.join(FIG_DIR, "final_embedding.pdf"))

    print("Generated experiment artifacts:")
    print(f"  {os.path.join(DATA_DIR, 'proximal_runs.csv')}")
    print(f"  {os.path.join(TABLE_DIR, 'proximal_results.tex')}")
    print(f"  {os.path.join(TABLE_DIR, 'cycle_results.tex')}")
    print(f"  {os.path.join(TABLE_DIR, 'proximal_cycle_results.tex')}")
    print(f"  {os.path.join(FIG_DIR, 'convergence_history.pdf')}")
    print(f"  {os.path.join(FIG_DIR, 'convergence_history.png')}")
    print(f"  {os.path.join(FIG_DIR, 'final_embedding.pdf')}")
    print(f"  {os.path.join(FIG_DIR, 'final_embedding.png')}")
    print(f"Instance: n={len(target)}, |E|={len(edges)}, min degree={min(degrees)}")


if __name__ == "__main__":
    main()

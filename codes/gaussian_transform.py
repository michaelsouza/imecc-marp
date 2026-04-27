"""Gaussian smoothing utilities.

The routines in this module implement the Gaussian transform used by
More and Wu's continuation approach, plus a one-dimensional helper that is
useful for experiments and figures.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.polynomial.hermite_e import hermegauss


def gaussian_transform_1d(
    f: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray | float,
    lam: float,
    q: int = 40,
) -> np.ndarray | float:
    """Approximate the one-dimensional Gaussian transform of ``f``.

    This computes

        <f>_lambda(x) = E[f(x + lambda Z / sqrt(2))],

    where Z is a standard normal random variable. This matches the
    normalization in More and Wu's Gaussian transform

        1 / (sqrt(pi) lambda) int f(y) exp(-(y-x)^2 / lambda^2) dy.

    Parameters
    ----------
    f:
        Vectorized scalar function.
    x:
        Evaluation points.
    lam:
        Smoothing parameter. ``lam = 0`` returns the original function.
    q:
        Number of Gauss-Hermite quadrature nodes.
    """

    values = np.asarray(x, dtype=float)
    scalar_input = values.ndim == 0

    if lam < 0:
        raise ValueError("lam must be nonnegative")
    if q < 1:
        raise ValueError("q must be positive")
    if lam == 0:
        result = f(values)
        return float(result) if scalar_input else result

    nodes, weights = hermegauss(q)
    normalized_weights = weights / np.sqrt(2.0 * np.pi)
    samples = values[..., np.newaxis] + lam * nodes / np.sqrt(2.0)
    result = np.sum(normalized_weights * f(samples), axis=-1)

    return float(result) if scalar_input else result


def bounded_distance_penalty(
    r: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
) -> np.ndarray:
    """Evaluate the bounded DGP radial penalty h_ij(r)."""

    r2 = np.square(r)
    lower2 = np.square(lower)
    upper2 = np.square(upper)
    lower_violation = np.minimum((r2 - lower2) / lower2, 0.0)
    upper_violation = np.maximum((r2 - upper2) / upper2, 0.0)
    return np.square(lower_violation) + np.square(upper_violation)


def gaussian_transform_dgp(
    x: np.ndarray,
    edges: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    lam: float,
    q: int = 10,
    eps: float = 1e-12,
) -> float:
    """Approximate the Gauss-Hermite transform of the bounded DGP objective.

    Parameters
    ----------
    x:
        Coordinates with shape ``(m, dim)``.
    edges:
        Integer array with shape ``(s, 2)``. Each row stores one constrained
        pair ``(i, j)``.
    lower, upper:
        Lower and upper distance bounds for each edge.
    lam:
        Smoothing parameter.
    q:
        Number of Gauss-Hermite quadrature nodes.
    eps:
        Small value used to avoid division by zero.
    """

    if lam < 0:
        raise ValueError("lam must be nonnegative")
    if q < 1:
        raise ValueError("q must be positive")

    x = np.asarray(x, dtype=float)
    edges = np.asarray(edges, dtype=int)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)

    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError("edges must have shape (s, 2)")
    if lower.shape != (edges.shape[0],) or upper.shape != (edges.shape[0],):
        raise ValueError("lower and upper must have shape (s,)")

    diffs = x[edges[:, 0]] - x[edges[:, 1]]
    r = np.linalg.norm(diffs, axis=1)

    if lam == 0:
        return float(np.sum(bounded_distance_penalty(r, lower, upper)))

    nodes, weights = hermegauss(q)
    normalized_weights = weights / np.sqrt(2.0 * np.pi)

    samples = r[:, np.newaxis] + lam * nodes[np.newaxis, :]
    h_values = bounded_distance_penalty(samples, lower[:, np.newaxis], upper[:, np.newaxis])
    edge_values = np.sum(
        normalized_weights[np.newaxis, :] * samples * h_values,
        axis=1,
    )
    return float(np.sum(edge_values / np.maximum(r, eps)))

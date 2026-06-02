#!/usr/bin/env python3
"""Verify a locally attracting Powell-type DGP cycle candidate.

This is a scalar DGP instance with five free vertices and, for each free
vertex, 30 coincident anchors at 0.  The free-free graph is K5.  The two
states

    A = ( 1, -1,  1, -1,  1),      B = -A

form a legal two-cycle for exact cyclic coordinate minimization if the ties
at the limit are resolved along the displayed branch.  Off the tie surface,
the branch has a certified open cone of perturbations whose minimizers are
unique and whose two-sweep linearization is contracting.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

import numpy as np


N = 5
ANCHORS_PER_VERTEX = 30
A = tuple(Fraction(v) for v in (1, -1, 1, -1, 1))
B = tuple(-v for v in A)
PAIRS = list(combinations(range(N), 2))


def dec(value: str) -> Fraction:
    return Fraction(value)


# Six independent free-free squared distances.  The remaining four are
# defined by the linear sign-balance equations that force the tied
# subproblems to have derivative proportional to t(t^2 - 1).
BETA: dict[tuple[int, int], Fraction] = {
    (0, 1): dec("11.667256816711"),
    (0, 2): dec("3.348281865181"),
    (0, 3): dec("0.026359083326"),
    (1, 2): dec("0.286099935630"),
    (1, 3): dec("2.086006466135"),
    (2, 3): dec("7.761224706209"),
}
BETA[(0, 4)] = BETA[(0, 1)] - BETA[(0, 2)] + BETA[(0, 3)]
BETA[(1, 4)] = BETA[(0, 1)] - BETA[(1, 2)] + BETA[(1, 3)]
BETA[(2, 4)] = BETA[(0, 2)] - BETA[(1, 2)] + BETA[(2, 3)]
BETA[(3, 4)] = BETA[(0, 3)] - BETA[(1, 3)] + BETA[(2, 3)]


def edge(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def beta(i: int, j: int) -> Fraction:
    return BETA[edge(i, j)]


SIGMA = tuple(
    (Fraction(ANCHORS_PER_VERTEX + 16) - sum(beta(i, j) for j in range(N) if j != i))
    / ANCHORS_PER_VERTEX
    for i in range(N)
)


def partial_state(start: tuple[Fraction, ...], target: tuple[Fraction, ...], block: int) -> tuple[Fraction, ...]:
    current = list(start)
    for i in range(block):
        current[i] = target[i]
    return tuple(current)


def half_derivative_coefficients(block: int, context: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    """Coefficients of half the scalar subproblem derivative.

    The polynomial is returned as coefficients of t^3, t^2, t, 1.
    """
    coeff = [Fraction(ANCHORS_PER_VERTEX), Fraction(0), -ANCHORS_PER_VERTEX * SIGMA[block], Fraction(0)]
    for j in range(N):
        if j == block:
            continue
        c = context[j]
        d2 = beta(block, j)
        coeff[0] += 1
        coeff[1] += -3 * c
        coeff[2] += 3 * c**2 - d2
        coeff[3] += d2 * c - c**3
    return tuple(coeff)


def objective(x: tuple[Fraction, ...]) -> Fraction:
    total = Fraction(15) * sum((x[i] ** 2 - SIGMA[i]) ** 2 for i in range(N))
    for i, j in PAIRS:
        total += Fraction(1, 2) * ((x[i] - x[j]) ** 2 - beta(i, j)) ** 2
    return total


def gradient(x: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    grad = [60 * x[i] * (x[i] ** 2 - SIGMA[i]) for i in range(N)]
    for i, j in PAIRS:
        diff = x[i] - x[j]
        residual = diff**2 - beta(i, j)
        contribution = 2 * residual * diff
        grad[i] += contribution
        grad[j] -= contribution
    return tuple(grad)


def as_float_matrix(values: list[list[Fraction]] | tuple[tuple[Fraction, ...], ...]) -> np.ndarray:
    return np.array([[float(v) for v in row] for row in values], dtype=float)


def branch_jacobian_and_preferences(start: tuple[Fraction, ...], target: tuple[Fraction, ...]) -> tuple[np.ndarray, np.ndarray]:
    current = np.array([float(v) for v in start], dtype=float)
    target_f = np.array([float(v) for v in target], dtype=float)
    beta_f = as_float_matrix([[Fraction(0) if i == j else beta(i, j) for j in range(N)] for i in range(N)])
    state_jac = np.eye(N)
    rows: list[np.ndarray] = []
    gt = 2.0 * (ANCHORS_PER_VERTEX + N - 1)

    for i in range(N):
        intended = target_f[i]
        other = -intended
        row = np.zeros(N)
        for j in range(N):
            if j == i:
                continue
            c = current[j]
            r_intended = (intended - c) ** 2 - beta_f[i, j]
            r_other = (other - c) ** 2 - beta_f[i, j]
            d_delta_dc = -2 * r_intended * (intended - c) + 2 * r_other * (other - c)
            row += d_delta_dc * state_jac[j, :]
        rows.append(row)

        new_row = np.zeros(N)
        for j in range(N):
            if j == i:
                continue
            c = current[j]
            dt_dc = (3 * (intended - c) ** 2 - beta_f[i, j]) / gt
            new_row += dt_dc * state_jac[j, :]
        state_jac[i, :] = new_row
        current[i] = intended

    return state_jac, np.vstack(rows)


def solve_cone_margin(k_matrix: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Maximize s subject to K v >= s and |v_i| <= 1 by vertex enumeration."""
    constraints: list[tuple[str, int, np.ndarray, float]] = []
    for r in range(k_matrix.shape[0]):
        row = np.zeros(N + 1)
        row[:N] = k_matrix[r]
        row[N] = -1
        constraints.append(("ineq", r, row, 0.0))
    for i in range(N):
        row = np.zeros(N + 1)
        row[i] = 1
        constraints.append(("upper", i, row, 1.0))
        constraints.append(("lower", i, row, -1.0))

    best_s = -float("inf")
    best_v: np.ndarray | None = None
    for active in combinations(range(len(constraints)), N + 1):
        matrix = []
        rhs = []
        for index in active:
            _, _, row, value = constraints[index]
            matrix.append(row)
            rhs.append(value)
        try:
            candidate = np.linalg.solve(np.vstack(matrix), np.array(rhs))
        except np.linalg.LinAlgError:
            continue

        v = candidate[:N]
        s = candidate[N]
        values = k_matrix @ v
        if np.any(v > 1 + 1e-8) or np.any(v < -1 - 1e-8):
            continue
        if np.min(values) < s - 1e-8 or np.min(values) < -1e-8:
            continue
        if s > best_s:
            best_s = float(s)
            best_v = v

    if best_v is None:
        raise RuntimeError("No feasible cone direction found.")
    return best_s, best_v, k_matrix @ best_v


def fmt(value: Fraction) -> str:
    return f"{float(value):.12f}"


def main() -> None:
    expected = (Fraction(34), Fraction(0), Fraction(-34), Fraction(0))
    for start, target in ((A, B), (B, A)):
        for block in range(N):
            context = partial_state(start, target, block)
            coeff = half_derivative_coefficients(block, context)
            assert coeff == expected

    print("Squared anchor distances sigma_i:")
    print([fmt(v) for v in SIGMA])
    print("Squared free-free distances beta_ij:")
    for i, j in PAIRS:
        print(f"  ({i + 1},{j + 1}) {fmt(beta(i, j))}")
    print()

    print(f"F(A) = {fmt(objective(A))}")
    print(f"F(B) = {fmt(objective(B))}")
    print("grad F(A) =", [fmt(v) for v in gradient(A)])
    print("grad F(B) =", [fmt(v) for v in gradient(B)])
    print()

    jacobian, preferences = branch_jacobian_and_preferences(A, B)
    jacobian_b, preferences_b = branch_jacobian_and_preferences(B, A)
    assert np.allclose(jacobian, jacobian_b)
    assert np.allclose(preferences, -preferences_b)

    monodromy = jacobian @ jacobian
    rho = max(abs(np.linalg.eigvals(monodromy)))
    print(f"spectral radius of two-sweep monodromy = {rho:.12f}")
    print("one-sweep branch Jacobian:")
    print(jacobian)

    cone_matrix = np.vstack([-preferences, preferences @ jacobian, -preferences @ monodromy])
    margin, direction, cone_values = solve_cone_margin(cone_matrix)
    print(f"cone LP margin = {margin:.12f}")
    print("cone direction v =", direction)
    print("cone values [-R_A v, R_A J v, -R_A J^2 v]:")
    print(cone_values.reshape(3, N))

    assert rho < 1
    assert margin > 0


if __name__ == "__main__":
    main()

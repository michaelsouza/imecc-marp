#!/usr/bin/env python3
"""Verify a Powell-type cycle for exact cyclic BCD on a 1D DGP instance.

The instance has one fixed anchor x0 = 0 and three free scalar vertices.
The graph is complete on vertices {0,1,2,3}.  Squared distances are

    d_{0i}^2 = 9/5,      i = 1,2,3,
    d_{ij}^2 = 18/5,     1 <= i < j <= 3.

With the cyclic order x1, x2, x3, the two states

    A = (-1,  1, -1),    B = ( 1, -1,  1)

form a legal two-cycle for the unregularized exact coordinate method if the
tie at each coordinate subproblem is resolved by choosing the sign in the
next state.  The points A and B are not stationary for the anchored problem.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


ANCHOR = Fraction(0)
A = (Fraction(-1), Fraction(1), Fraction(-1))
B = (Fraction(1), Fraction(-1), Fraction(1))

EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
D2 = {
    (0, 1): Fraction(9, 5),
    (0, 2): Fraction(9, 5),
    (0, 3): Fraction(9, 5),
    (1, 2): Fraction(18, 5),
    (1, 3): Fraction(18, 5),
    (2, 3): Fraction(18, 5),
}


def key(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def full_state(free: tuple[Fraction, Fraction, Fraction]) -> tuple[Fraction, ...]:
    return (ANCHOR, *free)


def objective(free: tuple[Fraction, Fraction, Fraction]) -> Fraction:
    x = full_state(free)
    total = Fraction(0)
    for i, j in EDGES:
        residual = (x[i] - x[j]) ** 2 - D2[key(i, j)]
        total += Fraction(1, 2) * residual**2
    return total


def gradient(free: tuple[Fraction, Fraction, Fraction]) -> tuple[Fraction, ...]:
    x = full_state(free)
    grad = [Fraction(0) for _ in x]
    for i, j in EDGES:
        diff = x[i] - x[j]
        residual = diff**2 - D2[key(i, j)]
        contribution = 2 * residual * diff
        grad[i] += contribution
        grad[j] -= contribution
    return tuple(grad[1:])


def local_value(t: Fraction, points: tuple[Fraction, ...], block: int) -> Fraction:
    total = Fraction(0)
    for i, j in EDGES:
        if block not in (i, j):
            continue
        other = j if i == block else i
        residual = (t - points[other]) ** 2 - D2[key(i, j)]
        total += Fraction(1, 2) * residual**2
    return total


def derivative_factor(points: tuple[Fraction, ...], block: int) -> tuple[Fraction, Fraction]:
    """Return c such that local derivative is c * t * (t^2 - 1).

    In every subproblem visited by the cycle, the other two free coordinates
    are +1 and -1.  The derivative simplifies exactly to 6 t (t^2 - 1).
    """
    coeff_t3 = Fraction(0)
    coeff_t = Fraction(0)
    coeff_t2 = Fraction(0)
    coeff_0 = Fraction(0)
    for i, j in EDGES:
        if block not in (i, j):
            continue
        other = j if i == block else i
        a = points[other]
        d2 = D2[key(i, j)]
        coeff_t3 += 2
        coeff_t2 += -6 * a
        coeff_t += 2 * (3 * a**2 - d2)
        coeff_0 += 2 * (-a**3 + d2 * a)
    assert coeff_t2 == 0
    assert coeff_0 == 0
    assert coeff_t == -coeff_t3
    return coeff_t3, coeff_t


def verify_transition(
    start: tuple[Fraction, Fraction, Fraction],
    target: tuple[Fraction, Fraction, Fraction],
) -> list[tuple[int, Fraction, Fraction, Fraction, Fraction]]:
    current = list(start)
    rows = []
    for block in range(1, 4):
        points = full_state(tuple(current))
        desired = target[block - 1]
        coeff_t3, coeff_t = derivative_factor(points, block)
        values = {t: local_value(t, points, block) for t in (Fraction(-1), Fraction(0), Fraction(1))}
        assert coeff_t3 == 6
        assert coeff_t == -6
        assert values[desired] == min(values.values())
        assert values[Fraction(-1)] == values[Fraction(1)]
        assert values[Fraction(0)] > values[desired]
        rows.append((block, desired, values[Fraction(-1)], values[Fraction(0)], values[Fraction(1)]))
        current[block - 1] = desired
    return rows


def fmt(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    print("Powell-type 1D DGP cycle")
    print("Edges:", EDGES)
    print("Squared distances:", {edge: fmt(D2[edge]) for edge in EDGES})
    for name, state in (("A", A), ("B", B)):
        print(f"{name} = {tuple(fmt(v) for v in state)}")
        print(f"  F({name}) = {fmt(objective(state))}")
        print(f"  grad F({name}) = {tuple(fmt(v) for v in gradient(state))}")

    for label, start, target in (("A -> B", A, B), ("B -> A", B, A)):
        print(label)
        for block, desired, v_minus, v_zero, v_plus in verify_transition(start, target):
            print(
                "  update x{} to {:>2}: local values phi(-1)={}, phi(0)={}, phi(1)={}".format(
                    block,
                    fmt(desired),
                    fmt(v_minus),
                    fmt(v_zero),
                    fmt(v_plus),
                )
            )

    for i, j in combinations(range(1, 4), 2):
        assert D2[(i, j)] == Fraction(2) * D2[(0, 1)]


if __name__ == "__main__":
    main()

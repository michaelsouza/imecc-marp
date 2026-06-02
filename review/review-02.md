Yes. A Powell-type **nonstationary cycle is possible** for anchored DGP, at least in the selection-dependent sense of exact BCD: if a vertex subproblem has several exact global minimizers, the algorithm can choose among them in a way that cycles. So a blanket impossibility theorem is false without extra assumptions such as unique vertex minimizers or a tie-breaking rule.

The anchored DGP objective in your note is

[
F(x)=\frac12\sum_{{i,j}\in E}
\left(|x_i-x_j|^2-d_{ij}^2\right)^2,
]

with one anchor fixed and cyclic exact minimization over free vertices. 

## A concrete anchored-DGP cycle in (K=1)

Take one anchor (x_0=0) and three free scalar vertices

[
x=(x_1,x_2,x_3)\in \mathbb R^3.
]

Use the complete graph on the three free vertices, plus edges from each free vertex to the anchor. Prescribe squared distances

[
d_{0i}^2=7,\qquad i=1,2,3,
]

and

[
d_{ij}^2=1,\qquad 1\le i<j\le 3.
]

Then the anchored DGP objective is

[
F(x)=
\frac12\sum_{i=1}^3 (x_i^2-7)^2
+
\frac12\sum_{1\le i<j\le 3}
\left((x_i-x_j)^2-1\right)^2 .
]

Use cyclic exact vertex updates in the order

[
x_1,\ x_2,\ x_3,\ x_1,\ x_2,\ x_3,\ldots
]

Now consider the six points

[
\begin{aligned}
p^0&=(-1,1,-1),\
p^1&=(1,1,-1),\
p^2&=(1,-1,-1),\
p^3&=(1,-1,1),\
p^4&=(-1,-1,1),\
p^5&=(-1,1,1),\
p^6&=p^0.
\end{aligned}
]

At (p^0), when updating (x_1), the other two free vertices are (1) and (-1). The one-variable subproblem is

[
\begin{aligned}
\phi(u)
&=
\frac12(u^2-7)^2
+\frac12((u-1)^2-1)^2
+\frac12((u+1)^2-1)^2\
&=
\frac32(u^2-1)^2+23.
\end{aligned}
]

Therefore the exact global minimizers are

[
u=-1 \quad\text{and}\quad u=1.
]

So moving from (x_1=-1) to (x_1=1) is an exact BCD update. The same calculation applies at every subsequent step, because the two non-updated active vertices are always (+1) and (-1). Thus, with the tie choice that flips the updated sign, exact cyclic vertex BCD follows

[
(-1,1,-1)
\to
(1,1,-1)
\to
(1,-1,-1)
\to
(1,-1,1)
\to
(-1,-1,1)
\to
(-1,1,1)
\to
(-1,1,-1).
]

The objective is constant on this cycle:

[
F(p^k)=\frac{127}{2},\qquad k=0,\ldots,5.
]

But the full anchored gradient is not small. For example,

[
\nabla F(-1,1,-1)=(0,12,0),
]

and

[
\nabla F(1,1,-1)=(0,0,-12).
]

Around the cycle, the nonzero gradient component simply rotates among the three free vertices, with sign changes. Hence

[
|\nabla F(p^k)|=12,\qquad k=0,\ldots,5.
]

So this is an exact anchored-DGP cyclic BCD orbit with

[
\nabla_zF \not\to 0.
]

## What this example does and does not prove

This proves that **unconditional impossibility is false**. Exact cyclic vertex BCD for anchored DGP can cycle at nonstationary points.

However, this is weaker than Powell’s strongest phenomenon. Powell’s examples show cycling behavior for cyclic exact coordinate minimization, including stable cycling under perturbations in one construction and smooth cycling in another; your note frames exactly that obstruction as the target phenomenon.  The DGP example above is more degenerate:

It relies on **ties** in the exact vertex subproblem.

The objective does **not strictly decrease** along the cycle.

A tie-breaking rule such as “stay put if the current block is already a minimizer” would not produce this particular orbit.

So the right conclusion is:

[
\boxed{
\text{A Powell-type obstruction exists for anchored DGP in the selection-dependent sense.}
}
]

But the stronger question remains:

[
\boxed{
\text{Can anchored DGP exhibit a robust Powell cycle with unique updates or strict descent?}
}
]

The example above does not settle that stronger version.

## Extension to any dimension (K)

The (K=1) example is already a valid anchored DGP example. If you want the same construction in (K\ge 2), one can add “guard” vertices to prevent off-axis minimizers.

Let the active vertices cycle on the line spanned by (e_1), taking values (\pm e_1). Add guard vertices

[
g_\ell^+=h e_\ell,\qquad g_\ell^-=-h e_\ell,
\qquad \ell=2,\ldots,K,
]

with (h>1). Connect every active vertex to every guard vertex with squared distance

[
1+h^2.
]

Also connect the guards to the anchor and to each other using the squared distances induced by their displayed positions, so the guard framework has zero residual at those positions.

When an active vertex (u\in\mathbb R^K) is updated while its two active neighbors are (+e_1) and (-e_1), write

[
u=t e_1+w,\qquad w\perp e_1,\qquad S=|u|^2.
]

The anchor and active-neighbor terms give

[
\phi_0(u)-\phi_0(e_1)
=====================

\frac32(S-1)^2-4|w|^2.
]

Each guard pair contributes

[
(S-1)^2+4h^2 w_\ell^2.
]

Summing over (\ell=2,\ldots,K),

[
\phi(u)-\phi(e_1)
=================

\left(K+\frac12\right)(S-1)^2
+
4(h^2-1)|w|^2.
]

Since (h>1), this is nonnegative and vanishes only when

[
w=0,\qquad S=1.
]

Thus the only active-vertex global minimizers are

[
u=e_1
\quad\text{and}\quad
u=-e_1.
]

The same six-cycle therefore persists in any dimension (K), with guard updates chosen to remain at their zero-residual positions.

## Conditional impossibility theorem

Although unconditional impossibility is false, there is a clean conditional theorem.

Suppose a cyclic exact vertex-BCD sequence has a compact level set and suppose that, on its cluster set, every vertex subproblem has a **unique** global minimizer. Then Powell-type cycling is impossible: every cluster point is stationary.

Proof sketch:

Since exact BCD is monotone,

[
F(z^{k+1})\le F(z^k),
]

so (F(z^k)) converges.

If some block step failed to vanish, then along a convergent subsequence we would get two distinct limiting values of the same block, with all other blocks fixed, giving the same limiting objective value. By continuity, both would be global minimizers of the same limiting vertex subproblem. That contradicts uniqueness.

Hence all block steps vanish.

Now take any cluster point (\bar z). Because each updated block is an exact minimizer and the intervening steps vanish, every block of (\bar z) minimizes its own vertex subproblem with all other vertices fixed.

For a differentiable objective, blockwise exact minimality implies

[
\nabla_{x_i}F(\bar z)=0
]

for every free vertex (i). Therefore

[
\nabla_zF(\bar z)=0.
]

This kind of conclusion is consistent with standard BCD convergence theory: Tseng proves stationarity of subsequential limits under compactness plus additional block-structure assumptions such as uniqueness of coordinate minimizers in enough blocks, and Grippo–Sciandrone prove related convergence results for block nonlinear Gauss–Seidel methods under two-block or quasiconvexity/pseudoconvexity-type hypotheses. ([Springer Nature Link][1])

So the sharp message is:

[
\boxed{
\text{No general impossibility theorem is true.}
}
]

But a useful theorem is:

[
\boxed{
\text{Powell-type nonstationary cycling requires limiting nonunique vertex minimizers.}
}
]

The explicit DGP example above realizes exactly that mechanism.

[1]: https://link.springer.com/article/10.1023/A%3A1017501703105 "https://link.springer.com/article/10.1023/A%3A1017501703105"

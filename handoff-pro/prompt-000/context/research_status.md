# Current Research Status

## Objective

We study exact cyclic block-coordinate descent (BCD) for the Euclidean distance-geometry least-squares objective

\[
F(x)=\frac12\sum_{\{i,j\}\in E}\left(\|x_i-x_j\|^2-d_{ij}^2\right)^2.
\]

The proximal variant in the paper has convergence guarantees. The unresolved question is how strong a Powell-type cycling pathology can be obtained for the unregularized exact cyclic method.

## Confirmed Fragile Counterexample

A valid one-dimensional DGP counterexample is already in the paper. It has one anchor \(x_1=0\), three free scalar variables \(z=(x_2,x_3,x_4)\), complete graph on four vertices, and squared distances

\[
d_{1i}^2=\frac95,\qquad i=2,3,4,
\]

\[
d_{ij}^2=\frac{18}{5},\qquad 2\le i<j\le4.
\]

With cyclic order \(x_2,x_3,x_4\), the states

\[
A=(-1,1,-1),\qquad B=(1,-1,1)
\]

form a legal two-cycle if global minimizer ties are resolved along the cycle. The gradients are

\[
\nabla F(A)=(0,8/5,0),\qquad \nabla F(B)=(0,-8/5,0),
\]

so the accumulation points are not stationary.

This example is fragile: it depends on exact ties. Perturbations under a natural global-minimizer rule typically leave the cycle and converge elsewhere.

## General Obstruction Already Identified

For any \(C^1\) objective bounded below, exact cyclic BCD with unique global block minimizers in a neighborhood of the omega-limit set cannot have a nonstationary cycle or nonstationary omega-limit set. Descent forces the objective to be constant along the limit set; uniqueness then implies each block is already a global minimizer, and differentiability gives block-gradient zero for all blocks.

Therefore any nonstationary exact-BCD cycling example must involve nonunique coordinate minimizers at the limiting points.

## Failed Stronger Candidate

The file `review/review-01.md` describes a stronger candidate in \(K=1\), with five free variables, complete free-free graph \(K_5\), and 30 anchors at zero for each free variable. It proposes a cycle

\[
A=(1,-1,1,-1,1),\qquad B=-A.
\]

The candidate has several correct formal properties:

- The ten coordinate subproblems visited by the proposed cycle satisfy
  \[
  \frac12\phi_i'(t)=34t(t^2-1),
  \]
  so \(t=\pm1\) are tied global minimizers at the limiting partial states.
- The gradients at \(A\) and \(B\) are nonzero.
- The intended branch map has a strongly contracting two-sweep linearization, with \(\rho(J^2)\approx0.01035725\).

However, this does **not** certify a basin of attraction for exact global BCD. The branch-selection cone claimed in `review-01.md` is not invariant under further iterates. For the displayed direction \(v\), the first tests pass, but the next odd test fails:

\[
R_AJ^3v\approx
(-2.998\cdot10^{-5},\;-1.030\cdot10^{-4},\;1.842\cdot10^{-4},\;6.496\cdot10^{-6},\;-3.172\cdot10^{-6}),
\]

which is not entrywise positive as required for the \(B\to A\) branch. Floating-point simulations with exact scalar global minimization also leave the proposed cycle after getting very close to it.

The lesson is that formal contraction of a selected branch is insufficient. A valid stronger example must prove that exact global minimization keeps selecting the desired branch for all iterates in a genuine open set.

## What We Need From Pro Extended

We need either:

1. A concrete DGP construction with a genuine open basin converging to a nonstationary tied cycle, including rigorous branch-selection/global-minimizer verification; or
2. A rigorous obstruction theorem for a relevant class, explaining why such a stable Powell-type DGP example cannot exist under those assumptions.

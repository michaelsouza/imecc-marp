# Pro Extended Handoff: Powell-Type Cycling for Exact BCD on DGP Objectives

Use extended reasoning internally, but return only a clear, verifiable mathematical answer: definitions, equations, proof sketches, construction parameters, computational checks, and failure modes. Do not provide private chain-of-thought.

## Core Question

We are studying exact cyclic block-coordinate descent (BCD) for the Euclidean distance-geometry least-squares objective

\[
F(x)=\frac12\sum_{\{i,j\}\in E}\left(\|x_i-x_j\|^2-d_{ij}^2\right)^2.
\]

Each BCD step updates one vertex/block by **global minimization** of the corresponding block subproblem. The proximal version converges under the assumptions in the attached paper, but the unregularized exact cyclic method can cycle because coordinate subproblems may have nonunique minimizers.

We already have a fragile DGP counterexample with nonstationary accumulation points. We want to know whether a stronger Powell-type example exists:

> Is there a finite DGP instance for which exact cyclic global BCD has a nonempty open set of initial conditions converging to a degenerate tied, nonstationary cycle?

If yes, please construct one explicitly and verify it. If no, please prove a rigorous obstruction for a clearly defined, relevant class.

Clarify the update rule in any proposed construction: a valid open-basin construction must either specify a deterministic global-minimizer selection rule and prove convergence for that single-valued map, or prove that along the proposed open basin every finite iterate has a unique global minimizer selecting the intended branch. Choosing favorable tied minimizers step-by-step is not sufficient.

## Required Deliverables

Please provide one of the following.

### Option A: A Concrete Stable/Open-Basin Example

Give a finite DGP instance with:

- dimension \(K\);
- graph \(G=(V,E)\);
- anchor coordinates, if any;
- squared distances \(d_{ij}^2\);
- cyclic block order;
- proposed cycle points;
- exact coordinate subproblem formulas at the cycle;
- proof that the cycle points are nonstationary;
- proof that an open set of initial conditions follows exact global minimizers and converges to the cycle.

Please state the modeling class explicitly: simple unweighted DGP is preferred; repeated coincident anchors and multiple anchors are allowed if essential; weighted DGP is acceptable only if you either explain how to replace weights by finite unweighted gadgets or clearly label the result as weighted. Squared distances should be nonnegative, and any zero-distance or coincident-anchor use should be stated explicitly.

The proof must address **global minimizer selection**, not only a locally selected branch. A contraction of a formal branch is not sufficient unless you also prove that the exact global BCD map keeps selecting that branch for all iterates in a genuine open set.

Useful checks include:

- exact polynomial identities for scalar subproblems;
- interval arithmetic or rational arithmetic for strict inequalities;
- an invariant cone, sector, or other invariant set for branch selection;
- spectral-radius computation for the return map, together with branch-selection invariance;
- verification that no lower scalar minimizer/root steals the coordinate update.

### Option B: A Rigorous Obstruction

Prove a theorem of the form:

\[
\text{No such stable/open-basin nonstationary cycle exists in class } \mathcal C.
\]

The class \(\mathcal C\) should be mathematically meaningful for DGP/BCD, for example:

- \(K=1\), one anchor, three free variables;
- scalar sign-flip cycles \(A\leftrightarrow -A\);
- DGP scalar quartic subproblems with two tied global minima at the limit;
- complete free-free graphs with separable anchored self-potentials;
- another clearly motivated finite-dimensional DGP subclass.

Please state all assumptions and prove the result. It is acceptable if the obstruction covers only a relevant subclass, as long as the boundaries of the theorem are clear.

## Known Results and Failed Attempts

### Fragile Valid Counterexample

We have a valid \(K=1\) DGP counterexample with one anchor \(x_1=0\), free variables \(z=(x_2,x_3,x_4)\), complete graph on four vertices, and

\[
d_{1i}^2=\frac95,\qquad d_{ij}^2=\frac{18}{5}.
\]

The two states

\[
A=(-1,1,-1),\qquad B=(1,-1,1)
\]

form a legal exact-BCD two-cycle under a tie selection rule. The gradients are

\[
\nabla F(A)=(0,8/5,0),\qquad \nabla F(B)=(0,-8/5,0),
\]

so the accumulation points are nonstationary. This example is fragile because it depends on exact ties and does not appear to have an open basin under natural global-minimizer selection.

### General Uniqueness Obstruction

For any \(C^1\) objective bounded below, exact cyclic BCD cannot have a nonstationary cycle or omega-limit set if global block minimizers are unique in a neighborhood of the omega-limit set. Descent forces constant objective on the limit set; uniqueness forces every block already to be a global minimizer; differentiability gives \(\nabla F=0\).

Therefore any nonstationary exact-BCD cycling example must involve nonunique coordinate minimizers at the limiting points.

### Failed K5 Candidate

A stronger \(K=1\), five-free-variable candidate was proposed in `context/review/review-01.md` and checked in `context/codes/powell_stable_dgp_candidate.py`.

It has:

- formal tied subproblems satisfying \(\frac12\phi_i'(t)=34t(t^2-1)\);
- nonzero gradients at the proposed cycle points;
- a contracting two-sweep branch linearization with \(\rho(J^2)\approx0.01035725\).

But it does **not** prove an open basin. The branch-selection cone is not invariant. See `context/branch_invariance_failure.md`. In particular, for the displayed direction \(v\), the test requiring \(R_AJ^3v>0\) fails.

Treat any open-basin, invariant-cone, or "certified" stability claim in `context/review/review-01.md` or `context/codes/powell_stable_dgp_candidate.py` as obsolete/false unless independently repaired. Use those files only for the candidate parameters and formal tied-subproblem/Jacobian calculations.

Please do not repeat this mistake: a formal branch contraction is not enough.

## Attached Context Files

- `context/research_status.md`: concise status summary.
- `context/branch_invariance_failure.md`: details of the failed branch-invariance check.
- `context/paper/proximal_bcd_graph_realization.tex`: paper draft with the fragile Powell-type counterexample and convergence framework.
- `context/paper/references.bib`: bibliography, including Powell 1973.
- `context/codes/powell_type_dgp_cycle.py`: exact rational verifier for the confirmed fragile example.
- `context/codes/powell_stable_dgp_candidate.py`: verifier for formal properties of the failed \(K_5\) candidate; note that it does not certify true branch invariance.
- `context/review/review-01.md`: write-up of the \(K_5\) candidate; useful but contains the overclaim about open basin.
- `context/references/powell1973search.md`: OCR/Markdown conversion of Powell's 1973 paper.

## Acceptance Criteria

A useful answer must include at least one of:

1. A complete construction with enough exact or interval-certified details that we can reproduce the DGP instance and verify global coordinate minimization; or
2. A theorem-quality obstruction with a precise class, assumptions, and proof.

Please also say which attached files you used and whether any claims depend on numerical evidence rather than proof.

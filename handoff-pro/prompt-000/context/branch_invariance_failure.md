# Branch-Invariance Failure in the K5 Candidate

The candidate in `review/review-01.md` proposes a \(K=1\) DGP instance with five free variables and a two-cycle

\[
A=(1,-1,1,-1,1),\qquad B=-A.
\]

It verifies tied coordinate subproblems at the cycle and obtains a contracting branch linearization. However, the open-basin claim is not established.

Let \(T\) be the intended one-sweep branch map. Near \(A\),

\[
T(A+\delta)=B+J\delta+O(\|\delta\|^2),
\]

and near \(B\),

\[
T(B+\eta)=A+J\eta+O(\|\eta\|^2).
\]

The claimed branch-selection inequalities are

\[
R_A\delta<0
\]

near \(A\), and

\[
R_B\eta<0,\qquad R_B=-R_A
\]

near \(B\).

The vector displayed in `review-01.md` is

\[
v=(-1,\;0.12887394,\;1,\;0.37927651,\;0.36789886).
\]

It satisfies the first checks:

\[
R_Av<0,\qquad R_BJv<0,\qquad R_AJ^2v<0.
\]

But the next odd branch-selection test fails. Since \(R_B=-R_A\), the \(B\to A\) condition at the next odd iterate requires

\[
R_AJ^3v>0.
\]

Numerically,

\[
R_AJ^3v\approx
(-2.99804879\cdot10^{-5},\;-1.02980268\cdot10^{-4},\;1.84160936\cdot10^{-4},\;6.49642319\cdot10^{-6},\;-3.17242783\cdot10^{-6}).
\]

This vector is not entrywise positive. Therefore the claimed cone is not invariant under the branch dynamics.

Consequently, the candidate should be treated as a formal tied-cycle and formal contraction example, not as a validated open-basin Powell-type example.

A valid stronger example must provide a genuinely invariant branch-selection region, or a different proof that the exact global minimizers continue to select the intended branch for all iterates in a nonempty open set.

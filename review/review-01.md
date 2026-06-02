# DGP Counterexample Candidate: Open Basin Converging to a Degenerate Nonstationary Cycle

## 1. Goal

We want a DGP instance for which exact cyclic BCD has the following behavior:

1. The limiting set is a nonstationary two-cycle.
2. At the limiting points, the active coordinate subproblems have global ties.
3. Away from the tie surface, the active coordinate minimizers are unique.
4. There is an open set of initial conditions whose iterates converge to the degenerate cycle.

The construction below targets exactly this mechanism. It does **not** require uniqueness at the limiting cycle, because uniqueness there is impossible for a nonstationary cycle under exact descent.

---

## 2. DGP Objective

Use five free scalar variables

\[
x=(x_1,x_2,x_3,x_4,x_5)\in\mathbb R^5.
\]

For each free variable \(x_i\), add \(30\) anchors fixed at coordinate \(0\), connected only to \(x_i\).

The equivalent DGP objective is

\[
F(x)
=
\frac12\sum_{i=1}^5\sum_{r=1}^{30}(x_i^2-\sigma_i)^2
+
\frac12\sum_{1\le i<j\le 5}
\bigl((x_i-x_j)^2-\beta_{ij}\bigr)^2.
\]

Equivalently,

\[
F(x)
=
15\sum_{i=1}^5(x_i^2-\sigma_i)^2
+
\frac12\sum_{1\le i<j\le 5}
\bigl((x_i-x_j)^2-\beta_{ij}\bigr)^2.
\]

The free-free graph is \(K_5\). The BCD order is

\[
x_1,x_2,x_3,x_4,x_5.
\]

The proposed two-cycle is

\[
A=(1,-1,1,-1,1),
\qquad
B=-A=(-1,1,-1,1,-1).
\]

---

## 3. Parameters

The anchor squared distances are

\[
\begin{aligned}
\sigma_1&=0.753758939998,\\
\sigma_2&=0.616449114477,\\
\sigma_3&=0.792699561907,\\
\sigma_4&=1.014161080698,\\
\sigma_5&=0.255417288626.
\end{aligned}
\]

The free-free squared distances are

\[
\begin{array}{c|c}
\text{edge} & \beta_{ij}=d_{ij}^2\\
\hline
12 & 11.667256816711\\
13 & 3.348281865181\\
14 & 0.026359083326\\
15 & 8.345334034857\\
23 & 0.286099935630\\
24 & 2.086006466135\\
25 & 13.467163347216\\
34 & 7.761224706209\\
35 & 10.823406635760\\
45 & 5.701577323401
\end{array}
\]

This is an explicit finite DGP instance in \(K=1\). It uses multiple anchors only to realize separable quartic terms of the form \(15(x_i^2-\sigma_i)^2\). In a weighted DGP formulation, the same effect could be represented more compactly.

---

## 4. Coordinate Subproblems

When updating \(x_i=t\), with the other coordinates fixed at \(c_j\), the scalar subproblem is

\[
\phi_i(t)
=
15(t^2-\sigma_i)^2
+
\frac12\sum_{j\ne i}
\bigl((t-c_j)^2-\beta_{ij}\bigr)^2.
\]

Its derivative satisfies

\[
\frac12\phi_i'(t)
=
30(t^2-\sigma_i)t
+
\sum_{j\ne i}
\bigl((t-c_j)^2-\beta_{ij}\bigr)(t-c_j).
\]

Consider the partial states along the sweep \(A\to B\). Before updating \(x_i\), the variables \(x_j\) with \(j<i\) have already been updated to \(B_j\), while variables \(x_j\) with \(j>i\) still have value \(A_j\).

For each of the five coordinate subproblems along \(A\to B\), direct substitution gives

\[
\boxed{
\frac12\phi_i'(t)=34t(t^2-1).
}
\]

The same identity holds for the five subproblems along \(B\to A\).

Therefore, at every partial state of the proposed cycle,

\[
\phi_i(t)
=
17(t^2-1)^2+C_i,
\]

for some constant \(C_i\) depending on the current context.

Hence the two global minimizers are exactly

\[
t=-1
\quad\text{and}\quad
t=1.
\]

The third critical point \(t=0\) is strictly worse, with objective gap \(17\). Thus, with appropriate tie choices, exact cyclic BCD has the legal two-cycle

\[
A\longleftrightarrow B.
\]

---

## 5. Nonstationarity of the Cycle

For the full objective, the gradient at \(A\) is approximately

\[
\nabla F(A)
\approx
(0,\;30.66902727,\;12.24872772,\;6.80630929,\;0).
\]

By symmetry,

\[
\nabla F(B)=-\nabla F(A).
\]

Thus

\[
\|\nabla F(A)\|
=
\|\nabla F(B)\|
\approx
33.71863595.
\]

Therefore the two-cycle is nonstationary, with gradient bounded well away from zero.

---

## 6. Local Branch Map

Let \(T\) denote one full cyclic BCD sweep.

Near \(A\), if the updates follow the branch \(A\to B\), then

\[
T(A+\delta)
=
B+J\delta+O(\|\delta\|^2).
\]

Near \(B\), the corresponding branch satisfies

\[
T(B+\eta)
=
A+J\eta+O(\|\eta\|^2).
\]

For the parameters above,

\[
J=
\begin{pmatrix}
0 & -0.17157731 & 0.12723115 & -0.00038763 & 0.05374509\\
0 & -0.00083958 & -0.00358477 & 0.14579213 & -0.19778353\\
0 & 0.00830374 & -0.00688232 & -0.08900195 & -0.01941437\\
0 & -0.02966847 & 0.02208420 & -0.01008859 & -0.06952601\\
0 & 0.01700533 & -0.01239620 & 0.01013377 & -0.00567816
\end{pmatrix}.
\]

The two-sweep monodromy is

\[
M=J^2.
\]

Its spectral radius is

\[
\boxed{
\rho(M)\approx 0.01035725<1.
}
\]

Thus the intended branch cycle is strongly linearly attracting.

---

## 7. Branch-Selection Inequalities

The crucial question is whether the intended branch is selected by global minimization for an open set of initial perturbations.

For perturbations \(x=A+\delta\), the first-order branch-preference inequalities for the five updates along \(A\to B\) are

\[
R_A\delta<0,
\]

where

\[
R_A=
\begin{pmatrix}
0 & -30.66902727 & 2.60687254 & 15.89456367 & -17.38133614\\
0 & -5.26210908 & -10.95354468 & -7.66786248 & 39.51696295\\
0 & -0.45975258 & 0.27842142 & -12.88007980 & -30.09171301\\
0 & 2.85850315 & -2.09838240 & -2.44904477 & 7.17419300\\
0 & 2.98933000 & -2.03816533 & -3.01636193 & 7.55873844
\end{pmatrix}.
\]

For perturbations \(x=B+\eta\), the branch-preference inequalities for \(B\to A\) are

\[
R_B\eta<0,
\qquad
R_B=-R_A.
\]

A concrete interior direction for the cone at \(A\) is

\[
v=
(-1,\;0.12887394,\;1,\;0.37927651,\;0.36789886).
\]

It satisfies

\[
R_Av
\approx
(-1.71170505,\;-0.00168784,\;-15.73664731,\;-0.01948354,\;-0.01610255)<0.
\]

After one branch sweep,

\[
R_BJv
\approx
(-0.49698401,\;-0.37453553,\;-0.39454115,\;-0.00419634,\;-0.00168784)<0.
\]

After two sweeps,

\[
R_AJ^2v
\approx
(-0.00168784,\;-0.01165420,\;-0.00820005,\;-0.00168782,\;-0.00168784)<0.
\]

Therefore the cone is mapped strictly into the opposite cone and then back into itself. Since

\[
\rho(J^2)<1,
\]

there is an open conic neighborhood of \(A\) such that all sufficiently small perturbations in that cone follow unique global coordinate minimizers and converge alternately to \(A\) and \(B\).

This gives the desired mechanism:

\[
\boxed{
\text{unique global minimizers off the tie surface, but convergence to a degenerate tied nonstationary two-cycle.}
}
\]

---

## 8. Why the Three-Variable Scalar Example Cannot Do This

The obstruction in the original \(K=1\), three-free-variable setting is not merely due to the symmetric \(\pm1\) choice.

Even if one allows arbitrary separable anchored self-potentials

\[
R_i(x_i)
\]

added to each coordinate, a period-two sign-flip cycle

\[
A=(s_1,s_2,s_3),
\qquad
B=-A,
\qquad
s_i\in\{\pm1\},
\]

forces equality of the incident free-free squared distances.

For example, for \(x_1\), the two contexts swap the signs of \(x_2\) and \(x_3\). Requiring both \(t=1\) and \(t=-1\) to be stationary tied candidates in both contexts gives conditions of the form

\[
Q_1(1)+2(4-\beta_{13})=0,
\qquad
Q_1(1)+2(4-\beta_{12})=0,
\]

where

\[
Q_1(t)=\frac12R_1'(t).
\]

Hence

\[
\beta_{12}=\beta_{13}.
\]

Doing the same for the other variables forces

\[
\beta_{12}=\beta_{13}=\beta_{23}.
\]

Then the branch-preference inequalities collapse to the same cyclic contradiction observed in the rectangular \(K_4\) setting.

Thus, in the scalar three-variable sign-flip class, arbitrary anchors do not rescue the construction.

The \(K_5\) construction above escapes this by using four free neighbors per variable. The stationarity and value-tie constraints cancel collectively rather than edge-by-edge.

---

## 9. Status and Recommended Use

This should be treated as a strong computationally verifiable counterexample candidate.

All checks reduce to polynomial identities and strict inequalities:

1. The ten tied coordinate subproblems satisfy

   \[
   \frac12\phi_i'(t)=34t(t^2-1).
   \]

2. The tied minimizers are global, with gap \(17\) to the third critical point.

3. The gradients at the two limit points are nonzero.

4. The branch monodromy satisfies

   \[
   \rho(J^2)\approx 0.01035725<1.
   \]

5. The strict cone inequalities certify an open basin of perturbations.

For a paper-ready version, the decimal parameters should be rationalized or certified by interval arithmetic. The construction itself is finite, explicit, and already in DGP form.

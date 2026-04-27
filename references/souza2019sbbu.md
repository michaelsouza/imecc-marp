

Advances and New Challenges on Branch-and-Prune Algorithm

Michael Souza
Carlile Lavor
Luiz Mariano Carvalho

DIMACS’19

Workshop on Optimization
in Distance Geometry

---

"Proteins are the machines and building blocks of living cells. If we compare a living body to our world, each cell corresponds to a town, and the proteins are the houses, bridges, cars, cranes, roads, airplanes, etc. There are huge numbers of different proteins, each one performing its specific task." [Neumaier97]

![img-0.jpeg](img-0.jpeg)

![img-1.jpeg](img-1.jpeg)

![img-2.jpeg](img-2.jpeg)

![img-3.jpeg](img-3.jpeg)

![img-4.jpeg](img-4.jpeg)

![img-5.jpeg](img-5.jpeg)

Geometry

Function

---

PDB Data Distribution by Experimental Method

|  X-Ray | 128,043  |
| --- | --- |
|  NMR | 11,109  |
|  Electron Microscopy | 2,416  |
|  Other | 257  |
|  Multi Method | 137  |
|  Total | 141,962  |

![img-6.jpeg](img-6.jpeg)

NMR does not give atoms location, but info about pairwise distances.

![img-7.jpeg](img-7.jpeg)

---

# Distance Geometry Problem (DGPk)

Given a simple weighted undirected graph  $G = (V, E, d)$ , where  $d: E \to \mathbb{R}_+$  find an embedding  $x: V \to \mathbb{R}^k$  such that

$$
\forall \{i, j \} \in E, | | x (i) - x (j) | | = d (i, j)
$$

If the vertices are ordered, we can represent a  $\mathrm{DGP}_{\mathrm{k}}$  instance using an Euclidean Distance Matrix (EDM)

![img-8.jpeg](img-8.jpeg)

![img-9.jpeg](img-9.jpeg)

---

7/45

- DGPₖ is a NP-Hard problem
&gt; Reduction to Subset-Sum [Saxe79]

Background

- MDGP can be reformulated as a global optimization problem. However, the number of local minima grows exponentially with the number of vertices.

---

8/45

# Background

- DGPₖ is a NP-Hard problem
&gt; Reduction to Subset-Sum [Saxe79]

- MDGP can be reformulated as a global optimization problem. However, the number of local minima grows exponentially with the number of vertices.

- As usually happens with NP-hard problems, not all instances of DGPₖ are in fact hard to solve.

![img-10.jpeg](img-10.jpeg)

When all distances are available:
[Blumenthal53] O(|V|³) operations using SVD;

---

9/45

- DGPₖ is a NP-Hard problem
&gt; Reduction to Subset-Sum [Saxe79]

## Background

- MDGP can be reformulated as a global optimization problem. However, the number of local minima grows exponentially with the number of vertices.
- As usually happens with NP-hard problems, not all instances of DGPₖ are in fact hard to solve.

![img-11.jpeg](img-11.jpeg)

## When all distances are available:

[Blumenthal53] O(|V|³) operations using SVD;
[Dong02] O(|V|) operations by setting one point at a time following a specific order.

---

10/45

☑ DGPₖ is a NP-Hard problem
&gt; Reduction to Subset-Sum [Saxe79]

Background

## Background

## MDGP can be reformulated as a global optimization problem. However the number of local minima grows exponentially with

## The power of Divide-and-conquer

Instead of solving the entire problem at one step taking O( (kn)³ ) operations., it's better to solve n small problems with O( k³ ) operations.

![img-12.jpeg](img-12.jpeg)

When all distances are available:

[Blumenthal53] O( |V|³ ) operations using SVD;

[Dong02] O( |V| ) operations by setting one point at a time following a specific order.

---

# Discretizable Molecular  $\mathrm{DGP}_{\mathrm{k}}$  (DMDGPk)

We say that a  $\mathrm{DGP_k}$  instance is discretizable molecular DGP if there is an order  $v_{1}, v_{2}, \ldots, v_{n}$  of  $V$  satisfying the following requirements:

1.  $E$  contains all cliques on  $k + 1$  consecutive vertices;
2. For all  $w \in \{ k + 1, \dots, n \}$ ,

$$
\operatorname {Vol} (\{d (i, j): w - k \leq i, j &lt;   w \}) &gt; 0
$$

(2D:  $x_{w-1}$  and  $x_{w-2}$  are not overlapped)

(3D:  $x_{w-1}, x_{w-2}, x_{w-3}$  are not colinear)

![img-13.jpeg](img-13.jpeg)
The other distances (gray) may be available

---

✓ DMDGPₖ is NP-hard (theoretical challenge).

[Lavor12]

12/45

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge)
$\checkmark$  The search space of  $\mathsf{DMDGP}_k$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )

![img-14.jpeg](img-14.jpeg)

![img-15.jpeg](img-15.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
$\checkmark$  The search space of  $\mathsf{DMDGP}_k$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )

![img-16.jpeg](img-16.jpeg)

![img-17.jpeg](img-17.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
$\checkmark$  The search space of  $\mathsf{DMDGP}_k$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )

![img-18.jpeg](img-18.jpeg)

![img-19.jpeg](img-19.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$ is NP-hard (theoretical challenge).
$\checkmark$ The search space of  $\mathsf{DMDGP}_k$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )

![img-20.jpeg](img-20.jpeg)

![img-21.jpeg](img-21.jpeg)

[Lavor12]

---

$\nu$  DMDGP $_k$  is NP-hard (theoretical challenge).
The search space of  $\mathsf{DMDGP}_k$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )

![img-22.jpeg](img-22.jpeg)

![img-23.jpeg](img-23.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
$\checkmark$  The search space of  $\mathsf{DMDGP}_k$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )

![img-24.jpeg](img-24.jpeg)

![img-25.jpeg](img-25.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$ is NP-hard (theoretical challenge).
$\checkmark$ The search space of DMDGP $_k$ is finite (discrete), but it can be large (up to $2^{|V| - k}$)

![img-26.jpeg](img-26.jpeg)

![img-27.jpeg](img-27.jpeg)

[Lavor12]

---

$\nu$ DMDGP $_k$ is NP-hard (theoretical challenge).
$\checkmark$ The search space of DMDGP $_k$ is finite (discrete), but it can be large (up to $2^{|V| - k}$)

## Globally not rigid

![img-28.jpeg](img-28.jpeg)

[Lavor12]

![img-29.jpeg](img-29.jpeg)

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
The search space of  $\mathrm{DMDGP_k}$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )
$\checkmark$  Almost surely, any distance  $d(i,j)$  implicitly defines the distances  $d(u,v) \forall u,v \in \{i+1,\dots,j-1\}$

![img-30.jpeg](img-30.jpeg)

![img-31.jpeg](img-31.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
The search space of  $\mathrm{DMDGP_k}$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )
$\checkmark$  Almost surely, any distance  $d(i,j)$  implicitly defines the distances  $d(u,v) \forall u,v \in \{i+1,\dots,j-1\}$

![img-32.jpeg](img-32.jpeg)

![img-33.jpeg](img-33.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
The search space of  $\mathrm{DMDGP_k}$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )
$\checkmark$  Almost surely, any distance  $d(i,j)$  implicitly defines the distances  $d(u,v) \forall u,v \in \{i+1,\dots,j-1\}$

![img-34.jpeg](img-34.jpeg)

![img-35.jpeg](img-35.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
The search space of  $\mathrm{DMDGP_k}$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )
$\checkmark$  Almost surely, any distance  $d(i,j)$  implicitly defines the distances  $d(u,v) \forall u,v \in \{i+1,\dots,j-1\}$

![img-36.jpeg](img-36.jpeg)

![img-37.jpeg](img-37.jpeg)

[Lavor12]

---

$\checkmark$ DMDGP $_k$  is NP-hard (theoretical challenge).
The search space of  $\mathrm{DMDGP_k}$  is finite (discrete), but it can be large (up to  $2^{|V| - k}$ )
$\checkmark$  Almost surely, any distance  $d(i,j)$  implicitly defines the distances  $d(u,v) \forall u,v \in \{i+1,\dots,j-1\}$

![img-38.jpeg](img-38.jpeg)

![img-39.jpeg](img-39.jpeg)

[Lavor12]

---

# Independent Components (IC) – Symmetries

![img-40.jpeg](img-40.jpeg)

The orientation of the structures associated to the distances  $d(i,j)$  and  $d(u,v)$ , where  $i &lt; j$ ,  $u &lt; v$  are independent iff

$$
j - u \leq k
$$

![img-41.jpeg](img-41.jpeg)

A set with more than  $k$  points in  $\mathbb{R}^k$  has internal orientation (3D chirality)

26/45

---

# Independent Components (IC) – Symmetries

![img-42.jpeg](img-42.jpeg)

The number of solutions of  $\mathrm{DMDGP}_k$  is  $2^{\lambda}$ , where  $\lambda$  is the number of IC's.

A set with more than  $k$  points in  $\mathbb{R}^k$  has internal orientation (3D chirality)

![img-43.jpeg](img-43.jpeg)

---

# Does  $\mathrm{DMDGP}_{\mathrm{k}}$  have any application?

Hand-crafted re-order of protein backbone

![img-44.jpeg](img-44.jpeg)

# Assumptions (knowledge):

$\checkmark$  Distances between close hydrogen atoms are known.

---



DMDGP_{k} re-order

Given a undirected graph G=( V, E ), a DMDGP_{k} re-order is a surjetive function r: { 1, …, m } → V, with length m in IN bounded by polynomial in |V|, such that

1. { r( 1 ), …, r( k ) } defines a k-clique in G;

2. For i in { k + 1, …, m }, { r(i – k + j ), r( i ) } in E for 0 < j < k.

3. For i in { k + 1, …, m }, { r(i – k ), r( i ) } in E or r(i – k) = r( i )

If there is a re-order, then the discrete search space and symmetry properties are valid

---

# How hard is to check if there is a re-order (k fixed)?

![img-45.jpeg](img-45.jpeg)

![img-46.jpeg](img-46.jpeg)

1. Creates G' from G.
2. Does G' have a component H whose the union of its vertices is V?

Decision problem

$$
O \left(k ^ {3} | V | ^ {k + 1}\right) \text{ steps}
$$

---

# How hard is to check if there is a re-order (k fixed)?

![img-47.jpeg](img-47.jpeg)
(f,g,a,c,e,c,b,d)

![img-48.jpeg](img-48.jpeg)

1. Creates G' from G.
2. Does G' have a component H whose the union of its vertices is V?
3. Creates a re-order from any walk on H.

Finding a re-order

$$
O \left(| V | ^ {2 k}\right) \text{ steps}
$$

---

![img-49.jpeg](img-49.jpeg)
Branch-and-Prune (BP) algorithm

---

# Branch-and-Prune (BP) algorithm

1. Explores all branches (one at time)
2. Prunes the infeasible ones ASAP
3. When a one solution is found, BP generates all the other by flipping the IC’s.

33/45

---

# Improving BP (SBP)

BP remembers its mistakes (prunning), but not its right choices.

![img-50.jpeg](img-50.jpeg)

- Amount of calculations (worst scenario)

Let say that solutions are

$$
C_{11} + C_{22} + C_{32} \text{ and } C_{12} + C_{21} + C_{31}
$$

BP will generate the following configurations

$$
\begin{array}{l}
\checkmark C_{11} + C_{21} + C_{31} \quad 2^{|C_1| - k} + 2^{|C_2| - k} + 2^{|C_3| - k} \\
\checkmark C_{11} + C_{21} + C_{32} \quad 2^{|C_1| - k} \\
\checkmark C_{11} + C_{22} + C_{31} \quad 2^{|C_2| - k} \\
\checkmark C_{11} + C_{22} + C_{32} \quad 2^{|C_3| - k} \\
\end{array}
\quad
2^{|C_1| - k} + 2^{|C_2| - k} + 2^{|C_3| - k + 1}
$$

This additional cost can be reduced

---

# Improving BP (SBP)

BP remembers its mistakes (prunning), but not its right choices.

![img-51.jpeg](img-51.jpeg)

- Deal with all components as rigid bodies (only reflections are needed)

Again, let say that solutions are

$$
C_{11} + C_{22} + C_{32} \text{ and } C_{12} + C_{21} + C_{31}
$$

SBP will generate the following configurations

- $C_{11} + C_{21} + C_{31}$ $2^{|C_1| - k} + 2^{|C_2| - k} + 2^{|C_3| - k}$
- $C_{11} + C_{21} + C_{32}$ $|C_3|$
- $C_{11} + C_{22} + C_{31}$ $|C_2|$
- $C_{11} + C_{22} + C_{32}$ $C_3$

![img-52.jpeg](img-52.jpeg)

![img-53.jpeg](img-53.jpeg)

---

# Improving BP (SBP)

BP remembers its mistakes (prunning), but not its right choices.

![img-54.jpeg](img-54.jpeg)

- Deal with all components as rigid bodies (only reflections are needed)

Again, let say that solutions are

$$
C_{11} + C_{22} + C_{32} \text{ and } C_{12} + C_{21} + C_{31}
$$

SBP will generate the following configurations

- $C_{11} + C_{21} + C_{31}$ $2^{|C_1| - k} + 2^{|C_2| - k} + 2^{|C_3| - k}$
- $C_{11} + C_{21} + C_{32}$ $|C_1| - |C_2| - |C_3|$
- $C_{11} + C_{21} + C_{31}$ $C_2$ $|C_2| + 2|C_3|$
- $C_{11} + C_{21} + C_{32}$ $C_3$ This additional cost can be reduced

---

# Improving BP (SBP)

BP remembers its mistakes (prunning), but not its right choices.

![img-55.jpeg](img-55.jpeg)

- Deal with all components as a rigid bodies (only reflections are needed)
- Only the possible locations of the last node need to be tested

37/45

---

# Improving BP (SBP)

BP remembers its mistakes (prunning), but not its right choices.

![img-56.jpeg](img-56.jpeg)

- Deal with all components as a rigid bodies (only reflections are needed)
- Only the possible locations of the last node need to be tested

38/45

---

# Improving BP (SBP)

BP remembers its mistakes (prunning), but not its right choices.

![img-57.jpeg](img-57.jpeg)

- Deal with all components as a rigid bodies (only reflections are needed)
- Only the possible locations of the last node need to be tested

$$
C_{11} + C_{21} + C_{31} \quad 2^{|C_1| - k} + 2^{|C_2| - k} + 2^{|C_3| - k}
$$

$$
C_{11} + C_{21} + C_{32} \quad 1 \quad 1
$$

$$
C_{11} + C_{22} + C_{31} \quad 1 \quad 1
$$

$$
C_{11} + C_{22} + C_{32} \quad 1
$$

$$
|C_2| + 2|C_3|
$$

39/45

---

# Improving BP (SBP)

BP remembers its mistakes (prunning), but not its right choices.

![img-58.jpeg](img-58.jpeg)

- Deal with all components as rigid bodies (only reflections are needed)
- Only the possible locations of the last node need to be tested
- Merge components

![img-59.jpeg](img-59.jpeg)

---

# Numerical Experiments

Implementations of SBP, BP, MDJEEP (V 2.0) [Mucherino19]

Language: C &amp; C++,

Compiler: GCC 7.4.0 (Ubuntu 7.4.0-1ubuntu1~18.04.1)

Intel(R) Core(TM) i3-7100U CPU @ 2.40GHz (16 GB RAM)

Only one solution

Default config of MDJEEP

Max time 500 secs

Instances:

- Artificial
- Only N, C, Ca
- {i, j} in E iff
((j - i) &lt; 4) or
d(i, j) &lt; 5 AA and atom(i) = atom(j) = C

41/45

---

# Numerical Experiments

![img-60.jpeg](img-60.jpeg)
data/1bdk.nmr, nnodes:15, nedges:40

![img-61.jpeg](img-61.jpeg)
data/1afo.nmr, nnodes:120, nedges:518

![img-62.jpeg](img-62.jpeg)
data/3fhh.nmr, nnodes:1863, nedges:7796

![img-63.jpeg](img-63.jpeg)

![img-64.jpeg](img-64.jpeg)

![img-65.jpeg](img-65.jpeg)

---

# Numerical Experiments

|   |   |   | TIME (secs) |   |   | LDE  |   |   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  fname | natoms | nedges | MDJEEP | BP | SBP | MDJEEP | BP | SBP  |
|  1bfw | 12 | 30 | 1.60E-05 | 2.60E-06 | 6.92E-07 | 2.23E-16 | 5.58E-16 | 6.09E-16  |
|  1bdk | 15 | 40 | 2.50E-05 | 3.51E-06 | 1.28E-06 | 2.53E-16 | 3.75E-16 | 4.07E-16  |
|  1f3r | 27 | 91 | 5.70E-05 | 1.58E-05 | 6.41E-06 | 2.46E-16 | 8.87E-16 | 8.06E-16  |
|  1bcv | 57 | 187 | 9.20E-05 | 2.03E-05 | 8.55E-06 | 3.75E-06 | 2.18E-13 | 2.96E-13  |
|  1afo | 120 | 518 | 1.99E-04 | 3.47E-05 | 2.42E-05 | 1.82E-06 | 2.23E-12 | 1.97E-12  |
|  1b4r | 240 | 1016 | 4.97E+00 | 6.80E-01 | 9.52E-04 | 5.79E-14 | 2.19E-12 | 3.39E-12  |
|  1a3s | 474 | 2043 | 1.84E+01 | 3.61E+00 | 1.74E-03 | 5.25E-07 | 9.28E-11 | 2.13E-11  |
|  1all | 480 | 2199 | 1.06E-01 | 1.28E-02 | 1.72E-04 | 1.98E-07 | 9.02E-13 | 2.45E-12  |
|  1bob | 918 | 3858 | 1.48E+02 | 7.36E+00 | 1.46E-02 | 1.48E-13 | 1.86E-11 | 2.38E-11  |
|  3fhh | 1863 | 7796 | -- | -- | 2.66E-01 | -1.00E+00 | -1.00E+00 | 6.50E-11  |

---

# References

[Saxe79] Saxe, J.B.: Embeddability of weighted graphs in k-space is strongly NP-hard. In: Proceedings of 17th Allerton Conference in Communications, Control and Computing, pp. 480–489 (1979)

[Neumaier97] Neumaier, A. (1997). Molecular modeling of proteins and mathematical prediction of protein structure. SIAM review, 39(3), 407–460.

[Liberti08] Liberti, L., Lavor, C., &amp; Maculan, N. (2008). A branch-and-prune algorithm for the molecular distance geometry problem. International Transactions in Operational Research, 15(1), 1–17.

[Lavor12] Lavor, C., Liberti, L., Maculan, N., &amp; Mucherino, A. (2012). The discretizable molecular distance geometry problem. Computational Optimization and Applications, 52(1), 115–146.

44/45

---

45/45

# References

[Mucherino12] Mucherino, A., Lavor, C., &amp; Liberti, L. (2012). Exploiting symmetry properties of the discretizable molecular distance geometry problem. Journal of Bioinformatics and Computational Biology, 10(03), 1242009.

[Lavor13] Lavor, C., Liberti, L., &amp; Mucherino, A. (2013). The interval Branch-and-Prune algorithm for the discretizable molecular distance geometry problem with inexact distances. Journal of Global Optimization, 56(3), 855-871.

[Cassioli15] Cassioli, A., Günlük, O., Lavor, C., &amp; Liberti, L. (2015). Discretization vertex orders in distance geometry. Discrete Applied Mathematics, 197, 27-41.

[Lavor19] Lavor, C., Souza, M., Carvalho, L. M., &amp; Liberti, L. (2019). On the polynomiality of finding $^{\kappa}$DMDGP re-orders, to appear.
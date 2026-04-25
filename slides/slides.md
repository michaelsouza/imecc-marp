---
marp: true
theme: default
style: |
  @import url('./theme.css');
paginate: true
math: katex
title: Continuous and Discrete Optimization Applied to Distance Geometry
---

<!-- _class: title -->
<!-- _paginate: false -->

# Continuous and Discrete Optimization Applied to Distance Geometry

Michael Souza  
IMECC/UNICAMP

<div class="logos">
  <img src="./images/logo-imecc_site-header.png" alt="IMECC">
  <img src="./images/logo-uec.png" alt="UNICAMP">
</div>

---

# Core Idea

Distance Geometry transforms partial metric information into geometric configurations.

<div class="cols">
<div class="panel">

### Input

A weighted graph

$$G = (V,E,d)$$

with distances known only for a subset of edges.

</div>
<div class="panel">

### Output

A realization

$$x: V \to \mathbb{R}^k$$

that satisfies the prescribed distances.

</div>
</div>

---

# Distance Geometry Problem

Given a simple weighted undirected graph $G=(V,E,d)$, find an embedding $x:V\to\mathbb{R}^k$ such that

$$
\forall \{i,j\}\in E,\quad \|x(i)-x(j)\| = d(i,j).
$$

When the data are uncertain, one may consider interval constraints:

$$
l_{ij}\leq \|x_i-x_j\|\leq u_{ij},\quad (i,j)\in E.
$$

<p class="note">This formulation naturally arises in experimental settings, such as NMR data.</p>

---

# Why It Matters

<div class="cols-3">
<div class="panel">

### Molecules

NMR does not provide atomic coordinates directly; it provides distance constraints between selected atom pairs.

</div>
<div class="panel">

### Sensors

Localization in networks without GPS can be formulated from relative distances.

</div>
<div class="panel">

### Data

MDS and dimensionality reduction use distances to construct geometric representations.

</div>
</div>

---

# Computational Complexity

- The DGP in $\mathbb{R}^k$ is NP-hard.
- Direct global optimization formulations typically exhibit many local minima.
- Nevertheless, instances with additional structure can often be solved much more efficiently.

> The practical question becomes: which geometric structure can be exploited?

---

# Two Complementary Strategies

<div class="cols">
<div class="panel">

### Continuous Optimization

Formulate the problem with real-valued variables and minimize distance violations.

$$
\min_x \sum_{(i,j)\in E}
\big(\|x_i-x_j\|^2-d_{ij}^2\big)^2
$$

</div>
<div class="panel">

### Discrete Optimization

Use geometric assumptions to transform continuous choices into finitely many alternatives.

$$
\text{tree search}
$$

</div>
</div>

---

# Continuous Formulation

A common approach is to minimize an error function:

$$
f(x)=
\sum_{(i,j)\in E}
\left(\|x_i-x_j\|^2-d_{ij}^{\,2}\right)^2.
$$

<div class="cols">
<div>

**Strengths**

- flexible under noisy data;
- compatible with classical optimization methods;
- suitable for interval constraints.

</div>
<div>

**Challenges**

- many variables: $kn$;
- nonconvexity;
- the number of local minima grows rapidly with $n$.

</div>
</div>

---

# Complete Distance Information

If the Euclidean distance matrix is complete, the geometry becomes substantially more rigid.

- Blumenthal: solution by SVD in $O(|V|^3)$.
- Dong and Wu: point-by-point construction in a suitable order, with linear cost in structured cases.

<p class="note">The key point is that vertex order and distance structure may radically change the computational cost.</p>

---

# Geometric build-up

Incremental construction places each new point from previously positioned points.

In $\mathbb{R}^3$, if a new atom has exact distances to three previously positioned atoms, its location is restricted to a small number of possibilities.

$$
(i,i_k)\in E,\quad k=1,2,3
$$

This observation is the starting point for discretization.

---

# DMDGP

A DGP instance is discretizable when there exists an order

$$
v_1,v_2,\dots,v_n
$$

such that each new vertex has known distances to the previous $k$ vertices, and these $k$ vertices form a nondegenerate simplex.

<div class="math-large">

$$
\operatorname{Vol}
\left(\{d(i,j): w-k\leq i,j<w\}\right)>0
$$

</div>

---

# Geometric Consequence

In the discretizable case:

- the search space is no longer continuous;
- each new point generically has two possible positions;
- the search can be represented by a binary tree;
- the maximum size is $2^{|V|-k}$.

<p class="note">The problem remains NP-hard, but its structure becomes algorithmically exploitable.</p>

---

# Branch-and-Prune

The Branch-and-Prune algorithm combines geometric construction with pruning by additional distances.

1. Position the first $k$ vertices.
2. For each new vertex, generate the possible positions.
3. Test all available distance constraints.
4. Eliminate infeasible branches as early as possible.

---

# Pruning

An additional distance $(i,j)\in E$ may eliminate choices before the search tree grows too large.

$$
\left|\|x_i-x_j\|-d_{ij}\right|\leq \varepsilon
$$

<div class="cols">
<div class="panel">

### Branch

Explore the geometric alternatives compatible with the $k$ predecessors.

</div>
<div class="panel">

### Prune

Use additional edges to discard incompatible configurations.

</div>
</div>

---

# Symmetries and Independent Components

In the DMDGP, successive reflections generate solutions related by symmetry.

A distance $d(i,j)$ can, almost surely, determine internal distances among intermediate vertices:

$$
d(u,v),\quad u,v\in\{i+1,\dots,j-1\}.
$$

The number of solutions can be described as

$$
2^\lambda,
$$

where $\lambda$ is the number of independent components.

---

# Improving BP: SBP

BP records pruning failures, but it does not necessarily reuse successful structural choices.

Sequential Branch-and-Prune treats components as rigid bodies:

- reflect entire components;
- avoid recomputing already determined substructures;
- test only the relevant positions of the last node;
- reduce the additional cost associated with symmetric solutions.

---

# Bridge to Dihedral Angles

For proteins, some assumptions are more natural:

- covalent bond lengths are known;
- planar angles are known;
- distances between nearby hydrogen atoms can be estimated by NMR.

The main degree of freedom then becomes the dihedral angle.

---

# Dihedral-Angle Parametrization

Each new point can be expressed through a rotation:

$$
x_k = R\left(\omega_k,x_k^0,x_{k-1},x_{k-2}\right),
$$

with

$$
\omega_k\in[0,2\pi].
$$

This parametrization reduces the number of variables from $3n$ to $n$ in $\mathbb{R}^3$.

---

# DGP in Angle Variables

Find $\omega\in[0,2\pi]^n$ such that

$$
l_{ij}
\leq
\|x_i(\omega)-x_j(\omega)\|
\leq
u_{ij},
\quad (i,j)\in E.
$$

<div class="cols">
<div class="panel">

### Continuous

Search over a hypercube of angles.

</div>
<div class="panel">

### Geometric

Preserve molecular structure within the model itself.

</div>
</div>

---

# Continuous and Discrete in One Problem

<div class="cols">
<div>

**Continuous**

- suitable for noise;
- allows relaxations;
- supports gradient, quasi-Newton, and global search methods.

</div>
<div>

**Discrete**

- exploits local rigidity;
- enumerates solutions;
- uses pruning by additional constraints;
- reveals symmetries.

</div>
</div>

---

# A Computational Narrative

1. Use molecular knowledge to order or reorder the vertices.
2. Identify when the instance is discretizable.
3. Apply Branch-and-Prune when the structure permits it.
4. Use a continuous dihedral-angle parametrization when discretization is too restrictive.
5. Combine geometric pruning, local optimization, and uncertainty handling.

---

# References

- Saxe, J. B. Embeddability of weighted graphs in $k$-space is strongly NP-hard, 1979.
- Blumenthal, L. M. Theory and Applications of Distance Geometry, 1953.
- Dong, Q.; Wu, Z. A linear-time algorithm for solving the MDGP with exact inter-atomic distances, 2002.
- Liberti, L. et al. A branch-and-prune algorithm for the MDGP, 2008.
- Lavor, C. et al. Artificial backbones of hydrogens, 2011.
- Souza, M.; Lavor, C.; Carvalho, L. M. DIMACS 2019 materials.

---

<!-- _class: title -->
<!-- _paginate: false -->

# Thank You

Questions?

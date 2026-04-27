Algorithmica (2021) 83:2400-2426 https://doi.org/10.1007/s00453-021-00835-6

Check for updates

# A New Algorithm for the  $^{K}$ DMDGP Subclass of Distance Geometry Problems with Exact Distances

Douglas S. Gonçalves $^{1}$ $\bullet$  Carlile Lavor $^{2}$ $\cdot$  Leo Liberti $^{3}$ $\cdot$  Michael Souza $^{4}$

Received: 12 October 2020 / Accepted: 13 May 2021 / Published online: 23 May 2021
© The Author(s), under exclusive licence to Springer Science+Business Media, LLC, part of Springer Nature 2021

# Abstract

The fundamental inverse problem in distance geometry is the one of finding positions from inter-point distances. The Discretizable Molecular Distance Geometry Problem (DMDGP) is a subclass of the Distance Geometry Problem (DGP) whose search space can be discretized and represented by a binary tree, which can be explored by a Branch-and-Prune (BP) algorithm. It turns out that this combinatorial search space possesses many interesting symmetry properties that were studied in the last decade. In this paper, we present a new algorithm for this subclass of the DGP, which exploits DMDGP symmetries more effectively than its predecessors. Computational results show that the speedup, with respect to the classic BP algorithm, is considerable for sparse DMDGP instances related to protein conformation.

Keywords Distance Geometry  $\cdot$  Discretization  $\cdot$  Symmetries  $\cdot$  DMDGP

This is an extended version of an abstract published at CTW2020.

Douglas S. Gonçalves
douglas.goncalves@ufsc.br

Carlile Lavor
clavor@unicamp.br

Leo Liberti
liberti@lix.polytechnique.fr

Michael Souza
michael@ufc.br

$^{1}$  CFM, Federal University of Santa Catarina, 88.040-900 Florianópolis, Brazil
2 IMECC, University of Campinas, 13081-970 Campinas, Brazil
3 LIX CNRS École Polytechnique, Institut Polytechnique de Paris, 91128 Palaiseau, France
4 Federal University of Ceará, 60440-900 Fortaleza, Brazil

Springer

---

# Introduction

Given a simple, undirected, weighted graph G = (V, E, d), with weight function d : E → (0, ∞) and an integer K > 0, the Distance Geometry Problem (DGP) [26] consists in finding a realization x : V → ℝ^{K} such that,$$\forall \left\{u,v\right\}\in E,\;{||x_{u}-x_{v}||}=d_{uv},$$where || · || denotes the Euclidean norm, x_{v} := x(v),∀v ∈ V and d_{uv} := d({u,v}),∀{u,v} ∈ E. Each equation in (1) is called a distance constraint. We say that a realization x satisfies d_{uv} if the corresponding distance constraint is verified. A realization satisfying all distance constraints in (1) is called a valid realization. We shall call a pair (G, K) a DGP instance.

There are many applications of Distance Geometry, mainly related to K ∈ {1,2,3} [4, 5, 34]. An application to Data Science can be found in [23], and a very recent survey on this subject is [22]. An important class of the DGP arises in the context of 3D protein structure calculations (K = 3), with distance information provided by Nuclear Magnetic Resonance (NMR) experiments [9, 31, 39].

Existence and uniqueness of DGP solutions, among other theoretical aspects of the problem, are discussed in [26]. Henceforth, we will consider that the DGP admits a solution.

### Assumption 1

The solution set of (1) is non-empty.

The DGP is naturally cast as a search in continuous space. Depending on the graph structure, however, combinatorial search algorithms can be defined, notably via the identification of appropriate vertex orders [7, 15, 21]. Although DGP is NP-hard [38], these combinatorial approaches allowed to show that it is Fixed Parameter Tractable (FPT) on certain graph structures, as those arising in protein conformation [27].

The aforementioned vertex orders define a DGP subclass, called the Discretizable Molecular Distance Geometry Problem (DMDGP) [16, 17], formally given as follows.

### Definition 1

A DGP instance (G, K) is a^{K}DMDGP if there is a vertex order v_{1}, ..., v_{n} ∈ V, such thatG[{v_{1}, ..., v_{K}}] is a clique;G(a) For every i > K, v_{i} is adjacent to v_{i-1}, , ..., v_{i-K},CM(v_{i-1}, ..., v_{i-K}) ≠ 0.

In the above definition, G[·] denotes the induced subgraph and CM(v_{i-1}, ..., v_{i-K}) is the Cayley-Menger determinant of v_{i-1}, ..., v_{i-K} [26, Sec. 2]. Its squared value is proportional to the (K - 1)-volume of a realization x_{i-1}, ... ,x_{i-K} for v_{i-1}, ..., v_{i-K}.

---

Algorithmica (2021) 83:2400-2426

Fig.1 Example of a  $^2$  DMDGP graph. The natural vertex order (1,2,3,4,5,6,7) satisfies the assumptions of Definition 1 for  $K = 2$ . In this case, the discretization edges are given by  $E_{D} = \{\{1,2\}, \{1,3\}, \{2,3\}, \{2,4\}, \{3,4\}, \{3,5\}, \{4,5\}, \{4,6\}, \{5,6\},$  and the pruning edges by  $E_{P} = \{\{1,5\}, \{1,7\}, \{4,7\}\}$

![img-0.jpeg](img-0.jpeg)

Condition  $CM(\nu_{i-1}, \dots, \nu_{i-K}) \neq 0$  means that the points  $x_{i-1}, \dots, x_{i-K}$  span an affine subspace of dimension  $K - 1$ .

Although Definition 1 applies to any dimension  $K$ , therefore covering other applications rather than molecular conformation where  $K = 3$ , the term "molecular" is commonly kept in the related literature [7, 16, 26], regardless of the dimension, to enforce the property that the adjacent predecessors of  $\nu_{i}$  are contiguous (the term "contiguous  $K$ -lateration order" to mean  $^K$  DMDGP is used in [7]), a desirable property when ordering atoms of a protein [15, 16].

When the dimension  $K$  is clear from the context, we shall simply use DMDGP rather than  $^K$  DMDGP. Moreover, without loss of generality, whenever we denote an edge by  $\{v_i, v_j\} \in E$ , we will assume that  $i &lt; j$ , i.e  $v_i$  precedes  $v_j$  in the vertex order of Definition 1.

Example 1 Consider a DGP instance, with  $K = 2$ , represented by the graph  $G = (V, E)$  of Fig. 1, and suppose  $d: E \to (0, \infty)$  is such that Assumption 1 holds. It is simple to verify that the vertex order  $\nu_1 = 1, \nu_2 = 2, \ldots, \nu_7 = 7$  satisfies the assumptions of Definition 1. Hence, this example is a  $^2$ DMDGP.

Properties 1 and 2(a) of Definition 1 says that  $G$  is composed by a chain of contiguous  $(K + 1)$ -cliques. Moreover, properties 1 and 2 allow us to turn the search space into a binary tree, in the following way. After fixing the positions for the first  $K$  vertices, for each new vertex  $\nu_{i}$ , with  $i &gt; K$ , property 2(a) ensures that the possible positions  $x_{i}$  for vertex  $\nu_{i}$  lie in the intersection of  $K$  spheres centered at  $x_{i - 1},\ldots ,x_{i - K}$  with radii  $d_{i - 1,i},\ldots ,d_{i - K,i}$ , respectively. Property 2(b) guarantees that there are at most two points, let us say  $\{x_i^+,x_i^-\}$ , in such intersection [30]. This sphere intersection can be computed in many different ways that we will not cover in this paper but are well studied in the literature [1, 6, 30].

Remark 1 The above process is known in the literature as  $K$ -lateration [26].

Thus, following the vertex order, after fixing the first  $K$  vertices, each new vertex has at most 2 possible positions, which of course depend on the position of its  $K$  immediate adjacent predecessors, leading to a binary tree of possible positions,

Springer

---

Algorithmica (2021) 83:2400-2426

where each path, from the root to a leaf node, corresponds to a possible realization for the graph  $G$ .

However, not all of these possible realizations (paths on the tree) are valid, because  $G$  may contain other edges  $\{h,i\}$ , with  $|h - i| &gt; K$ , associated to distance constraints that are not satisfied by such realizations. The edges given in Definition 1 are called discretization edges and the others, that may be (or not) available, are called pruning edges [16]. In Example 1, the discretization edges for the vertex order (1,2,3,4,5,6,7) are shown in Fig. 1 with solid lines whereas pruning edges are represented by dashed lines.

Henceforth, let us partition  $E = E_{D} \cup E_{P}$ , where  $E_{D}$  is the set of discretization edges and  $E_{P}$  the set of pruning edges. Clearly, we can also partition the equations in (1) in discretization edge constraints and pruning edge constraints. We remark that, according to Definition 1,  $E_{D} = \{\{i,j\} \in E \mid |i - j| \leq K\}$  and therefore  $E_{P} = \{\{i,j\} \in E \mid |i - j| &gt; K\}$ .

The Branch-and-Prune (BP) algorithm [25] explores the DMDGP binary tree in a depth-first manner and validates possible positions for vertices as soon as a pruning edge appears. A pseudo-code is given in Algorithm 1.

In Algorithm 1, the phrase "  $x_{i}^{+}$  is feasible" means that the equations

$$
\forall h: h &lt;   i, \left\{v _ {h}, v _ {i} \right\} \in E _ {P}, \quad \left\| x _ {i} ^ {+} - x _ {h} \right\| = d _ {h i},
$$

are satisfied up to a certain tolerance  $\varepsilon &gt; 0$ . In Step 5 positions  $x_{i}^{+}$  and  $x_{i}^{-}$  are computed via  $K$ -lateration. See [1, 26, 30] for details.

Algorithm 1 BP
1: BP(i,n,G,x) # (i&gt;K)
2: if (i&gt;n) then
3: return x
4: else
5: Find solutions  $\{x_i^+, x_i^-\}$  for the system:  $\| x_\ell - x_i \|^2 = d_{\ell,i}^2, \ell = i - K, \dots, i - 1$ .
6: if  $x_i^+$  is feasible then
7: Set  $x_i = x_i^+$  and call BP(i+1,n,G,x). # 1st candidate position
8: end if
9: if  $x_i^-$  is feasible then
10: Set  $x_i = x_i^-$  and call BP(i+1,n,G,x). # 2nd candidate position
11: end if
12: end if

Computational experiments in [16] showed that BP outperforms methods based on global continuation [32] and semidefinite programming [14] on instances of the DMDGP subclass, suggesting BP as the method of choice for this subclass of DGPs.

We remark that if the requirements of Definition 1 are strengthened, replacing  $K$  by  $K + 1$  (except for the dimension), then we have a  $(K + 1)$ -lateration order and the  $^K$ DMDGP instance can be solved in linear time by a Geometric Build-up algorithm [11].

In addition to the discretization of the DGP search space, the DMDGP order also implies symmetry properties of such discrete space [19, 24, 28, 37]. From the

Springer

---

computational point of view, one of the most important of such properties, in the context of this paper, is that all DMDGP solutions can be determined from just one solution. This property is related to the DMDGP symmetry vertices, which can be identified a priori, based on the input graph (see next section). Once a first solution is found, the others can be obtained by partial reflections of the first, based on symmetry hyperplanes associated to these vertices.

Previous works [27, 33] exploited symmetry to reconstruct all valid realizations from the first one found and to prove that the BP algorithm is fixed-parameter tractable. Others [12, 13], considered decomposition-based variants of BP which leverage DMDGP symmetry information.

In this work, we exploit DMDGP symmetry in order to find the first valid realization more quickly. We handle the DMDGP as a sequence of nested subproblems, each one defined by a pruning edge {i,j} ∈ E_{P}. For each subproblem, we can exploit any realization x (valid or not) for building the symmetry hyperplanes (which will define partial reflections). Once we have them, we apply compositions of such partial reflections only to x_{j} to find its correct position. Only after finding the correct combination of partial reflections do we use it to obtain the positions of other vertices. After a subproblem is solved, the set of valid partial reflections is reduced and a single symmetry hyperplane is enough to handle positions x_{i+K}, ... ,x_{j} in the next subproblem.

In terms of the system of nonlinear equations (1), we solve a subset of equations and then gradually include new equations to this subset: the new equations are solved subject to the original equations in the subset. This process is repeated until all equations in (1) are satisfied.

These ideas lead to a new algorithm which deals with pruning edges, one-by-one, and takes advantage of a valid realization for already solved subproblems. Computational results illustrate the advantage of the new algorithm, compared with the classic BP.

This paper is organized as follows. Section 2 briefly reviews the main results about DMDGP symmetries and Sect. 3 explains how they can be used to solve a sequence of nested subproblems. The new algorithm, its correctness and implementation details are presented in Sect. 4, and comparisons with the classic BP in protein-like instances are given in Sect. 5. Concluding remarks are given in Sect. 6.

## DMDGP Symmetries

Before discussing the new algorithm, we shall present a theoretical background on DMDGP symmetries and recall some results from the literature [26, 27, 29].

Given a realization x satisfying (1), it is clear that there are uncountably many others, which satisfy the same set of distances, and which can be obtained by translations, rotations or reflections of x (because these transformations preserve all pairwise distances). Since the assumptions of Definition 1 ensure that the first K vertices form a clique, a valid realization for G[v_{1}, ... ,v_{K}] in ℝ^{K} can be found by matrix decomposition methods [10] or a sequence of sphere intersections [1], for example. Once the positions of these first K vertices are fixed, the degrees of freedom of

---

translations and rotations are removed. We also remark that these first K positions define a hyperplane and, given a realization x, the total reflection of x through this hyperplane is also a solution of (1).

From here on, we say that two realizations are incongruent (modulo translations and rotations) if they are not translations, rotations or total reflections of each other. For technical reasons, we only allow the total reflections through the hyperplane defined by the positions of the first K vertices in the vertex order (so two realizations, one of which is a reflection of the other through this hyperplane, will both be considered members of any set of “incongruent realizations”).

###### Definition 2

Let $\hat{X}$ be the set of all incongruent realizations satisfying distance constraints associated to discretization edges in E_{D}, i.e. ∥x_{i} - x_{j}∥ = d_{ij} such that |i - j| ≤ K. A realization x ∈ $\hat{X}$ is called a possible realization.

As discussed in Sect. 1, each x ∈ $\hat{X}$ corresponds to a path from the root to a leaf node in the binary tree of a DMDGP instance. Notice that |$\hat{X}$| = 2^{|V|-K}.

###### Definition 3

A realization x ∈ $\hat{X}$ is said to be valid if x is a solution of (1). Let X ⊂ $\hat{X}$ denote the set of all incongruent valid realizations of a DMDGP instance.

The computational experiments in [16] suggested that |X| is always a power of 2. A conjecture was formulated and quickly disproved using some instances constructed by hand, until the conjecture was shown to be true with probability one in [29].

Given x ∈ $\hat{X}$, for i > K, let R^{i}_{x}(y) be the reflection of y ∈ ℝ^{K} through the hyperplane defined by x_{i-K}, ... ,x_{i-1}, with normal p_{i}:

R^{i}_{x}(y) = (I - 2p_{i}p_{i}^{T})(y - x_{i-1}) + x_{i-1} ,

assuming ∥p_{i}∥ = 1. Let us also define, for all i > K and x ∈ $\hat{X}$, partial reflection operators:

g_{i}(x) = (x_{1}, x_{2}, ... , x_{i-1}, R^{i}_{x}(x_{i}), R^{i}_{x}(x_{i+1}), ... , R^{i}_{x}(x_{n})).

###### Remark 2

Some direct but useful properties of reflections and partial reflections are in order:A reflection R^{i}_{x}(y) preserves the distance from y to any point in the hyperplane defined by x_{i-K}, ... ,x_{i-1}.The pairwise distances for R^{i}_{x}(x_{i}),R^{i}_{x}(x_{i+1}),... ,R^{i}_{x}(x_{n}) are the same as those for x_{i}, ... ,x_{n}. As a consequence of this, and the fact that R^{i}_{x}(x_{ϕ}) = x_{ϕ}, for ϕ = i - K, ... ,i - 1, all pairwise distances for x_{i-K}, ... ,x_{n} from x ∈ $\hat{X}$ are preserved in g_{i}(x).Partial reflections preserve distances related to discretization edges E_{D}, so that g_{i}(x) ∈ $\hat{X}$

---

Algorithmica (2021) 83:2400-2426

![img-1.jpeg](img-1.jpeg)
Fig.2 The leftmost path/realization  $x$  is represented by a straight line whereas the rightmost  $x'$  by a dashed line. All 4 possible positions for the fourth vertex (denoted by  $a, b, c$  and  $d$ ) can be generated by  $x$  and its induced reflections  $R_x^3$  and  $R_x^4$ . Illustration for  $K = 2$

4. All realizations in  $\hat{X}$  can be generated from a single  $x\in \hat{X}$  by the composition of partial reflection operators  $g_{i}$  [26, Sec. 3.3.8].

Let us now recall one of the main results about  $^K$  DMDGP symmetries.

Theorem 1 (Theorem 3.2 in [26]) With probability 1, for all  $j &gt; K + i$ , there is a set  $H^{ij}$  of  $2^{j - i - K}$  real positive values such that for each  $x \in \hat{X}$ , we have  $\| x_j - x_i\| \in H^{ij}$ . Furthermore, for all  $x', x \in \hat{X}$  such that  $x' \neq x$  and  $x_t' = x_t$ , for  $t \leq i + K - 1$ ,  $\| x_j - x_i\| = \| x_j' - x_i\|$  if and only if  $x_j' = R_x^{i + K}(x_j)$ .

In Theorem 1, "with probability 1" means that the set of  $^K$  DMDGP instances for which the statements do not hold has Lebesgue measure zero in the set of all  $^K$  DMDGP instances [29].

The first part of Theorem 1 says that, for  $j &gt; K + i$ , the possible realizations  $x \in \hat{X}$  yield a set of  $2^{j - i - K}$  distinct values for  $\| x_{i} - x_{j}\|$ . Let  $\hat{X}_{i + K - 1}(x)$  be the subset of possible realizations  $x' \in \hat{X}$  that agree with  $x$  in the first  $i + K - 1$  positions. Given a possible realization  $x' \in \hat{X}_{i + K - 1}(x)$ , each of these  $2^{j - i - K}$  distinct values is associated to a pair of  $2^{j - i - K + 1}$  possible positions for  $x_{j}$  from realizations in  $\hat{X}_{i + K - 1}(x)$  (see Fig. 2 where possible values for  $\| x_{1} - x_{3}\|$  and  $\| x_{1} - x_{4}\|$  are represented by the radii of gray and, respectively, black arcs centered at  $x_{1}$ ).

Since  $j - i &gt; K$ , if the distance  $d_{ij}$  is available, it must be a pruning distance. In view of Assumption 1, then  $\| x_i - x_j\| = d_{ij}$  for some  $x\in \hat{X}$ . Let such a  $x$  define the set  $\hat{X}_{i + K - 1}(x)$ . Now, from the second part of Theorem 1, we have that among the possible realizations  $x^{\prime}\in \hat{X}_{i + K - 1}(x)$ , only those such that  $x_{j}^{\prime}\in \{x_{j},R_{x}^{i + K}(x_{j})\}$  are feasible with respect to  $d_{ij}$ . If  $\nu_{j}$  is the last vertex in the order, then only two realizations in  $\hat{X}_{i + K - 1}(x)$  are feasible.

For every DMDGP solution, there is another one symmetric to the hyperplane defined by the positions of the first  $K$  vertices. Moreover, as a consequence of

Springer

---

Algorithmica (2021) 83:2400-2426

Theorem 1, the number of solutions doubles for every other symmetry vertex belonging to the following set [29]:

$$
S := \left\{v _ {\ell} \in V \mid \not \exists \left\{v _ {i}, v _ {j} \right\} \in E \text { with } i + K &lt;   \ell \leq j \right\}. \tag {3}
$$

The vertex  $\nu_{K + 1}$  is always in  $S$ , because the first  $K$  vertices define a symmetry hyperplane. The other symmetry hyperplanes are given by the positions of  $\nu_{i - K},\ldots ,\nu_{i - 1}$ , if  $\nu_{i}\in S$ , for  $i &gt; K + 1$ . As mentioned in the Sect. 1,  $S$  can be computed before solving a  $^K$ DMDGP instance, which implies that the number of solutions is known a priori, and given by  $2^{|S|}$ , with probability one.

Theorem 2 (Theorem 3.4 in [27]) Let  $(G, K)$  be a feasible  $^K$  DMDGP and  $S$  its set of symmetry vertices. Then, with probability  $1, |X| = 2^{|S|}$ .

The  $2^{|S|}$  valid realizations are incongruent modulo translations and rotations, meaning that they differ one from another only by partial reflections (or a total reflection through the first symmetry hyperplane, as explained above).

It is important to notice from (3) that the addition of new pruning edges in  $E$  may reduce the number of elements (symmetry vertices) in  $S$ . Recalling Example 1, since  $K = 2$  and  $\{1,7\} \in E$ , it follows that  $S = \{3\}$ . But, if edge  $\{1,7\}$  is removed, then  $S = \{3,4\}$ .

A direct consequence of Theorem 2 is the following corollary.

Corollary 1 Let  $(G, K)$  be a feasible  $^K$  DMDGP instance where  $|V(G)| = n &gt; K$ . If  $\{v_1, v_n\} \in E$ , then  $(G, K)$  has only two incongruent solutions which are reflections of each other through the symmetry hyperplane defined by the position of the first  $K$  vertices.

Proof If  $\{v_1, v_n\} \in E$ , then  $S = \{v_{K+1}\}$ , which implies that the number of solutions is  $2^{|S|} = 2^1 = 2$ . If one of these solutions is  $x$ , then the other is  $x'$ , the reflection of  $x$  through the hyperplane defined by  $x_1, \ldots, x_K$ .

A result that will be useful ahead is given in Proposition 1 and illustrated in Fig. 2.

Proposition 1 (Lemma 4.2 in [27]) Let  $x \in \hat{X}$ ,  $k &gt; i + 1$  and  $p_i \neq p_k$  be the normals to the hyperplanes defining  $R_{\downarrow}^{i}(\cdot)$  and  $R_{\downarrow}^{k}(\cdot)$ , respectively. If  $y \in \mathbb{R}^K$  is not in the hyperplanes containing the origin and normal to  $p_i$  and  $p_k$ , then

$$
R _ {g _ {i} (x)} ^ {k} \left(R _ {\downarrow} ^ {i} (y)\right) = R _ {\downarrow} ^ {i} \left(R _ {\downarrow} ^ {k} (y)\right).
$$

Proposition 1 tells us that compositions of partial reflections that depend on more than one realization (e.g.  $x$  and  $x' := g_k(x)$ ) can be described in terms of reflections based on a single realization. For example, for  $k &gt; i$ , we have

Springer

---

Algorithmica (2021) 83:2400-2426

$$
\begin{array}{l}
(g_{k} \circ g_{i})(x) := g_{k}(g_{i}(x)) \\
= g_{k}(x_{1}, \dots, x_{i-1}, R_{x}^{i}(x_{i}), \dots, R_{x}^{i}(x_{n})) \\
= (x_{1}, \dots, x_{i-1}, R_{x}^{i}(x_{i}), \dots, R_{x}^{i}(x_{k-1}), R_{x'}^{k}(R_{x}^{i}(x_{k})), \dots, R_{x'}^{k}(R_{x}^{i}(x_{n}))) \\
= (x_{1}, \dots, x_{i-1}, R_{x}^{i}(x_{i}), \dots, R_{x}^{i}(x_{k-1}), R_{x}^{i}(R_{x}^{k}(x_{k})), \dots, R_{x}^{i}(R_{x}^{k}(x_{n})), \\
\end{array}
$$

where the last equality follows from Proposition 1.

Therefore, for a DMDGP, given $x \in \hat{X}$, problem (1) can be cast as finding a binary vector $s \in \{0,1\}^{n - K}$, such that

$$
x(s) := U(x, s) = g_{K+1}^{s_{1}} \circ \dots \circ g_{n}^{s_{n-K}}(x) \tag{4}
$$

satisfies $\| x_{i}(s) - x_{j}(s)\| = d_{ij}$, for all $\{i,j\} \in E$. Here, $g_{i}^{1}(\cdot) = g_{i}(\cdot)$ and $g_{i}^{0}(\cdot) = I(\cdot)$, where $I(x) = x$. In Sect. 3 we shall explain how to efficiently perform the search of this binary vector taking into account DMDGP symmetry information.

To close this section, let us describe how to generate other valid realization $x(s') \in X$ from a given one $x(s) \in X$. Let $x(s)$ be a valid realization for $(G, K)$. The vertices in the set $S$ determine which components of the binary vector $s \in \{0, 1\}^{n-K}$ from (4) are allowed to change in order to obtain another valid realization for $(G, K)$. In other words, the search space for the new $s' \in \{0, 1\}^{n-K}$ is reduced to

$$
s' \in B := \{s' \in \{0, 1\}^{n-K} \mid s'_{\ell} = s_{\ell} \text{ if } v_{K+\ell} \notin S\}. \tag{5}
$$

Recall that for Example 1 without edge $\{1,7\}$, we have $S = \{3,4\}$. Suppose a valid realization $x$ corresponds to $s = (0,0,0,0,0)$. Then, the other three valid realizations correspond to $s' = (1,0,0,0,0)$, $s'' = (0,1,0,0,0)$ and $s''' = (1,1,0,0,0)$.

**Lemma 1** Let $x(s)$ be a valid realization for a feasible $^K$ DMDGP instance $(G, K)$. For every $s' \in B$, $x(s') \in X$.

**Proof** Since $x(s')$ from Eq. (4) involves only partial reflections, in view of Property 3 in Remark 2, $x(s') \in \hat{X}$, i.e. $\| x_i(s') - x_j(s') \| = d_{ij}, \forall \{i,j\} \in E_D$.

It remains to show that $x(s')$ does not violate distance constraints associated to pruning edges $\{i,j\} \in E_P$. Since the reflections are applied to positions $x_\ell$ such that $\ell \geq K + 1$, edges $\{i,j\} \in E$ with $i &lt; j \leq K$ are not affected. Thus, assume that $K + 1 \leq j \leq n$.

We have that $v_{i + K + 1}, \ldots, v_j \notin S$, and from (4) and (2), the positions $x_\ell, \ldots, x_{i + K + 1}, \ldots, x_j$ are updated by reflections $R_x^\ell(x_\ell), \ldots, R_x^\ell(x_{i + K + 1}), \ldots, R_x^\ell(x_j)$, for $K + 1 \leq \ell \leq i + K$ such that $v_\ell \in S$. Since either $i \leq \ell - 1$, i.e. $x_i$ is in the hyperplane associated to $v_\ell$, or $i \geq \ell$, i.e. $x_i$ comes after this hyperplane, in view of Remark 2, Property 2, these reflections are such that $\| x_i(s') - x_j(s') \| = d_{ij}$.

Springer

---

Algorithmica (2021) 83:2400-2426

# 3 Nested DMDGP Subproblems

Given a DMDGP instance, properties 1 and 2 of Definition 1 give rise to a rich symmetric structure for the corresponding DGP problem, as discussed in Sect. 2.

It is interesting to observe that, on one hand, the absence of pruning edges turns the DMDGP into a trivial problem, because any path from the root to a leaf node of the search tree corresponds to a valid realization, i.e  $\hat{X} = X$ , and all other solutions can be built by partial reflections. On the other hand, one of the most challenging DMDGP instances to solve with BP is the one where the only pruning edge is  $\{v_1, v_n\}$ . In that case, feasibility can only be verified at a leaf node, and for a standard depth-first search (DFS), it may represent a costly backtracking process until the first valid realization is found.

Differently, given  $x \in \hat{X}$ , the present proposal is to iteratively handle the pruning edge constraints following a given order  $&lt;$  on  $E_P$ .

As mentioned in Sect. 2 (after Theorem 2), each pruning edge  $\{i,j\}$  may reduce the set of valid partial reflection operations that can be applied to realizations of the vertices  $\nu_{i + K},\ldots ,\nu_{j}$ . Thus, by keeping track of valid partial reflections (or equivalently their corresponding symmetry vertices), it is possible to consistently modify a given realization satisfying a subset of distance constraints to also satisfy a new pruning edge constraint. This process is repeated until all distance constraints are satisfied.

For this, we enumerate edges in  $E_P$  as  $e_1, e_2, \ldots, e_m$ , with  $m = |E_P|$ , and use  $e_k &lt; e_\ell$  to mean that edge  $e_k$  precedes  $e_\ell$  in this order. We define the set of pruning edges preceding edge  $\{i,j\}$  by

$$
P ^ {i j} := \{\{u, w \} = e ^ {\prime} \in E _ {P} \mid e ^ {\prime} &lt;   e = \{i, j \} \}. \tag {6}
$$

In Example 1, we may consider the following order for the pruning edges

$$
e _ {1} = \{1, 5 \}, \quad e _ {2} = \{4, 7 \}, \quad e _ {3} = \{1, 7 \}, \tag {7}
$$

giving rise to

$$
P ^ {1 5} = \varnothing , \quad P ^ {4 7} = \{\{1, 5 \} \}, \quad P ^ {1 7} = \{\{1, 5 \}, \{4, 7 \} \}.
$$

Following the pruning edge order  $e_1, \ldots, e_m$ , we can define a sequence of subproblems spanned by  $e_k = \{i,j\} \in E_P$ .

Definition 4 Let  $(G, K)$  be a feasible  $^K DMDGP$  with  $G = (V, E, d)$ . Let  $G^{ij} = (V, E^{ij}, d_{|E^{ij}})$ , where  $E^{ij} = E_D \cup P^{ij} \cup \{i, j\}, \{i, j\} \in E_P$  and  $d_{|E^{ij}}$  is the restriction of  $d$  to  $E^{ij}$ . We say that  $(G^{ij}, K)$  is a  $^K DMDGP$  subproblem of  $(G, K)$  spanned by pruning edge  $\{i, j\}$ .

Springer

---

Algorithmica (2021) 83:2400-2426

It is clear that  $(G^{ij},K)$  is itself a DMDGP. Let us denote by  $X(G)$  the solution set of  $(G,K)$  and by  $X^{ij} := X(G^{ij})$  the solution set of  $(G^{ij},K)$ .

Proposition 2 Let  $G = (V, E, d)$  and  $H = (V, F, \hat{d})$  such that  $(G, K)$  and  $(H, K)$  are feasible  $^K$  DMDGPs. If  $E \subset F$  and  $d(\{i, j\}) = \hat{d}(\{i, j\}), \forall \{i, j\} \in E$ , then  $X(G) \supset X(H)$ .

Proof It follows directly from the fact that if  $x \in X(H)$ , then  $\| x_i - x_j \| = \hat{d}_{ij}, \forall \{i, j\} \in F$ , and since  $F \supset E$ , we have  $\| x_i - x_j \| = \hat{d}_{ij} = d_{ij}, \forall \{i, j\} \in E$ , which implies that  $x \in X(G)$ .

Let  $(G^{uw},K)$  and  $(G^{ij},K)$  be DMDGP subproblems spanned by edges  $\{u,w\}$  and  $\{i,j\}$ , respectively, such that  $\{u,w\} &lt; \{i,j\}$ . In view of Proposition 2, we have  $X^{uw} \supset X^{ij}$ .

Moreover, in this sequence of DMDGP subproblems, each time a new pruning edge is included, e.g.  $E^{ij} = E^{uw} \cup \{i,j\}$ , the set of symmetry vertices (see Eq. (3)) for  $(G^{uw},K)$  may be affected. This motivates us to define the set of necessary symmetry vertices for subproblem  $(G^{ij},K)$  as:

$$
S ^ {i j} = \left\{v _ {\ell} \in \left\{v _ {i + K + 1}, \dots , v _ {j} \right\} \mid \not \exists \{u, w \} \in P ^ {i j}, u + K &lt;   \ell \leq w \right\}. \tag {8}
$$

Remark 3 We remark that differently from  $S$ ,  $S^{ij}$  may be empty. For instance, if we consider Example 1 with the pruning edge order given by (7), then

$$
S ^ {1 5} = \{4, 5 \}, \quad S ^ {4 7} = \{7 \}, \quad S ^ {1 7} = \{6 \}.
$$

However, if the edges in  $E_P$  follow the order  $e_1 = \{1, 5\}$ ,  $e_2 = \{1, 7\}$  and  $e_3 = \{4, 7\}$ , then

$$
S ^ {1 5} = \{4, 5 \}, \quad S ^ {1 7} = \{6, 7 \}, \quad S ^ {4 7} = \varnothing .
$$

Furthermore, for  $e_1 = \{1,7\}$ ,  $e_2 = \{1,5\}$  and  $e_3 = \{4,7\}$ , we have

$$
S ^ {1 7} = \{4, 5, 6, 7 \}, \quad S ^ {1 5} = \varnothing , \quad S ^ {4 7} = \varnothing .
$$

The meaning of  $S^{ij} = \varnothing$  will become clear ahead (see Proposition 5).

Let  $x(s)$  be the current realization which is valid for  $(G^{uw}, K)$  and let  $e_{k+1} = \{i, j\} &gt; \{u, w\} = e_k$ . The vertices in the set  $S^{ij}$  determine which components of the binary vector  $s \in \{0, 1\}^{n-K}$  from Eq. (4) are allowed to change in order to obtain a valid realization for  $(G^{ij}, K)$ . In other words, the search space for the new  $s' \in \{0, 1\}^{n-K}$  is reduced to

$$
s ^ {\prime} \in B ^ {i j} := \left\{s ^ {\prime} \in \{0, 1 \} ^ {n - K} \mid s _ {\ell} ^ {\prime} = s _ {\ell} \text { if } v _ {i + K + \ell} \notin S ^ {i j} \right\}. \tag {9}
$$

Springer

---

Algorithmica (2021) 83:2400-2426

Lemma 2 Let  $S^{ij} \neq \emptyset$  and  $e_{k+1} = \{i,j\} &gt; \{u,w\} = e_k$ . Let  $x(s)$  be a valid realization for  $(G^{uw},K)$ . For every  $s' \in B^{ij}$ ,  $x(s') \in X^{uw}$ .

Proof The proof of Lemma 2 is similar to that of Lemma 1 (see "Appendix").

Remark 4 From Proposition 2 and Lemma 2, if  $e_{k+1} = \{i,j\} &gt; \{u,w\} = e_k$ , then, given  $x(s) \in X^{uw} \supset X^{ij}$ , to obtain  $x(s') \in X^{ij}$  it suffices to find  $s' \in B^{ij}$  such that  $\| x_i(s') - x_j(s') \| = d_{ij'}$ .

Furthermore, in the following we show that there is a unique  $s' \in B^{ij}$  satisfying such condition. For this, let us recall a simple fact that follows from Definition 1.

Proposition 3 If  $(G, K)$  is a  $^K$  DMDGP instance, so is  $(G[v_i, \ldots, v_j], K)$ , for  $j &gt; K + i$ .

Thus, given a  $^K$  DMDGP instance  $(G, K)$ , any subgraph induced by at least  $K + 1$  consecutive (w.r.t. the vertex order) vertices of  $V(G)$  is a  $^K$  DMDGP itself. Proposition 3 implies that each  $\{v_i, v_j\} \in E_P$  defines a DMDGP instance based on the subgraph  $G[v_i, \ldots, v_j]$ .

Proposition 4 Any DMDGP instance  $(G[v_i,\dots ,v_j],K)$  spanned by  $\{v_{i},v_{j}\} \in E_{P}$  has only two solutions.

Proof It follows from Proposition 3 and Corollary 1.  $\square$

Proposition 4 says that each DMDGP instance  $(G[v_i, \ldots, v_j], K)$  spanned by a pruning edge  $\{i, j\}$  has only two solutions, which are reflections of each other through the hyperplane defined by  $x_i, \ldots, x_{i + K - 1}$ . These two solutions correspond to a particular configuration of the components  $s_{i + K}^{\prime}, \ldots, s_{j}^{\prime}$ . The only difference between the two is the first component  $s_{i + K}^{\prime}$ . Since  $v_{i + K} \notin S^{ij}$  and the components of  $s_{\ell}^{\prime}$  with  $\ell \leq i + K$  or  $\ell &gt; j$  are kept fixed, we conclude that  $s^{\prime} \in B^{ij}$  is unique.

# 4 New Algorithm

First, we present a conceptual algorithm (Algorithm 2) which summarizes the ideas discussed in the previous sections. It requires a possible realization  $x \in \hat{X}$  and an arbitrary ordering  $e_1, \ldots, e_m$  for edges in  $E_P$ . Then, we describe a specialization of such algorithm (Algorithm 4) based on a specific ordering for the pruning edges which allows us to specify an efficient procedure for computing the sets  $S^{ij}$  along with some improvements in the computation of vertex positions.

# 4.1 The Conceptual Algorithm

Algorithm 2 describes our framework called Symmetry-based Build-up (SBBU). Given  $x \in \hat{X}$ , it handles the subproblems  $(G^{ij}, K)$  induced by pruning edges

Springer

---

Algorithmica (2021) 83:2400-2426

sequentially, according to a given order  $&lt;$  in  $E_P$ , until a solution to  $(G, K)$  is found.

Algorithm 2 SBBU
1: SBBU(G, K, (e1, ..., em), x ∈ X)
2: Set s = 0, x(0) = x
3: for k = 1, 2, ..., m do
4: {i, j} = e_k
5: if |S^ij| &gt; 0 then
6: Find s' ∈ B^ij : ||xi(s') - x_j(s')|| = dij
7: Update s = s' and x(s) = U(x, s)
8: end if
9: end for
10: return a valid realization x(s)

When solving subproblem  $(G^{ij},K)$ , if  $|S^{ij}| = 0$  then this subproblem has already been solved implicitly, as the following proposition shows.

Proposition 5 Let  $x(s)$  be a valid realization for  $(G^{uw}, K)$ , for all  $\{u, w\} \in P^{ij}$ . If  $S^{ij} = \emptyset$ , then  $x(s)$  is valid for  $(G^{ij}, K)$ .

Proof If  $S^{ij} = \emptyset$ , then for every  $\nu_{\ell} \in \{\nu_{i + K + 1}, \dots, \nu_j\}$ ,  $\exists \{u, w\} \in P^{ij}$  such that  $u + K + 1 \leq \ell \leq w$ . Suppose that  $x(s)$  is such that  $\| x_i(s) - x_j(s) \| \neq d_{ij}$ . (a) By Theorem 1 and Assumption 1,  $\exists s' \in \{0, 1\}^{n - K}$  with some  $s_{\ell}' \neq s_{\ell}$ , for  $u + K + 1 \leq \ell \leq w$ , such that  $\| x_i(s') - x_j(s') \| = d_{ij}$ . (b) Since  $\ell \geq u + K + 1$ , it follows that  $x_w(s') \notin \{x_w(s), R_s^{u + K}(x_w(s))\}$ . Thus, by Theorem 1,  $\| x_w(s') - x_u(s') \| \neq d_{uw}$ . But (a) and (b) together contradict Assumption 1. Hence  $\| x_i(s) - x_j(s) \| = d_{ij}$  and the assertion follows from Lemma 2.

Otherwise, for  $|S^{ij}| &gt; 0$ , in Step 6 we perform an exhaustive search to find  $s' \in B^{ij}$  such that  $\| x_i(s') - x_j(s') \| = d_{ij}$ . In Step 7, we update the current realization to  $x(s')$  according to Eq. (4).

Theorem 3 Let  $(G, K)$  be a feasible  $^K$  DMDGP instance. Considering exact arithmetic, Algorithm 2 finds  $x \in X$ .

Proof Since  $x(0) = x \in \hat{X}$ , due to Assumption 1 and Lemma 2, Step 6 is well-defined. From Remark 4 and Step 6, it follows that  $x(s') \in X^{ij}$ , for every  $e_k = \{i,j\}$ . Thus, since for the last pruning edge  $e_m = \{i_m,j_m\}$ , we have  $E^{i_m,j_m} = E$ , i.e.  $G^{i_m,j_m} = G$  after this last subproblem is solved,  $x(s) \in X^{i_m,j_m} = X(G) = X$ .

# 4.2 A Practical Algorithm

In this section, based on a particular pruning edge order, we introduce a practical version of Algorithm 2 which:

Springer

---

does not require an initial realization x ∈ X̂;avoids the computation and storage of unnecessary reflectors R_{x}^{i}(·);may result in fewer operations in the update step (Step 7) of Algorithm 2;allows us to discuss a concrete implementation for the sets S^{ij}.

For this, instead of working with a full realization x(s) ∈ X̂, which is updated through the binary vector s by Eq. (4), and computing and storing reflectors R_{x}^{i}(·) based on x = x(0) ∈ X̂, the idea is to grow a partial realization x_{1}, ... ,x_{t}, where t = arg max{w | {u, w} ∈ P^{ij}}, and compute the necessary reflectors on the fly based on the current partial realization and S^{ij}. This way, for each subproblem (G^{ij}, K), we do not compute valid full realizations x_{1}, ... ,x_{n} but valid partial realizations x_{1}, ... ,x_{j}, with j ≤ t.

### Assumption 2

Pruning edges {i,j}, with i < j, are sorted in increasing order of j, followed by a decreasing order of i.

Under this order, we can re-write the set of pruning edges preceding {i,j} as$$P^{ij} = \left{ {\left{ u,w} \right} \in E_{P} \mid u < w < j \vee (w = j \wedge w > u > i)} \right}.$$

### Definition 5

We say that x_{1}, ... ,x_{t}, with j ≤ t, is a valid partial realization for (G^{ij}, K), if x_{1}, ... ,x_{t} satisfies all distance constraints associated to edges in {{u, w} ∈ E(G^{ij}) | w ≤ j}.

### Remark 5

Recall that E(G^{ij}) = E^{ij} = E_{D} ∪ P^{ij} ∪ {i,j} and thanks to Assumption 2, there is no {u, w} ∈ P^{ij}, with w > j. This allows us to extend a valid partial realization x_{1}, ... ,x_{t} for (G^{ij}, K) to a valid full realization x_{1}, ... ,x_{n} for (G^{ij}, K), i.e x ∈ X^{ij}, by simply growing x_{1}, ... ,x_{t} to x_{1}, ... ,x_{n} using discretization distances (see Sect. 4.2.1), because no distance constraint {u, w} ∈ P^{ij} ∪ {i,j} is affected by this operation.

Assumption 2, along with (10), will be used in the results that follow. Using this concepts we will show in the next subsections that when dealing with subproblem (G^{ij}, K):given a partial realization x_{1}, ... ,x_{t} satisfying discretization distances and distances corresponding to pruning edges in P^{ij}, it can be extended to x_{1}, ... ,x_{t}, ... ,x_{j} keeping feasibility of such distance constraints and new discretization constraints;it is possible to apply partial reflections to this extended partial realization in order to fulfill ||x_{i} - x_{j}|| = d_{ij} without violating the distance constraints considered so far.

---

Algorithmica (2021) 83:2400-2426

## 4.2.1 Initialization of Candidate Positions

In Sect. 4.2.2 we shall explain how to find a valid partial realization for DMDGP subproblems $(G^{ij}, K)$ by composing reflections through symmetry hyperplanes and applying them to positions $x_{i+K+1}, \ldots, x_j$. This procedure assumes that candidate positions for $x_i, \ldots, x_j$ are available when we start to solve $(G^{ij}, K)$. In this section, we describe how to initialize these positions.

From now on, we assume that initialization of candidate positions must follow the vertex order from $v_{1}$ to $v_{n}$, meaning that if $(G^{ij}, K)$ is the current subproblem, and $x_{t}$ is the last initialized position, such that $t &lt; j$, then we initialize positions from $x_{t+1}$ to $x_{j}$, whereas $x_{1}, \ldots, x_{t}$ remain unchanged. In other words, the candidate positions $x_{t+1}, \ldots, x_{j}$ are grown from the current partial realization $x_{1}, \ldots, x_{t}$ using only distance constraints associated to discretization edges. Moreover, each position $x_{i}$ is initialized only once, although it can be modified later (see Sect. 4.2.2) in order to satisfy a distance constraint corresponding to a pruning edge $e' \geq e = \{i, j\}$. This is formalized in Proposition 6.

**Proposition 6** Assume edges in $E_P$ are ordered as $e_1, \ldots, e_m$. Then, before solving $(G^{ij}, K)$, positions $x_1, \ldots, x_j$ can be initialized such that

$$
\forall \{\ell , k \} \in E _ {D} \cap E (G [ v _ {1}, \dots , v _ {j}]), \quad \| x _ {\ell} - x _ {k} \| = d _ {\ell k}, \tag {11}
$$

and

$$
\forall \{\ell , k \} \in P ^ {i j}, \quad \| x _ {\ell} - x _ {k} \| = d _ {\ell k}. \tag {12}
$$

**Proof** We prove this by induction on the edge order. In the base case we consider $e_1 = \{i_1, j_1\} \in E_P$ spanning the first subproblem to be solved. The positions $x_i$, for $i = 1, \ldots, j_1$ are initialized right away. From Definition 1, $x_1, \ldots, x_K$ can be localized uniquely (up to rotations and translations) by different methods [1, 10]. Hence, $\| x_\ell - x_k \| = d_{\ell k}, \forall \{\ell, k\} \in E(G[v_1, \ldots, v_K])$. Then, by $K$-lateration (see Remark. 1), there are at most two positions $\{x_i^+, x_i^-\}$ for $v_i \in V$ for each $K &lt; i \leq j_1$. Notice that any partial realization $x_1, \ldots, x_{j_1}$ is enough to build partial reflections. The correct alternative will be chosen later by the appropriate partial reflection composition which satisfies constraints defined by pruning edges (Sect. 4.2.2 gives more details). Thus, without loss of generality, let $x_1, \ldots, x_{j_1}$ be the partial realization obtained by choosing $x_i^-$, for every $i = K + 1, \ldots, j_1$. Since $P^{i_1 j_1} = \emptyset$, this partial realization satisfies (11) and (12) for $\{i, j\} = \{i_1, j_1\}$.

The induction hypothesis is that (11) and (12) hold for pruning edges $e_1, \ldots, e_k$, i.e. $x_1, \ldots, x_t$ is a valid partial realization for all subproblems spanned by these edges, where

$$
t := \max  \{w \mid \{u, w \} \in P ^ {i j} \} = \max  \{w \mid \{u, w \} \in \{e _ {1}, \dots , e _ {k} \} \}.
$$

In the inductive step, let us prove that (11) and (12) also hold for pruning edge $e_{k+1} = \{i,j\}$ spanning subproblem $(G^{ij}, K)$.

Springer

---

Algorithmica (2021) 83:2400-2426

Since subproblems spanned by edges in  $P^{ij}$  are solved, positions  $x_{1},\ldots ,x_{t}$  are already initialized and satisfy (12), and (11) with  $\nu_{j} = \nu_{t}$ .

If  $j \leq t$ , then there is nothing left to do. Thus, suppose  $j &gt; t$ . Then, positions  $x_{t+1}, \ldots, x_j$  can be initialized by  $K$ -lateration (see Remark 1 and Algorithm 3) based on discretization distances such that (11) holds.

Remark 6 The proof of Proposition 6 describes a procedure for initialization of  $x_{1}, \ldots, x_{j}$  before solving  $(G^{ij}, K)$ . It is important to notice that such initialization is done sequentially and depends on previously computed positions which are not recomputed in this step. Thus, after the initialization, the current partial realization continues to be valid for all already solved subproblems.

Algorithm 3 InitializePositions
1: InitializePositions(x,t,j)
2: if  $t \geq j$  then
3: return  $x$  and  $t$ .
4: end if
5: if  $t = 0$  then
6: Initialize  $x_{1}, \ldots, x_{K}$  as a solution of  $\| x_{\ell} - x_{i} \|^{2} = d_{\ell i}^{2}, \forall \{\ell, i\} \in E(G[v_{1}, \ldots, v_{K}])$
7: Set  $t = K$
8: end if
9: for  $i = t + 1, \ldots, j$  do
10: Find solutions  $\{x_{i}^{+}, x_{i}^{-}\}$  for the system:  $\| x_{\ell} - x_{i} \|^{2} = d_{\ell, i}^{2}, \ell = i - K, \ldots, i - 1$ .
11: Set  $x_{i} = x_{i}^{-}$
12: end for
13: Set  $t = j$
14: return  $x$  and  $t$ .

Algorithm 3 gives a pseudocode for the function InitializePositions which receives a current realization  $x$  (actually, the current partial realization  $x_{1},\ldots ,x_{t}$ ), the index of the last initialized position  $t$ , and the index  $j$  of the last vertex whose position needs initialization. Updated  $x$  and  $t$  are returned by this function.

# 4.2.2 Solving DMDGP Subproblems

Now we explain how to find a valid partial realization for a DMDGP subproblem  $(G^{ij},K)$ , given a valid partial realization  $x_{1},\ldots ,x_{t}$  for  $(G^{uw},K)$ ,  $\forall \{u,w\} \in P^{ij}$ .

Recall from Proposition 1 that given positions  $x_{i+1}, \ldots, x_{i+K+1}, \ldots, x_j$ , valid or not, we can build all necessary symmetry hyperplanes and their corresponding reflection operators  $R_x^{i+K+1}(\cdot), \ldots, R_x^j(\cdot)$ . Then, based on Theorem 1 and Proposition 4, we can apply compositions of such reflection operators only to  $x_j$  until we find its correct position, as illustrated in Fig. 2 for the 2D case.

The only decision to be taken is whether each of the reflectors  $R_{x}^{i + K + 1}(\cdot), \dots, R_{x}^{j}(\cdot)$  should be applied or not to  $x_{j}$  in order to fulfill  $\| x_{i} - R(x_{j},\bar{x})\| = d_{ij}$ , where  $R(.,\bar{x})$  is a composition of the chosen reflectors:

Springer

---

Algorithmica (2021) 83:2400-2426

$$
R(y, \bar{s}) := \left(R_{x}^{i + K + 1}\right)^{\bar{s}_{1}} \left(R_{x}^{i + K + 2}\right)^{\bar{s}_{2}} \dots \left(R_{x}^{i}\right)^{\bar{s}_{j - i - K}}(y). \tag{13}
$$

In (13), the binary vector $\bar{s}$ is of size $j - i - K$, and $(R_x^\ell)^0 := I$, the identity operator in $\mathbb{R}^K$, whereas $(R_x^\ell)^1 := R_x^\ell$, for $\ell = i + K + 1, \ldots, j$.

In contrast to Algorithm 2, where $s$ is the global binary decision variable and all the reflectors are computed based on the first realization $x(0) = x \in \hat{X}$, now the reflectors $R_{x}^{i + K + \ell}(\cdot)$ for which $v_{i + K + \ell} \in S^{ij}$ are computed based on the current partial realization $x_{1}, \ldots, x_{t}$, and $\bar{s}$ is a local binary decision variable belonging to

$$
\tilde{B}^{ij} := \left\{\bar{s} \in \{0, 1\}^{j - i - K} \mid \bar{s}_{\ell} = 0 \text{ if } v_{i + K + \ell} \notin S^{ij} \right\}. \tag{14}
$$

Thus, we look for a binary vector $\bar{s} \in \tilde{B}^{ij}$ such that

$$
x_{j}' = R(x_{j}, \bar{s}) = \left(R_{x}^{i + K + 1}\right)^{\bar{s}_{1}} \left(R_{x}^{i + K + 2}\right)^{\bar{s}_{2}} \dots \left(R_{x}^{i}\right)^{\bar{s}_{j - i - K}}(x_{j}) \tag{15}
$$

satisfies $\| x_{i} - x_{j}'\| = d_{ij}$. We remark that this search is exhaustive: we test all $|\tilde{B}^{ij}| = 2^{|S^0|}$ possible choices for $\bar{s}$ (recall that there is a unique $\bar{s}$ that works, as discussed after Proposition 4).

Algorithm 4 SBBU
1: SBBU(G,K)
2: Order edges $\{v_{i},v_{j}\} \in E_{P}$ in increasing order of $j$ and decreasing order of $i$ obtaining a sequence $(e_1,\ldots ,e_m)$, with $m = |E_P|$. Set $t = 0$, $n = |V|$
3: Set $\mathcal{C} = \{\{v_i\} \}_{i = K + 1}^n$
4: InitializePositions(x,t,K) # positions for the initial clique
5: for $k = 1,2,\dots,m$ do
6: $\{i,j\} = e_k$
7: if $\rho_{\mathcal{C}}(i + K) \neq \rho_{\mathcal{C}}(j)$ then
8: InitializePositions(x,t,j)
9: Set $C^0 = \rho_{\mathcal{C}}(i + K)$ and $\mathcal{D} = \rho_{\mathcal{C}}(\{i + K + 1,\dots,j\}) \setminus \{C^0\}$
10: Let $S^{ij} = \cup_{C\in \mathcal{D}}\mathrm{first}(C)$ # local symmetry vertices
11: Compute $R_x^\ell (\cdot)$ for each $v_{\ell}\in S^{ij}$
12: Find $s\in \tilde{B}^{ij}: \| x_i - R(x_j,\bar{s})\| = d_{ij}$ # find position $x_j$
13: Update $x_{\ell} = \tilde{U} (x_{\ell},\bar{s})$, for $\ell = i + K + 1,\ldots ,j$
14: Set $C^+ = (\cup_{C\in \mathcal{D}}C)\cup C^0$ and update $\mathcal{C} = (\mathcal{C}\setminus (\mathcal{D}\cup \{C^0\}))\cup C^+$
15: end if
16: end for
17: if $t &lt; n$ then
18: InitializePositions(x,t,n)
19: end if
20: return a valid realization $x$

Once $\bar{s}$ is found, the positions of $v_{i + K + 1},\ldots ,v_{t}$ are updated according to:

$$
x_{\ell}' = \tilde{U}(x_{\ell}, \bar{s}) := \left(\prod_{t = 1}^{\ell - i - K} \left(R_{x}^{i + K + t}\right)^{\bar{s}_{t}}\right)(x_{\ell}), \quad \ell = i + K + 1, \dots , t. \tag{16}
$$

This update maintain feasibility of $x_{1}, \ldots, x_{t}$ with respect to $(G^{uv}, K)$, for every $\{u, w\} \in P^{ij}$, because positions $x_{u + K}, \ldots, x_{w}$ are only updated simultaneously by partial reflections $R_{x}^{v}(x_{u + K}), \ldots, R_{x}^{v}(x_{w})$, for $v \leq u + K$.

Springer

---

Algorithmica (2021) 83:2400-2426

# 4.2.3 Symmetry Vertex Sets

The ideas of the Sects. 4.2.1 and 4.2.2 lead to Algorithm 4. This algorithm makes use of  $\mathcal{C}$ , a partition of  $\{v_{K + 1},\ldots ,v_n\}$  used to obtain the sets  $S^{ij}$ . At the beginning, we set  $\mathcal{C} = \{\{v_i\} \}_{i = K + 1}^n$ . This partition is updated in Step 14 taking into account already solved subproblems. Assume that subsets of vertices in  $\mathcal{C}$  are ordered according to the vertex order of Definition 1. Let us denote the first vertex of  $C^0\in \mathcal{C}$  by first  $(C^0)$ .

We also introduce a function  $\rho_{\mathcal{C}}: \{K + 1, \ldots, n\} \to \mathcal{C}$ , parametrized by  $\mathcal{C}$ , such that  $\rho_{\mathcal{C}}(\ell)$  returns the unique element of  $\mathcal{C}$  containing vertex  $\nu_{\ell}$ . The next proposition shows that this function is well-defined at every iteration of Algorithm 4.

Proposition 7 At every iteration of Algorithm 4,  $\mathcal{C}$  is a partition of the subset of vertices  $\{v_{K + 1},\ldots ,v_n\}$ .

Proof At the first iteration  $\mathcal{C} = \{\{v_i\}\}_{i = K + 1}^n$ . Assume that, at the beginning of iteration  $k$ ,  $\mathcal{C}$  is a partition of  $\{v_{K + 1},\ldots ,v_n\}$ . If  $\rho_{\mathcal{C}}(i + K) = \rho_{\mathcal{C}}(j)$ , then we go to the next iteration with  $\mathcal{C}$  unchanged. Otherwise, from Step 9,  $\mathcal{D}$  and  $C^0$  are subsets of  $\mathcal{C}$ . Then, Step 14 updates  $\mathcal{C}$  by removing these subsets and including their union, hence, the updated  $\mathcal{C}$  is still a partition of  $\{v_{K + 1},\ldots ,v_n\}$ .

Proposition 8 In Algorithm 4, after Step 14, there exists a unique  $C \in \mathcal{C}$  such that  $C \supset \{v_{i + K}, \ldots, v_j\}$ .

Proof Follows directly from Steps 9 and 14 of Algorithm 4.

We remark that  $\rho_{\mathcal{C}}(\{i + K + 1,\dots ,j\})$  denotes the image of  $\{i + K + 1,\ldots ,j\}$  by  $\rho_{\mathcal{C}}$  in the definition of  $\mathcal{D}$  (Step 9), i.e it returns elements of  $\mathcal{C}$  whose union contains  $v_{i + K + 1},\ldots ,v_{j}$  and  $C^0 = \rho_{\mathcal{C}}(i + K)$  is the element of  $\mathcal{C}$  containing  $v_{i + K}$ .

Proposition 9 In Algorithm 4,  $\rho_{\mathcal{C}}(i + K) = \rho_{\mathcal{C}}(j)$  if and only if  $S^{ij} = \emptyset$

Proof If  $S^{ij} = \emptyset$ , then

$$
\forall v _ {\ell} \in \left\{v _ {i + K + 1}, \dots , v _ {j} \right\}, \quad \exists \{u, w \} \in P ^ {i j}: u + K &lt;   \ell \leq w \leq j. \tag {17}
$$

In particular, for  $\ell = i + K + 1$ , there exists  $\{r,z\} \in P^{ij}$  such that  $r + K &lt; i + K + 1 \leq z \leq j$ . Clearly,  $r \leq i$ . From Proposition 8 there exists a unique  $C^1 \in \mathcal{C}$  such that  $C^1 \supset \{v_{r + K}, \ldots, v_{i + K}, v_{i + K + 1}, \ldots, v_z\}$ .

Thus, if  $z = j$ , then  $\rho_{\mathcal{C}}(i + K) = \rho_{\mathcal{C}}(j)$ .

Otherwise, for  $z &lt; j$ , because  $i + K + 1 \leq z$ , it follows that  $v_{z+1} \in \{v_{i+K+1}, \ldots, v_j\}$  and from (17), there exists  $\{u, w\} \in P^{ij}$  such that  $u + K &lt; z + 1 \leq w \leq j$  (clearly,  $u + K \leq z$ ). Thus, from Proposition 8:

$$
\exists ! C ^ {2} \supset \left\{v _ {u + K}, \dots , v _ {z}, v _ {z + 1}, \dots , v _ {w} \right\}. \tag {18}
$$

If  $u\leq r\leq i$  , then from (18), we obtain  $v_{i + K}\in C^2$

Springer

---

Algorithmica (2021) 83:2400-2426

Otherwise, for  $r &lt; u \leq z - K$ , then  $r + K &lt; u + K \leq z$ , implying that  $\nu_{u + K} \in C^1$ . In either case, we have  $\rho_{\mathcal{C}}(u + K) = \rho_{\mathcal{C}}(i + K)$ . From Proposition 7, we conclude that  $\rho_{\mathcal{C}}(i + K) = \rho_{\mathcal{C}}(w)$ .

Hence, if  $w = j$ ,  $\rho_{\mathcal{C}}(i + K) = \rho_{\mathcal{C}}(j)$ .

Otherwise  $(w &lt; j)$ , in view of (17), we can apply the same argument to  $\nu_{w+1}$ , and repeat until we find  $\{h, p\} \in P^{ij}$  with  $p = j$ .

On the other hand, to prove that  $\rho_{\mathcal{C}}(i + K) = \rho_{\mathcal{C}}(j)$  implies  $S^{ij} = \emptyset$ , we use the counter-positive. Suppose there exists  $\nu_{\ell} \in \{\nu_{i + K + 1}, \dots, \nu_j\}$  such that  $\nexists \{u, w\} \in P^{ij}$  such that  $u + K &lt; \ell \leq w \leq j$ . This means that  $\forall \{u, w\} \in P^{ij}$  either (i)  $w &lt; \ell$  or (ii)  $\ell \leq w \leq j$  and  $u + K \geq \ell$ . If  $w &lt; \ell, \forall \{u, w\} \in P^{ij}$ , then  $\rho_{\mathcal{C}}(j) = \{\nu_j\} \neq \rho_{\mathcal{C}}(i + K)$ , because  $w &lt; \ell \leq j$  ( $\nu_{\ell}$  and  $\nu_j$  were never reached).

Thus, let us consider  $\{u,w\} \in P^{ij}$  such that  $\ell \leq w\leq j$  and  $u + K\geq \ell$ . Without loss of generality, assume  $w = j$ . Since  $\ell \geq i + K + 1$ , then  $u + K\geq i + K + 1$  (or  $u\geq i + 1$ ), implying that  $\rho_{\mathcal{C}}(i + K)\neq \rho_{\mathcal{C}}(u + K) = \rho_{\mathcal{C}}(j)$ , where the last equality follows from Proposition 8.

Proposition 9 shows that if  $\nu_{i + K}$  and  $\nu_{j}$  are in the same subset, i.e.,  $\rho_{\mathcal{C}}(i + K) = \rho_{\mathcal{C}}(j)$ , then this subproblem was already solved implicitly (see Proposition 5). This is equivalent to condition  $|S^{ij}| = 0$  in Algorithm 2.

Otherwise, we need to obtain the set  $S^{ij}$  of symmetry vertices for  $(G^{ij}, K)$ . This is accomplished in Step 10.

Theorem 4 If  $\rho_{\mathcal{C}}(i + K) \neq \rho_{\mathcal{C}}(j)$ , then  $\bigcup_{C \in \mathcal{D}} \mathrm{first}(C) = S^{ij}$ .

Proof Let  $\nu_{\ell} \in \cup_{C \in \mathcal{D}} \mathrm{first}(C)$ . From Proposition 7,  $\exists! \hat{C} \supset \{\nu_{\ell}\}$  such that  $\hat{C} \in \mathcal{D}$  and  $\mathrm{first}(\hat{C}) = \nu_{\ell}$ . Suppose  $\nu_{\ell} \notin S^{ij}$ , i.e. there exists  $\{u, w\} \in P^{ij}$  such that  $u + K &lt; \ell \leq w$ . From Proposition 8,  $\exists! C \in \mathcal{C}$  such that  $C \supset \{\nu_{u + K}, \dots, \nu_w\}$  and since  $\mathcal{C}$  is a partition, it follows that  $C = \hat{C}$ . But  $\mathrm{first}(\hat{C}) = \mathrm{first}(C) = \nu_{u + K} \neq \nu_{\ell}$  contradicting  $\mathrm{first}(\hat{C}) = \nu_{\ell}$ . Therefore,  $\nexists \{u, w\} \in P^{ij}$  such that  $u + K &lt; \ell \leq w$ . Thus,  $\nu_{\ell} \in S^{ij}$ .

Conversely, let  $\nu_{\ell} \in S^{ij}$ . Then, for every  $\{u, w\} \in P^{ij}$  either (i)  $w &lt; \ell$  or (ii)  $\ell \leq w \leq j$  and  $u + K \geq \ell$ . If  $\forall \{u, w\} \in P^{ij}$ , we have  $w &lt; \ell$ , then  $C = \rho_{\mathcal{C}}(\ell) = \{\nu_{\ell}\} \in \mathcal{D}$  and  $\mathrm{first}(C) = \nu_{\ell}$ . Otherwise, there are  $\{u, w\} \in P^{ij}$  such that  $\ell \leq w \leq j$ . For all of these,  $u + K \geq \ell$ . Recall from Algorithm 4 that  $\mathrm{first}(C^{+}) = \mathrm{first}(C^{0}) = \mathrm{first}(\rho_{\mathcal{C}}(u + K))$ . We split the analysis in three cases.

Case 1:  $\nu_{\ell} \notin \rho_{\mathcal{C}}(u + K)$ . Then,  $\ell &lt; u + K$ , implying that  $\nu_{\ell} &lt; \mathrm{first}(u + K)$ . Thus, after iteration  $k$  with  $e_k = \{u, w\}$ , we have  $\rho_{\mathcal{C}}(\ell) = \{\nu_{\ell}\}$ .

Case 2:  $\nu_{\ell} \in \rho_{\mathcal{C}}(u + K)$  and  $\mathrm{first}(\rho_{\mathcal{C}}(u + K)) = \nu_{\ell}$ . In this case, after iteration  $k$  with  $e_k = \{u, w\}$ ,  $\rho_{\mathcal{C}}(\ell) = \rho_{\mathcal{C}}(u + K) = \{\nu_{\ell}, \dots, \nu_{u + K}, \dots, \nu_p\}$ . Thus  $\mathrm{first}(\rho_{\mathcal{C}}(\ell)) = \nu_{\ell}$ .

Case 3:  $\nu_{\ell} \in \rho_{\mathcal{C}}(u + K)$  but  $\mathrm{first}(\rho_{\mathcal{C}}(u + K)) &lt; \nu_{\ell}$ . In this case,  $C^0 = \rho_{\mathcal{C}}(\ell) = \rho_{\mathcal{C}}(u + K) = \{\nu_q, \dots, \nu_\ell, \dots, \nu_{u + K}, \dots, \nu_p\}$ . The set  $C^0$  is the result of iteration  $k'$ , with  $e_{k'} = \{h, p\} &lt; \{u, w\} = e_k$ . Clearly  $\{h, p\} \in P^{ij}$ . From Proposition 8,  $\exists! C \supset \{\nu_{h + K}, \dots, \nu_p\}$ . Notice that  $u + K \leq p \leq w$ . Since  $q &gt; h + K$  contradicts the fact that  $\mathrm{first}(\rho_{\mathcal{C}}(u + K)) = \nu_q$ , then  $q \leq h + K &lt; \ell$ . This leads to  $h + K &lt; \ell \leq u + K \leq p$  which implies that  $\nu_{\ell} \notin S^{ij}$ , a contradiction.

Springer

---

Algorithmica (2021) 83:2400-2426

Therefore, only cases 1 and 2 can happen and both imply in  $\nu_{\ell} \in \bigcup_{c \in \mathcal{D}} \mathrm{first}(C)$ .

Corollary 2 Let  $(G, K)$  be a feasible  $^K$  DMDGP instance. Considering exact arithmetic, Algorithm 4 finds a valid realization  $x$  for  $(G, K)$ .

Proof From Theorem 4 and Remark 5 correctness of Algorithm 4 follows from Theorem 3.

In the end, we obtain a valid partial realization  $x_{1}, \ldots, x_{t}$  for all subproblems  $(G^{ij}, K)$ , with  $\{i, j\} \in E_{P}$ . If  $t = n$ , we are done. Otherwise, in view of Remark 5,  $x_{1}, \ldots, x_{t}$  can be extended to a valid realization  $x \in X$ . This explains Step 18.

Even under Assumption 1, due to floating point arithmetic, in Step 12 we may not be able to find  $s$  such that  $\| x_{i} - R(x_{j},\bar{s})\| = d_{ij}$ . Thus, instead of stopping as soon as we find a  $s$  such that  $\| x_{i} - R(x_{j},\bar{s})\| -d_{ij}\| \leq \varepsilon$ , for a prescribed tolerance  $\varepsilon$ , we actually consider all  $2^{|S^ij|}$  possibilities and choose  $s$  for which  $\| x_{i} - R(x_{j},\bar{s})\| -d_{ij}\|$  is minimum. In case  $\| x_{i} - R(x_{j},\bar{s})\| -d_{ij}\| &gt;\varepsilon$  for every  $\bar{s}\in \tilde{B}^{ij}$ , then we actually interrupt the algorithm and return "failure". However, this never happened in the numerical experiments of Sect. 5.

We remark that  $|S^{ij}|$  is a good indicator of the computational cost for solving subproblem  $(G^{ij}, K)$ , because it determines the number of  $2^{|S^ij|}$  reflection compositions that we need to apply to  $x_j$  in order to find its correct position. Thus, we define the corresponding total work to solve a  $^K$ DMDGP instance as

$$
W := \sum_ {\{i, j \} \in \hat {E}} 2 ^ {| S ^ {i j} |}, \tag {19}
$$

where  $\hat{E} = \{\{v_i,v_j\} \in E_P\mid |S^{ij}| &gt; 0\}$ . Let us also denote by  $\hat{W} = \max_{\hat{E}}2^{|S^ij|}$ , the maximum work per pruning edge. We also remark that  $|S^{ij}|$  depends on the order in which the pruning edges are handled and, in this paper, we consider only the order described in Assumption 2.

Nevertheless, it is easy to illustrate worst case and best case scenarios. For example, if  $\{1,n\} \in E_P$  and it is the first pruning edge to be handled, then  $W = 2^{|S^{1n}|} = 2^{n - K}$  incurring in an exponential cost in  $n$ . On the other hand, if  $\{i - K - 1,i\} \in E_P$ , for  $i = K + 2,\ldots ,n$ , and they are handled in the order  $\{1,K + 2\},\ldots ,\{n - K - 1,n\}$ , then  $W = \sum_{i = K + 2}^{n}2^{1} = 2(n - K - 1)$ , which is linear in  $n$ .

# 5 Experimental Results

An efficient implementation of the function  $\rho_{\mathcal{C}}$  needs to deal with its evaluation and the update of the subsets of  $\mathcal{C}$ . We adopted the structure proposed by Newman and Ziff [36], which allows the evaluation of  $\rho_{\mathcal{C}}$  and the subsets update in time  $O(\log_2(|V|))$  and memory  $O(|V|)$ .

Springer

---

In order to validate Algorithm 4 and assess its performance, we generate a set of protein-like instances (K = 3) whose data were extracted from Protein Data Bank (PDB) [3], and compare the results with those of the classic BP [16, 25]. For each protein, we consider only the backbone composed by the sequence of atoms N - C_{α} - C and include an edge in the corresponding graph:either when the atoms are separated by at most three covalent bondsor the distance between pairs of atoms is smaller than a certain cut-off value.

The resolution of NMR experiments usually varies between 5 Å and 6 Å. The smaller the cut-off value, the sparser the DMDGP instance. Each instance was generated by the first model and first chain of the PDB file.

The natural backbone order for instances generated in this way provides a vertex order satisfying the assumptions of Definition 1, implying we are working with DMDGP instances.

In our experiments we consider two test sets: one using cut-off 6 Å and other using 5 Å . In Tables 1 and 2, we present the results obtained by both algorithms: BP is the classic Branch-and-Prune implementing a depth-first search [16, 25], whereas SBBU (Symmetry-based Build-up) corresponds to Algorithm 4.

The algorithms were implemented in C++ and the experiments carried out in a machine Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz, 8G RAM, running Linux Ubuntu 20.04 LTS, GCC compiler version 9.3.0. Our codes and datasets can be found in https://github.com/michaelsouza/sbbu, whereas an implementation of the Branch-and-Prune algorithm is available at https://github.com/mucherino/mdjeep.

Both tables bring the ID of the protein in PDB, the number of atoms |V|, number of edges (available distances) |E|, CPU time in seconds for the two algorithms and the normalized Mean Distance Error (MDE):$$MDE\left(X,E,d\right)=\frac{1}{\left|E\right|}\sum_{\left\{i,j\right]\in E}\frac{\left| {\left\| {x}_{i}-{x}_{j}\right\|}_{2}-{d}_{ij}\right|}{{d}_{ij}}.$$

Both algorithms were stopped as soon as the first solution is found and a “--” symbol means that the algorithm was not able to find a solution in less than 300 seconds. For each instance, we also present the total and maximum works W and $\tilde{W}$, respectively. The last column, called “Speed-up”, contains the ratio time BP / time SBBU.

From the tables, we observe that the new algorithm provides a non-trivial speed-up in most of the instances. In particular, the new algorithm is considerably faster than the classic BP for the sparser instances where it was up to 100 times faster.

Concerning the estimated total work of SBBU, it seems that the time varies linearly with W as depicted in Fig. 3. The relationship between BP time and W and/or $\tilde{W}$ is not so clear. However, we argue that while the most costly subproblem (G^{ij}, K) represents

---

Algorithmica (2021) 83:2400-2426

Table 1 Computational results in some protein-like instances (cut-off:  $6\mathring{\mathrm{A}}$  )

|  ID | |V| | |E| | % | |S| | BP |   | SBBU |   |   |   | Speed-up  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|   |   |   |   |   |  Time | MDE | Time | MDE | W | W  |   |
|  1N6T | 30 | 236 | 0.543 | 1 | 1.48E-04 | 1.74E-11 | 1.89E-05 | 8.09E-12 | 2 | 52 | 7.83  |
|  1FW5 | 60 | 558 | 0.325 | 1 | 1.83E-04 | 4.68E-12 | 2.41E-05 | 3.47E-12 | 2 | 112 | 7.60  |
|  1ADX | 120 | 1008 | 0.141 | 1 | 4.34E-04 | 4.17E-12 | 7.59E-05 | 6.11E-12 | 2 | 232 | 5.72  |
|  1BDO | 241 | 2167 | 0.075 | 1 | 5.88E-04 | 8.13E-12 | 9.24E-05 | 1.32E-11 | 2 | 474 | 6.37  |
|  1ALL | 480 | 4932 | 0.043 | 1 | 1.04E-03 | 9.45E-13 | 1.90E-04 | 3.44E-12 | 2 | 952 | 5.49  |
|  6S61 | 522 | 5298 | 0.039 | 1 | 9.91E-04 | 7.19E-13 | 2.37E-04 | 7.10E-12 | 2 | 1036 | 4.18  |
|  1FHL | 1002 | 9811 | 0.019 | 1 | 1.80E-03 | 6.94E-12 | 4.00E-04 | 1.71E-11 | 2 | 1996 | 4.50  |
|  4WUA | 1033 | 9727 | 0.018 | 1 | 2.07E-03 | 7.74E-12 | 5.90E-04 | 6.10E-12 | 8 | 2060 | 3.51  |
|  6CZF | 1494 | 14163 | 0.013 | 1 | 2.65E-03 | 1.49E-12 | 6.08E-04 | 4.14E-12 | 2 | 2980 | 4.36  |
|  5IJN | 1950 | 18266 | 0.009 | 1 | 3.63E-03 | 1.07E-12 | 8.66E-04 | 1.79E-11 | 16 | 3908 | 4.19  |
|  6RN2 | 2052 | 19919 | 0.009 | 1 | 4.42E-03 | 8.98E-13 | 9.28E-04 | 1.16E-11 | 16 | 4104 | 4.76  |
|  1CZA | 2694 | 26452 | 0.007 | 1 | 4.56E-03 | 3.63E-12 | 1.10E-03 | 6.11E-11 | 2 | 5380 | 4.15  |
|  6BCO | 2856 | 27090 | 0.007 | 5 | 5.03E-03 | 1.11E-12 | 1.14E-03 | 1.45E-11 | 16 | 5730 | 4.42  |
|  1EPW | 3861 | 35028 | 0.005 | 1 | 6.70E-03 | 1.81E-11 | 1.49E-03 | 1.94E-10 | 2 | 7714 | 4.48  |
|  5NP0 | 7584 | 80337 | 0.003 | 1 | 3.48E-02 | 4.80E-12 | 3.50E-03 | 1.47E-10 | 256 | 15562 | 9.95  |
|  5NUG | 8760 | 82717 | 0.002 | 1 | 2.27E-02 | 3.99E-12 | 4.64E-03 | 6.07E-10 | 16 | 17592 | 4.90  |
|  4RH7 | 9015 | 85831 | 0.002 | 1 | 2.04E-02 | 1.80E-12 | 3.64E-03 | 2.15E-10 | 16 | 18054 | 5.60  |
|  3VKH | 9126 | 87621 | 0.002 | 1 | 4.99E-01 | 3.52E-11 | 3.72E-03 | 8.92E-10 | 256 | 18556 | 134.42  |

Springer

---

Algorithmica (2021) 83:2400-2426

Table 2 Computational results in some protein-like instances (cut-off:  $5\AA$

|  ID | |V| | |E| | % | |S| | BP |   | SBBU |   |   |   | Speed-up  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|   |   |   |   |   |  Time | MDE | Time | MDE | W | W  |   |
|  1N6T | 30 | 176 | 0.405 | 1 | 1.36E-04 | 9.08E-12 | 1.01E-05 | 5.64E-12 | 2 | 52 | 13.51  |
|  1FW5 | 60 | 417 | 0.236 | 1 | 2.33E-04 | 2.44E-12 | 2.21E-05 | 2.48E-12 | 2 | 112 | 10.51  |
|  1ADX | 120 | 659 | 0.092 | 1 | 2.99E-04 | 1.87E-12 | 3.76E-05 | 2.93E-12 | 2 | 232 | 7.96  |
|  1BDO | 241 | 1345 | 0.047 | 1 | 5.68E-04 | 3.39E-12 | 9.74E-05 | 8.49E-12 | 2 | 474 | 5.83  |
|  1ALL | 480 | 3443 | 0.029 | 1 | 1.87E-03 | 2.51E-13 | 2.66E-04 | 1.23E-12 | 2 | 952 | 7.03  |
|  6S61 | 522 | 3699 | 0.027 | 1 | 9.13E-04 | 2.20E-13 | 2.35E-04 | 2.96E-12 | 2 | 1036 | 3.89  |
|  1FHL | 1002 | 6378 | 0.013 | 1 | 2.16E-03 | 2.47E-12 | 3.82E-04 | 1.02E-11 | 2 | 1996 | 5.65  |
|  4WUA | 1033 | 6506 | 0.012 | 1 | 2.82E-03 | 2.73E-12 | 5.09E-04 | 2.91E-12 | 16 | 2066 | 5.54  |
|  6CZF | 1494 | 9223 | 0.008 | 1 | 2.87E-03 | 5.55E-13 | 4.56E-04 | 2.35E-12 | 2 | 2980 | 6.28  |
|  5IJN | 1950 | 11981 | 0.006 | 1 | 8.67E-03 | 3.50E-13 | 6.13E-04 | 7.24E-12 | 16 | 3908 | 14.14  |
|  6RN2 | 2052 | 13710 | 0.007 | 1 | 6.91E-03 | 3.06E-13 | 7.68E-04 | 6.22E-12 | 16 | 4112 | 9.01  |
|  1CZA | 2694 | 17451 | 0.005 | 1 | 5.22E-03 | 1.16E-12 | 8.19E-04 | 3.22E-11 | 2 | 5380 | 6.38  |
|  6BCO | 2856 | 18604 | 0.005 | 12 | 5.29E-03 | 4.62E-13 | 9.21E-04 | 8.95E-12 | 16 | 5706 | 5.75  |
|  1EPW | 3861 | 23191 | 0.003 | 1 | 6.90E-03 | 6.67E-12 | 1.86E-03 | 9.52E-11 | 8 | 7716 | 3.72  |
|  5NP0 | 7584 | 59478 | 0.002 | 1 | 2.57E-01 | 1.84E-12 | 2.88E-03 | 4.37E-11 | 256 | 16138 | 89.08  |
|  5NUG | 8760 | 56979 | 0.001 | 1 | 5.58E-02 | 1.27E-12 | 2.74E-03 | 1.78E-10 | 128 | 17700 | 20.39  |
|  4RH7 | 9015 | 59346 | 0.001 | 1 | 2.02E-02 | 6.30E-13 | 2.94E-03 | 8.69E-11 | 32 | 18068 | 6.89  |
|  3VKH | 9126 | 59592 | 0.001 | 1 | - | - | 2.38E-02 | 1.72E-09 | 65536 | 84066 | -  |

Springer

---

Algorithmica (2021) 83:2400-2426

![img-2.jpeg](img-2.jpeg)
Fig. 3 Scatter plots  $W \times$  time and linear regressions for the two datasets (Table 1 on the left, Table 2 on the right but not considering the last row)

![img-3.jpeg](img-3.jpeg)

a cost of  $\tilde{W}$  in the total cost  $W$  for SBBU, for the usual DFS recursive implementation of BP, it may contribute much more to the BP total cost because such subproblem may have to be solved several times in the occasion of backtrackings.

# 6 Conclusion

We propose a new algorithm for the DMDGP which leverages symmetry information to find the first solution quickly. By efficiently exploiting symmetries of subproblems defined by pruning edges, and reducing the degrees of freedom by taking into account already solved subproblems, the resulting algorithm appears to be quite efficient in sparse DMDGP instances, sometimes giving a significant speed-up with respect to the classic BP algorithm.

We reinforce that the numerical experiments in [16] showed that the classic BP has a better performance than continuous methods [14, 32] in problems of the DMDGP subclass. Thus, the new proposed algorithm can do even better and, to the best of our knowledge, figures as one of the fastest algorithms for this class of Distance Geometry problems.

In the proposed version of SBBU algorithm we consider a specific order for the pruning edges. In future works we expect to generalize SBBU in order to handle different pruning edges orderings and study the impact of these in the total cost  $W$ .

Another line for future investigation is the study of how imprecise distances affect the symmetry properties used in this paper and the extension of the proposed algorithm to tackle problems with noisy or interval distances. Some results considering interval distances are given in [2, 8, 18, 20, 35].

# Appendix: Proof of Lemma 2

Proof Since  $x(s')$  from Eq. (4) involves only partial reflections, in view of Property 3 in Remark 2,  $x(s') \in \hat{X}$ , i.e.  $\| x_u(s') - x_w(s') \| = d_{uw}, \forall \{u, w\} \in E_D$ .

Springer

---

Algorithmica (2021) 83:2400-2426

It remains to show that  $x(s')$  does not violate distance constraints associated to pruning edges in  $P^{ij}$ . Since the reflections are applied to positions  $x_{\ell}$  such that  $\ell \geq i + K + 1$ , pruning edges  $\{u, w\} \in P^{ij}$  with  $u &lt; w \leq i + K$  are not affected. Thus, assume that  $i + K + 1 \leq w \leq n$ . If  $u \leq i$ , then for  $\ell = i + K + 1, \ldots, w$  there exists  $\{u, w\}$  such that  $u + K + 1 \leq \ell \leq w$ , which implies that  $v_{i + K + 1}, \ldots, v_w \notin S^{ij}$ , meaning that the first symmetry vertex  $v_{\ell}$  in  $S^{ij}$  is such that  $\ell \geq w + 1$ . Thus, according to (2) and (4), partial reflections are not applied to  $x_{i + K + 1}, \ldots, x_w$ , i.e.  $x_{\ell}(s') = x_{\ell}(s)$ , for  $\ell = i + K + 1, \ldots, w$  and  $\| x_u(s') - x_w(s') \| = d_{uw}$  holds. Otherwise, for  $u \geq i + 1$ , we have that  $v_{u + K + 1}, \ldots, v_w \notin S^{ij}$ , and from (4) and (2), positions  $x_{\ell}, \ldots, x_{u + K + 1}, \ldots, x_w$  are updated by reflections  $R_s^\ell(x_\ell), \ldots, R_s^\ell(x_{u + K + 1}), \ldots, R_s^\ell(x_w)$ , for  $i + K + 1 \leq \ell \leq u + K$  such that  $v_\ell \in S^{ij}$ . Since either  $u \leq \ell - 1$ , i.e.  $x_u$  is in the hyperplane associated to  $v_\ell$ , or  $u \geq \ell$ , i.e.  $x_u$  comes after this hyperplane, in view of Remark 2, Property 2, these reflections are such that  $\| x_u(s') - x_v(s') \| = d_{uw}$ .

Acknowledgements We would like to thank Prof. Luiz M. Carvalho for valuable discussions. This work was partly supported by: (a) the Brazilian research agencies CNPq, CAPES, and FAPESP; (b) the French research agency ANR under grant ANR-19-CE45-0019 "multiBioStruct"; (c) the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie Grant Agreement No. 764759 ETN "MINOA". Part of this work was done during the visit of DG to LL at École Polytechnique, supported by CAPES/Print Process 88887.465828/2019-00.

# Declarations

Conflict of interest The authors declare that they have no conflict of interest.

# References

1. Alencar, J., Lavor, C., Liberti, L.: Realizing Euclidean distance matrices by sphere intersection. Discret. Appl. Math. 256, 5-10 (2019)
2. Baez-Sanchez, A., Lavor, C.: On the estimation of unknown distances for a class of Euclidean distance matrix completion problems with interval data. Linear Algebra Appl. 592, 287-305 (2020)
3. Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., Shindyalov, I.N., Bourne, P.E.: The protein data bank. Nucl. Acids Res. 28, 235-242 (2000)
4. Billinge, S., Duxbury, P., Gonçalves, D., Lavor, C., Mucherino, A.: Assigned and unassigned distance geometry: applications to biological molecules and nanostructures. 4OR 14, 337-376 (2016)
5. Billinge, S., Duxbury, P., Gonçalves, D., Lavor, C., Mucherino, A.: Recent results on assigned and unassigned distance geometry with applications to protein molecules and nanostructures. Ann. Oper. Res. 271, 161-203 (2018)
6. Camargo, V.S., Castelani, E.V., Fernandes, L.A.F., Fidalgo, F.: Geometric algebra to describe the exact discretizable molecular distance geometry problem for an arbitrary dimension. Adv. Appl. Clifford Algebras 29(75), (2019)
7. Cassioli, A., Günlük, O., Lavor, C., Liberti, L.: Discretization vertex orders in distance geometry. Discret. Appl. Math. 197, 27-41 (2015)
8. Costa, T., Bouwmeester, H., Lodwick, W., Lavor, C.: Calculating the possible conformations arising from uncertainty in the molecular distance geometry problem using constraint interval analysis. Inf. Sci. 415, 41-52 (2017)
9. Crippen, G., Havel, T.: Distance Geometry and Molecular Conformation. Wiley, New York (1988)
10. Dokmanic, I., Parhizkar, R., Ranieri, J., Vetterli, M.: Euclidean distance matrices: essential theory, algorithms, and applications. IEEE Signal Proc. Magaz. 32(6), 12-30 (2015)

Springer

---

Algorithmica (2021) 83:2400-2426

11. Dong, Q., Wu, Z.: A geometric build-up algorithm for solving the molecular distance geometry problem with sparse distance data. J. Global Optim. 26(3), 321–333 (2003)
12. Fidalgo, F., Gonçalves, D., Lavor, C., Liberti, L., Mucherino, A.: A symmetry-based splitting strategy for discretizable distance geometry problems. J. Global Optim. 71, 717–733 (2018)
13. Gramacho, W., Mucherino, A., Lavor, C., Maculan, N.: A parallel BP algorithm for the discretizable distance geometry problem. In: Proceedings of the Workshop on Parallel Computing and Optimization, pp. 1756–1762. IEEE, Piscataway (2012)
14. Krislock, N., Wolkowicz, H.: Explicit sensor network localization using semidefinite representations and facial reductions. SIAM J. Optim. 20, 2679–2708 (2010)
15. Lavor, C., Liberti, L., Donald, B., Worley, B., Bardiaux, B., Malliavin, T., Nilges, M.: Minimal NMR distance information for rigidity of protein graphs. Discret. Appl. Math. 256, 91–104 (2019)
16. Lavor, C., Liberti, L., Maculan, N., Mucherino, A.: The discretizable molecular distance geometry problem. Comput. Optim. Appl. 52, 115–146 (2012)
17. Lavor, C., Liberti, L., Maculan, N., Mucherino, A.: Recent advances on the discretizable molecular distance geometry problem. Eur. J. Oper. Res. 219, 698–706 (2012)
18. Lavor, C., Liberti, L., Mucherino, A.: The interval BP algorithm for the discretizable molecular distance geometry problem with interval data. J. Global Optim. 56, 855–871 (2013)
19. Lavor, C., Oliveira, A., Rocha, W., Souza, M.: On the optimality of finding DMDGP symmetries. Comput. Appl. Math. 40, 98–107 (2021)
20. Lavor, C., Souza, M., Carvalho, L., Gonçalves, D., Mucherino, A.: Improving the sampling process in the interval branch-and-prune algorithm for the discretizable molecular distance geometry. Appl. Math. Comput. 389, 125586–125597 (2021)
21. Lavor, C., Souza, M., Carvalho, L.M., Liberti, L.: On the polynomiality of finding  $^K$ DMDGP reorders. Discret. Appl. Math. 267, 190–194 (2019)
22. Liberti, L.: Distance geometry and data science. TOP 28, 271–339 (2020)
23. Liberti, L., Lavor, C.: Euclidean Distance Geometry: An Introduction. Springer, New York (2017)
24. Liberti, L., Lavor, C., Alencar, J., Abud, G.: Counting the number of solutions of  $^k$ DMDGP instances. In: Nielsen, F., Barbaresco, F. (eds.) Geometric Science of Information. Lecture Notes in Computer Science, vol. 8085, pp. 224–230. Springer, Berlin (2013)
25. Liberti, L., Lavor, C., Maculan, N.: A branch-and-prune algorithm for the molecular distance geometry problem. Int. Trans. Oper. Res. 15, 1–17 (2008)
26. Liberti, L., Lavor, C., Maculan, N., Mucherino, A.: Euclidean distance geometry and applications. SIAM Rev. 56, 3–69 (2014)
27. Liberti, L., Lavor, C., Mucherino, A.: The Discretizable Molecular Distance Geometry Problem Seems Easier on Proteins. In: Mucherino, A., Lavor, C., Liberti, L., Maculan, N. (eds.) Distance Geometry, pp. 47–60. Springer, New York (2013)
28. Liberti, L., Masson, B., Lee, J., Lavor, C., Mucherino, A.: On the number of solutions of the discretizable molecular distance geometry problem. In: Combinatorial Optimization. Constraints and Applications (COCOA11), LNCS, vol. 6831, pp. 322–342. Springer, New York (2011)
29. Liberti, L., Masson, B., Lee, J., Lavor, C., Mucherino, A.: On the number of realizations of certain Henneberg graphs arising in protein conformation. Discret. Appl. Math. 165, 213–232 (2014)
30. Maioli, D., Lavor, C., Gonçalves, D.S.: A note on computing the intersection of spheres in  $\mathbb{R}^n$ . ANZIAM J. 59(2), 271-279 (2017)
31. Malliavin, T., Mucherino, A., Lavor, C., Liberti, L.: Systematic exploration of protein conformational space using a distance geometry approach. J. Chem. Inf. Model. 59, 4486–4503 (2019)
32. Moré, J.J., Wu, Z.: Distance geometry optimization for protein structures. J. Global Optim. 15, 219–234 (1999)
33. Mucherino, A., Lavor, C., Liberti, L.: Exploiting symmetry properties of the discretizable molecular distance geometry problem. J. Bioinform. Comput. Biol. 10(3), 1242009(1-15) (2012)
34. Mucherino, A., Lavor, C., Liberti, L., Maculan, N. (eds.): Distance Geometry: Theory, Methods, and Applications. Springer, Berlin (2013)
35. Neto, L.S., Lavor, C., Lodwick, W.: A note on the Cayley-Menger determinant and the molecular distance geometry problem. Inf. Sci. 559, 1-7 (2021)
36. Newman, M.E., Ziff, R.M.: Fast Monte Carlo algorithm for site or bond percolation. Phys. Rev. E 64(1), 016706 (2001)
37. Nucci, P., Nogueira, L., Lavor, C.: Solving the discretizable molecular distance geometry problem by multiple realization trees. In: Mucherino et al. [34], pp. 161-176

Springer

---

Algorithmica (2021) 83:2400-2426

38. Saxe, J.B.: Embeddability of weighted graphs in $k$-space is strongly NP-hard. In: Proceedings of $17^{th}$ Allerton Conference in Communications, Control and Computing, pp. 480-489. Monticello, IL (1979)

39. Wütrich, K.: Protein structure determination in solution by nuclear magnetic resonance spectroscopy. Science 243, 45-50 (1989)

Publisher's Note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.

Springer
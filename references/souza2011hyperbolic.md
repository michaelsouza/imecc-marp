Operations Research Letters 39 (2011) 461-465

ELSEVIER

Contents lists available at SciVerse ScienceDirect

Operations Research Letters

journal homepage: www.elsevier.com/locate/orl

Operations Research Letters

# Hyperbolic smoothing and penalty techniques applied to molecular structure determination

Michael Souza $^{a,\ast}$, Adilson Elias Xavier $^{a}$, Carlile Lavor $^{b}$, Nelson Maculan $^{a}$

$^{a}$ Federal University of Rio de Janeiro (COPPE/UFRJ), 21945-970, Rio de Janeiro - RJ, Brazil
$^{b}$ State University of Campinas (IMECC-UNICAMP), 13081-970, Campinas - SP, Brazil

## ARTICLE INFO

Article history:
Received 29 March 2010
Accepted 19 July 2011
Available online 29 July 2011

Keywords:
Protein structure
NMR experiments
Global optimization
Smoothing techniques
Molecular distance geometry problem

## ABSTRACT

This work considers the problem of estimating the relative positions of all atoms of a protein, given a subset of all the pair-wise distances between the atoms. This problem is NP-hard, and the usual formulations are nonsmoothed and nonconvex, having a high number of local minima. Our contribution is an efficient method that combines the hyperbolic smoothing and the penalty techniques that are useful in obtaining differentiability and reducing the number of local minima.

© 2011 Elsevier B.V. All rights reserved.

## 1. Introduction

The problem of determining molecular structures has attracted great interest due to its application in relevant areas such as medicine, pharmacy, biology, design of materials, and chemistry [25]. Until 1984, the X-ray crystallography technique used to be the ultimate tool for obtaining information about protein structure at atomic resolution. The introduction of nuclear magnetic resonance (NMR) as a technique for protein structure determination made it possible to obtain structures with high precision in an environment (solution) that is much closer to the natural surroundings of a living organism than the crystals used in crystallography [25].

NMR experiments yield valuable information: a network of bounds for distances involving pairs of atoms spatially close, at distances less than $5 - 6\AA$ ($1\AA = 10^{-10}\mathrm{m}$). In the molecular distance geometry problem (MDGP), the challenge is to obtain a valid three-dimensional configuration for the atoms, i.e., a configuration that is compatible with the network of bounds obtained in the NMR experiment. The MDGP can be precisely defined as

$$
\text{determine } x_1, x_2, \dots, x_m
$$

$$
\text{subject to } l_{ij} \leq \| x_i - x_j \| \leq u_{ij}, \quad \forall (i,j) \in K,
$$

$$
x_i \in \mathbb{R}^3, \; i = 1, \dots, m,
$$

where $K \subset \{1, \ldots, m\} \times \{1, \ldots, m\}$ is a set that identifies the available lower $(l_{ij})$ and upper $(u_{ij})$ bounds for the pair-wise Euclidean distances, and $m$ is the number of atoms in the protein.

The MDGP can be reformulated as a mathematical programming problem.

## Definition 1 (MDGP)

$$
(P) \min f(x) = \sum_{(i,j) \in K} \max \{l_{ij} - \| x_i - x_j \|, 0\}
$$

$$
+ \max \{\| x_i - x_j \| - u_{ij}, 0\}
$$

$$
\text{where } x = (x_1, \dots, x_m) \in \mathbb{R}^{3m},
$$

$$
x_k \in \mathbb{R}^3, \; \forall k \in \{1, \dots, m\}. \tag{1}
$$

It is easy to see that $f(x_1, \ldots, x_m) = 0$ if and only if all the restrictions $l_{ij} \leq \| x_i - x_j \| \leq u_{ij}$ are satisfied.

The contribution of this work is a new algorithm for the MDGP. Our proposal combines the hyperbolic smoothing and penalty techniques [33]. This approach allows the application of classical optimization methods by introducing differentiability into formulation (1) and, more importantly, it generates trajectories that have a great chance to converge to the global minimum of $(P)$.

In Section 2, we point out the complexity of the MDGP and the different approaches found in the literature. In Section 3, we propose a new algorithm to solve the MDGP by using hyperbolic functions. Section 4 presents computational experiments carried out with instances generated from real proteins. Finally, in Section 5, we summarize the contributions of this work.

0167-6377/$ - see front matter © 2011 Elsevier B.V. All rights reserved.
doi:10.1016/j.orl.2011.07.007

---

M. Souza et al./ Operations Research Letters 39 (2011) 461-465

# 2. The MDGP

The MDPG can be posed as a nonlinear global optimization problem, such as $(P)$, for instance. Unfortunately, the problem of finding the ideal formulation for the MDGP remains unsolved [7,22,24]. There are no results of convexity. Besides, in general, the formulations using penalty functions are not differentiable and have many local minima, which precludes the direct application of classical and more robust optimization methods. The nondifferentiability arises from the inclusion of Euclidean norm in the definition of the constraints and also from the use of the $\max\{\cdot, 0\}$ as the penalty function. The large number of local minima is connected to the strong combinatorial appeal of the MDGP, which gets worse with the increasing number of atoms considered [8,15,18].

Particular cases of the MDGP are solved in a relatively easy way. For instance, when we know all distances $d_{ij} = \| x_i - x_j\|$, i.e., $d_{ij} = l_{ij} = u_{ij}$ and $K = \{1,2,\dots,m\}^2$, a solution can be obtained by factoring the matrix $D = \{d_{ij}\}$. In fact, a very simple algorithm presented by Dong and Wu generates a solution to the MDGP in linear time, in the case where all distances are known [8].

However, in practice, the NMR experiments just provide the subset of the distances between atoms that are spatially close [25], and the data accuracy is limited. Thus, in a real scenario, the set $K$ is sparse, and $l_{ij} &lt; u_{ij}$. This form of the MDGP seems to be easier to solve, since the constraints are less restrictive. However, in practice, the lower and upper limits are close, and thus the problem is still difficult to solve. In [21], Moré and Wu showed that if the upper and lower bounds are close, then the MDGP with relaxed distances belongs to the NP-hard class.

Different approaches to the MDGP have been explored: smoothing techniques [19,22], semidefinite programming [5], differences of convex functions [3], geometric inequalities [7], graph rigidity [15], geometric build up algorithms [6,31,32], alternating projections [13], multi-scaling [16,27], stochastic perturbation [37] and branch-and-prune [18] (for a review of the MDGP methods, see [17]).

# 3. The hyperbolic functions

The most common penalty function used in the MDGP is the function

$$
p _ {\lambda} (y) = \lambda \max  \{y, 0 \},
$$

where $\lambda$ is a parameter that indicates the penalty intensity. We alternatively propose the use of the penalty function

$$
\phi_ {\lambda , \tau} (y) = \lambda y + \sqrt {\lambda^ {2} y ^ {2} + \tau^ {2}},
$$

where $\lambda &gt; 0$, $\tau &gt; 0$ and $y \in \mathbb{R}$.

Note that $\lim_{\tau \to 0} \phi_{1/2,\tau}(y) = \max \{y, 0\}$. In other words, if $\lambda = 1/2$ and $\tau$ is small, the function $\phi_{\lambda,\tau}$ is a good approximation to $\max \{\cdot, 0\}$. In Fig. 1, we see that the parameter $\tau$ controls the smoothness degree and $\lambda$ controls the intensity (weight) of the penalty function $\phi_{\lambda,\tau}$.

Finally, in order to obtain the differentiability in problem (1), we propose to replace the function $\| \cdot \|$ by the hyperbolic smoothing function

$$
\theta_ {\tau} \left(x _ {i}\right) = \sqrt {\tau^ {2} + \sum_ {k = 1} ^ {3} x _ {i k} ^ {2}},
$$

where $x_{i} = (x_{i1},x_{i2},x_{i3})\in \mathbb{R}^{3}$, $\tau &gt;0$. It needs to be pointed out that the hyperbolic smoothing $\theta_{\tau}$ of the Euclidean norm is a natural approach that has been frequently used not only in location theory [11,12,30] but also in other fields such as geomechanics [1,36], PDE theory [9,29], and packing problems [35].

![img-0.jpeg](img-0.jpeg)
Fig. 1. In the graph of the penalty function $\phi_{\lambda,\tau}$, the parameter $\lambda$ controls the intensity (slope) of penalty, and $\tau$ controls the smoothness.

![img-1.jpeg](img-1.jpeg)
Fig. 2. The graph of the hyperbolic smoothing $\theta_{\tau}(x)$ is an equilateral hyperbola.

In Fig. 2, we can see clearly that the maximum distance between the original function and the smooth function $\theta_{\tau}$ is equal to the parameter $\tau$. The graphs of $\theta_{\tau}$ and $\phi_{\lambda,\tau}$ are equilateral hyperbolas, which motivates the inclusion of the hyperbolic term used in the nomenclature of these functions.

Inserting the functions $\theta_{\tau}^{ij}(x) = \theta_{\tau}(x_i - x_j)$ and $\phi_{\lambda,\tau}$ in problem (1), we obtain the smoothed problem:

$$
\left(P _ {\lambda , \tau}\right) \min  f _ {\lambda , \tau} (x) = \sum_ {(i, j) \in K} \phi_ {\lambda , \tau} \left(l _ {i j} - \theta_ {\tau} ^ {i j} (x)\right) + \phi_ {\lambda , \tau} \left(\theta_ {\tau} ^ {i j} (x) - u _ {i j}\right),
$$

where $K\subset \{1,\ldots ,m\} ^2,x = (x_1,\dots,x_m),x_k\in \mathbb{R}^3,$

$$
\forall k = 1, \dots , m. \tag {2}
$$

Problem $(P_{\lambda,\tau})$ is infinitely differentiable with respect to $x$ and therefore allows the application of classical optimization methods. However, it should be noted that problem $(P_{\lambda,\tau})$ is not exactly equal to problem (1). Our approach to solving $(P)$ is similar to that of the homotopy continuation methods, where the idea is to start by solving the easiest problem $(P_{\lambda,\tau})$, and then the most difficult one $(P)$ [2,28]. We propose to solve an infinite sequence of smooth problems $(P_{\lambda,\tau_k})$ that are parameterized by a descending sequence of parameters $\tau_k$, $k = 1, 2, \ldots$, tending to zero, i.e., $\tau^{k+1} &lt; \tau_k$ with $\lim_{k \to +\infty} \tau_k = 0$. Through this procedure, the sequence of smooth problems $(P_{\lambda,\tau_k})$ gradually converges to the original problem.

Note that in the definition of the sequence $(P_{\lambda,\tau_k})$, we use the same parameter $\lambda$ in all terms of the objective function. This simply reflects our choice about not privileging any particular constraint. However, in an application with data from different sources, we may wish to influence the order in which the constraints are satisfied, and this could be done by setting different values for the parameter $\lambda$ in each constraint.

---

M. Souza et al. / Operations Research Letters 39 (2011) 461-465

## 3.1. Convex properties

Note that problem (2) is nonconvex, and the classical methods might stagnate at local optima. We claim that the use of hyperbolic smoothing and penalty functions is not justified solely by the introduction of differentiability. The great advantage in using these functions is related to the convexity, as originally pointed out by Xavier in the context of the MDGP in [34] and explored by Macambira in [23]. The following results prove that problem $(P_{\lambda,\tau})$ can be made convex by taking a smoothing parameter value $\tau$ sufficiently large.

**Proposition 1.** The function $\theta_{\tau}^{ij} : \mathbb{R}^{3m} \to \mathbb{R}$, given by $\theta_{\tau}^{ij}(x) = \theta_{\tau}(x_i - x_j)$, is convex for all $\tau &gt; 0$.

**Proof.** It can be checked directly that the Hessian matrix of $\theta_{\tau}^{ij}$ has the following structure:

$$
\nabla^ {2} \theta_ {\tau} ^ {i j} (x) = \left( \begin{array}{c c c c c} 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\ 0 &amp; H _ {\tau} ^ {i j} (x) &amp; 0 &amp; - H _ {\tau} ^ {i j} (x) &amp; 0 \\ 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\ 0 &amp; - H _ {\tau} ^ {i j} (x) &amp; 0 &amp; H _ {\tau} ^ {i j} (x) &amp; 0 \\ 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \end{array} \right) _ {3 m \times 3 m},
$$

where

$$
H _ {\tau} ^ {i j} (x) = \frac {l _ {3 \times 3}}{\theta_ {\tau} ^ {i j} (x)} - \frac {(x _ {i} - x _ {j}) (x _ {i} - x _ {j}) ^ {T}}{(\theta_ {\tau} ^ {i j} (x)) ^ {3}},
$$

$l_{3\times 3}$ is the identity matrix, and the zeros in the matrix $\nabla^2\theta_\tau^{ij}(x)$ are in fact zero block matrices.

For any $z \neq 0$ in $\mathbb{R}^3$, it follows that

$$
z ^ {T} H _ {\tau} ^ {i j} (x) z \geq \frac {\| z \| ^ {2}}{\theta_ {\tau} ^ {i j} (x)} - \frac {\| z \| ^ {2} \| x _ {i} - x _ {j} \| ^ {2}}{(\theta_ {\tau} ^ {i j} (x)) ^ {3}} &gt; 0.
$$

Thus, we have $H_{\tau}^{ij}(x) \succeq 0$.

Finally, for any $v = (v_{1}, v_{2}, \ldots, v_{m}) \neq 0$ in $\mathbb{R}^{3m}$,

$$
v ^ {T} \nabla^ {2} \theta_ {\tau} ^ {i j} (x) v = \left(v _ {i} - v _ {j}\right) ^ {T} H _ {\tau} ^ {i j} (x) \left(v _ {i} - v _ {j}\right) \geq 0.
$$

Therefore, $\nabla^2\theta_\tau^{ij}(x) \succeq 0$ and, consequently, $\theta_\tau^{ij}$ is a convex function.

**Corollary 1.** The objective function $f_{\lambda,\tau}$ given in (2) is convex for all $\tau &gt; \max \{u_{ij} : (i,j) \in K\}$.

**Proof.** By definition, each term of the objective function of problem $(P_{\lambda,\tau})$ is given by

$$
\begin{array}{l} f _ {\lambda , \tau} ^ {i j} (x) = \phi_ {\lambda , \tau} \left(l _ {i j} - \theta_ {\tau} ^ {i j} (x)\right) + \phi_ {\lambda , \tau} \left(\theta_ {\tau} ^ {i j} (x) - u _ {i j}\right) \\ = \left(l _ {i j} - u _ {i j}\right) + \sqrt {\lambda^ {2} \left(\theta_ {\tau} ^ {i j} (x) - l _ {i j}\right) ^ {2} + \tau^ {2}} \\ + \sqrt {\lambda^ {2} \left(\theta_ {\tau} ^ {i j} (x) - u _ {i j}\right) ^ {2} + \tau^ {2}}. \\ \end{array}
$$

Note that we can rewrite $f_{\lambda,\tau}^{ij}(x)$ as

$$
f _ {\lambda , \tau} ^ {i j} (x) = \left(l _ {i j} - u _ {i j}\right) + h _ {\lambda , \tau} ^ {i j} \left(\theta_ {\tau} ^ {i j} (x), l _ {i j}\right) + h _ {\lambda , \tau} ^ {i j} \left(\theta_ {\tau} ^ {i j} (x), u _ {i j}\right),
$$

where $h_{\lambda,\tau}^{ij}(x,\alpha) = \sqrt{\lambda^2(\theta_\tau^{ij}(x) - \alpha)^2 + \tau^2}$.

First, we shall prove that the function $h_{\lambda,\tau}^{ij}(x,\alpha)$ is convex if $\alpha \leq \tau$. In fact, if $\alpha \leq \tau$, we have

$$
\begin{array}{l} \nabla_ {\tau} ^ {2} h _ {\lambda , \tau} ^ {i j} (x, \alpha) = \frac {\lambda^ {2} \tau^ {2} \nabla \theta_ {\tau} ^ {i j} (x) \nabla^ {\prime} \theta_ {\tau} ^ {i j} (x)}{h _ {\lambda , \tau} ^ {2} (x , \alpha)} + \frac {(\theta_ {\tau} ^ {i j} (x) - \alpha) \nabla^ {2} \theta_ {\tau} ^ {i j} (x)}{h _ {\lambda , \tau} (x , \alpha)} \\ \succeq 0, \\ \end{array}
$$

because $\alpha \leq \tau \leq \theta_{\tau}^{ij}(x), \forall x \in \mathbb{R}^{3m}$, and, by Proposition 1, $\nabla^2\theta_{\tau}^{ij}(x) \succeq 0$.

Finally, by hypothesis, $\tau &gt; \max \{u_{ij} : (i,j) \in K\} \geq \max \{l_{ij} : (i,j) \in K\}$; thus each term of the objective function of problem $(P_{\lambda,\tau})$ is a convex function, which concludes the proof.

In practice, the convexification of $(P_{\lambda,\tau})$ obtained by choosing high values for the parameter $\tau$ should be avoided, because the objective function $f_{\lambda,\tau}$ tends to the constant $2\tau\sqrt{\lambda^2 + 1}$ as $\tau$ grows. Thus, the parameter $\tau$ should express a harmonious compromise between problem simplification and solution consistency. That is, we should make the terms $f_{\lambda,\tau}^{ij}$ convex by choosing an appropriate value for the parameter $\tau$, but avoiding excessively high values.

The properties of the hyperbolic smoothing and penalty techniques are used in Algorithm 1 for solving problem $(P)$, Definition 1. Following the notation proposed in [22], we denote by $\mathbf{locmin}(f, x, \mathcal{M})$ the minimizer of the function $f$ generated by a local minimization algorithm $\mathcal{M}$, taking $x$ as the starting point.

Despite the theory presented, it needs to be pointed out that the routine $sph$ is heuristic. This means that there is no guarantee that a global minimum will be reached for all possible starting points. Part of this uncertainty is strongly related to the lack of a precise way of setting the smoothing parameter $\tau$ and also to the application of the local minimization procedures at each step.

**Algorithm 1 Routine sph**

```
sph(x,λ,τ,ρ) : ρ ∈ (0,1);
while(fλ,τ(x) &gt; εf) do
x = locmin(fλ,τ, x, M);
τ = ρτ;
end do
return (x, fλ,τ(x));
```

## 4. Computational experiments

Our procedure $sph$ was implemented in $C$, except for local minimization routine $v_035$, encoded in FORTRAN, and available from the Harwell Subroutine Library. The routine $v_035$ implements the method BFGS with limited memory [20] (For additional information on this routine, see http://www.hsl.rl.ac.uk/catalog.html). We used a computer with an Intel Core 2 Duo CPU T5670 1.80 GHz, 3.0 GB of RAM, and Linux OS-64 bit.

The instances considered were derived from the three-dimensional structure of the fragments made up of the first 100 and 200 atoms of the chain $A$ of the protein 1 GPV [14,26]. For each fragment, we generated a set of constraints considering only atoms in the same residue or the neighboring residues. Formally,

$$
K = \{(i, j): x _ {i} \in R (k), x _ {j} \in R (k) \cup R (k + 1) \},
$$

where $R(k)$ represents the $k$-th residue.

As Moré and Wu in [22], we consider as a solution any set of coordinates $x \in \mathbb{R}^{3m}$ such that

$$
(1 - \tau_ {d}) l _ {i j} \leq \| x _ {i} - x _ {j} \| \leq (1 + \tau_ {d}) u _ {i j}, \quad (i, j) \in K, \tag {3}
$$

where $\tau_d = 10^{-2}$. This tolerance reflects the precision for the bond lengths [10,22].

The bounds $l_{ij}$ and $u_{ij}$ were given by the equations

$$
l _ {i j} = (1 - \varepsilon_ {b}) \| \hat {x} _ {i} - \hat {x} _ {j} \|, \quad u _ {i j} = (1 + \varepsilon_ {b}) \| \hat {x} _ {i} - \hat {x} _ {j} \|,
$$

where $\hat{x} = (\hat{x}_1, \hat{x}_2, \dots, \hat{x}_m)$ represents the known structure of protein 1GP available in the PDB (Protein Data Bank [4]). We generated four instances for each fragment by taking $\epsilon_b$ equal to 0.04, 0.08, 0.12, and 0.16.

The parameters of the routine $sph$ were set to $\lambda = 0.5$ and $\rho = 0.99$, and the initial value of the smoothing parameter $\tau$ was defined at run time as the third quartile of the set $\{(l_{ij} + u_{ij}) / 2 :$

---

M. Souza et al./ Operations Research Letters 39 (2011) 461-465

Table 1 The number of solutions (ns) obtained by the routines va35, dgsol and sph and the average CPU time (cput) in seconds.

|  ε | va35 |   | dgsol |   | sph  |   |
| --- | --- | --- | --- | --- | --- | --- |
|   |  ns | cput | ns | cput | ns | cput  |
|  100 atoms  |   |   |   |   |   |   |
|  0.04 | 0 | 0.26 | 73 | 7.33 | 100 | 1.58  |
|  0.08 | 0 | 0.26 | 72 | 8.13 | 100 | 1.46  |
|  0.12 | 3 | 0.26 | 100 | 8.17 | 100 | 1.40  |
|  0.16 | 27 | 0.23 | 100 | 9.28 | 100 | 1.32  |
|  200 atoms  |   |   |   |   |   |   |
|  0.04 | 0 | 0.54 | 44 | 29.25 | 83 | 8.48  |
|  0.08 | 0 | 0.54 | 36 | 30.36 | 88 | 7.06  |
|  0.12 | 0 | 0.54 | 100 | 32.57 | 100 | 6.50  |
|  0.16 | 7 | 0.51 | 100 | 33.03 | 100 | 6.01  |

$(i,j)\in K\}$ . With this choice of the value of the parameter  $\tau$ , by Corollary 1, about  $75\%$  of the terms of the objective function of the initial problem  $(P_{\lambda ,\tau})$  are convex. The stopped parameter  $\lambda_{j}$  of  $sph$ , Algorithm 1, were set to  $\lambda_{j} = \lambda_{d} = 10^{-2}$ .

The starting points of our experiments were generated by the procedure struct proposed by Moré and Wu in [22], which is described in Algorithm 2 (the source code of the routine struct is available at http://www.mcs.anl.gov/~more/dgsol).

Algorithm 2 Routine struct
```txt
struct(K) : K ⊂ {1,2,...,m}^2, x ∈ R^3m;
set x = (x1,...,xm) = 0;
L = {1,2,...,m};
while(L is not empty) do
choose i ∈ L;
set Mi = {j: (i, j) ∈ K, j ∈ L};
for each j ∈ Mi
generate y ∈ R^3 such that ||xi - y|| = δij;
set xj = y;
end for
remove i from L;
end do
return x = (x1,...,xm);
```

This routine generates starting points  $x = (x_{1},\ldots ,x_{m})$  that satisfy at least  $m - 1$  constraints  $\| x_{i} - x_{j}\| = \delta_{ij}$ , where  $m$  is the number of atoms and  $\delta_{ij} = (l_{ij} + u_{ij}) / 2$ .

In the first experiment, we compared the performance of the routine sph with va35 and with the routine dgsol proposed by Moré and Wu in [22]. We decided to compare our approach with that of the dgsol, whose importance is justified in the MDGP literature, because it also uses a smoothing technique and the instances used to test the algorithm are clearly explained in [22]. Besides that, its code is available on the Internet (http://www.mcs.anl.gov/~more/dgsol). Therefore, a fair comparison between the methods could be made.

In each test, we generated 100 random points using the routine struct. The three routines sph, va35, and dgsol were initialized with these 100 random starting points. The results appear in Table 1, where ns and cput represent the number of solutions and the average CPU time computed in seconds, respectively.

Table 1 shows that the routine sph found more global minima than the routines dgsol and va35. The routine sph also demanded lower computational effort (CPU time) than the one required by the routine dgsol, where multidimensional integrals need to be evaluated in the calculation of the objective function and its derivatives (see [22] for more details).

In the second experiment, we verified computationally the relation between the number of local minimizers of the smoothed function  $f_{\lambda,\tau}$  and the value of the parameter  $\tau$ . Following the

![img-2.jpeg](img-2.jpeg)
Fig. 3. The number of local minimizers of  $f_{\lambda, \tau}$  as a function of  $\tau$  for  $\varepsilon = 0.04$  and 100 different starting points.

strategy adopted by Moré and Wu, we estimated the number of local minimizers of  $f_{\lambda,\tau}$  for each  $\tau$  by taking 100 random starting points and by minimizing  $f_{\lambda,\tau}$  using the routine va35.

The number of distinct local minimizers found by  $va35$  is plotted in Fig. 3. For these results, minimizers  $u$  and  $v$  of  $f_{\lambda,\tau}$  are declared to be the same if

```txt
$|f_{\lambda ,\tau}(u) - f_{\lambda ,\tau}(v)|\leq \tau_{\tau}\max \{|f_{\lambda ,\tau}(u)|,|f_{\lambda ,\tau}(v)|\} ,$  where  $\tau_t = 10^{-6}$  or if
max  $\{f_{\lambda ,\tau}(u),f_{\lambda ,\tau}(v)\} \leq \tau_o$  where  $\tau_o = 10^{-2}$
```

The results in Fig. 3 show that the number of local minimizers of  $f_{\lambda, \tau}$  decreases as  $\tau$  increases, in accordance with the theory presented.

# 5. Conclusion and proposals for future work

The problem of determining molecular structures of proteins has attracted great interest due to its application in relevant areas such as medicine, pharmacy, biology, design of materials, and chemistry. The alternatives that have been used in this endeavor are experiments involving nuclear magnetic resonance (NMR), where the distances between some of the atoms of the protein can be estimated. In the molecular distance geometry problem (MDGP), the goal is to determine the three-dimensional conformation of a protein, based on estimates of the distances provided by NMR experiments.

We have presented an algorithm for the MDGP in its most realistic version, i.e., with a sparse set of inaccurate distances. Our proposal, the algorithm sph, starts with a formulation of least squares and tries to convexify the associated objective function through hyperbolic smoothing and penalty functions. Our approach was based on the proved property which states that appropriate values for the parameters of our smoothing function generate a convex problem.

We performed computational experiments with real data obtained from the PDB. The results indicate that the algorithm sph is more efficient than the multistart approach and the Gaussian smoothing technique proposed in the literature.

As a proposal for future research, we are interested in improving the MDGP formulation. In this work, we consider only distance constraints. However, by inserting restrictions arising from theoretical models of protein conformation, we can, in theory, get more significant results. Another point that deserves attention is the use of more realistic instances. For that purpose, we intend to establish partnerships with laboratories that produce real data from NMR spectroscopy.

---

M. Souza et al. / Operations Research Letters 39 (2011) 461-465

# Acknowledgments

The authors would like to thank the Brazilian research agencies FAPESP and CNPq for financial support and the referee for the valuable comments.

# References

[1] A.J. Abbo, S.W. Sloan, A smooth hyperbolic approximation to the Mohr-Coulomb yield criterion, Computers &amp; Structures 54 (1995) 427-441.
[2] E. Allgower, K. Georg, Numerical Continuation Methods: An Introduction, Springer-Verlag, New York, 1990.
[3] L.T.H. An, P.D. Tao, Large-scale molecular optimization from distance matrices by d. c. optimization approach, SIAM Journal on Optimization 14 (2003) 77-114.
[4] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne, The protein data bank, Nucleic Acids Research 28 (2000) 235-242.
[5] P. Biswas, K.C. Toh, Y. Ye, A distributed sdp approach for large-scale noisy anchor-free graph realization with applications to molecular conformation, SIAM Journal on Scientific Computing 30 (2008) 1251-1277.
[6] R.S. Carvalho, C. Lavor, F. Protti, Extending the geometric buildup algorithm for the molecular distance geometry problem, Information Processing Letters 108 (2008) 234-237.
[7] G.M. Crippen, T.F. Havel, Distance Geometry and Molecular Conformation, John Wiley &amp; Sons, 1988.
[8] Q. Dong, Z. Wu, A linear-time algorithm for solving the molecular distance geometry problem with exact inter-atomic distances, Journal of Global Optimization 22 (2002) 365-375.
[9] L. Eldén, Hyperbolic approximation for a Cauchy problem for the heat equation, Inverse Problems 4 (1988) 59-70.
[10] R.A. Engh, R. Huber, Accurate angles and bond parameters for X-ray protein structure refinement, Acta Crystallographica. A 47 (1991) 392-400.
[11] J.W. Eyster, J.A. White, W.W. Wierwille, On solving multifacility location problems using a hyperboloid approximation procedure, AIIE Transaction 5 (1973) 1-6.
[12] R.L. Francis, F. McGinnis, J.A. White, Facility Layout and Location: An Analytical Approach, Prentice-Hall, 1992.
[13] W. Glunt, T.L. Hayden, S. Hong, J. Wells, An alternating projection algorithm for computing the nearest euclidean distance matrix, SIAM Journal on Matrix Analysis and Applications 11 (1990) 589-600.
[14] Y. Guan, H. Zhang, R.N.H. Konings, C.W. Hilbers, T.C. Terwilliger, A.H.J. Wang, Crystal structure of y41h and y41f mutants of gene v suggest possible protein-protein interactions in the gvp-vsdna complex, Biochemistry 33 (1994) 7768.
[15] B. Hendrickson, The molecule problem: exploiting structure in global optimization, SIAM Journal on Optimization (5) (1995) 835-857.
[16] A.J. Kearsley, R.A. Tapsa, M. Trosset, The solution of the metric stress and stress problems in multidimensional scaling by newton's method, Computational Statistics 13 (1998) 369-396.
[17] C. Lavor, L. Liberti, N. Maculan, Molecular distance geometry problem, in: Encyclopedia of Optimization, 2nd edition, Springer, 2009, pp. 2305-2311.

[18] L. Liberti, C. Lavor, N. Maculan, A branch-and-prune algorithm for the molecular distance geometry problem, International Transactions in Operational Research 15 (2008) 1-17.
[19] L. Liberti, C. Lavor, N. Maculan, F. Marinelli, Double variable neighbourhood search with smoothing for the molecular distance geometry problem, Journal of Global Optimization 43 (2009) 207-218.
[20] D. Liu, Nocedal, On the limited memory BFGS method for large scale optimization, Technical Report NA-03, Department of Electrical Engineering and Computer Science Northwestern University, 1988.
[21] J.J. Moré, Z. Wu, Global continuation for distance geometry problems, SIAM Journal on Optimization 7 (1997) 814-836.
[22] J.J. Moré, Z. Wu, Distance geometry optimization for protein structures, Journal of Global Optimization 15 (1999) 219-234.
[23] A.F.U. dos Santos Macambira, Determinação de Estruturas de Proteínas via Suavização e Penalização Hyperbólica, Master's Thesis, Universidade Federal do Rio de Janeiro, 2003.
[24] J.B. Saxe, Embeddability of weighted graphs in k-space is strongly np-hard, in: Proc. 17th Allerton Conf. in Communications, Control, and Computing, pp. 480-489.
[25] T. Schlick, Molecular Modeling and Simulation: an Interdisciplinary Guide, Springer, 2002.
[26] M. Skinner, H. Zhang, D. Leschnitzer, Y. Guan, H. Bellamy, R. Sweet, C. Gray, R. Konings, A. Wang, T. Terwilliger, Structure of the gene v protein of bacteriophage f1 determined by multi-wavelength X-ray diffraction on the selenomethionyl protein, Proceedings of the National Academy of Sciences of the United States of America 91 (1994) 2071.
[27] M. Trosset, Applications of multidimensional scaling to molecular conformation, Computing Science and Statistics 29 (1998) 148-152.
[28] T. Watson Raphael, T. Layne, Modern homotopy methods in optimization, Computer Methods in Applied Mechanics and Engineering 74 (1989) 289-305.
[29] C.F. Weber, Analysis and solution of the ill-posed inverse heat conduction problem, International Journal of Heat and Mass Transfer 24 (1981) 1783-1792.
[30] G.O. Wesolowsky, R.F. Love, The optimal location of new facilities using a generalized rectangular distance weber problem, Management Science 18 (1972) 656-663.
[31] D. Wu, Z. Wu, An updated geometric build-up algorithm for solving the molecular distance geometry problems with sparse distance data, Journal of Global Optimization 37 (2007) 661-673.
[32] D. Wu, Z. Wu, Y. Yuan, Rigid versus unique determination of protein structures with geometric buildup, Optimization Letters 2 (2008) 319-331.
[33] A.E. Xavier, Hyperbolic penalty: a new method for nonlinear programming with inequalities, International Transactions in Operational Research 8 (2001) 659-671.
[34] A.E. Xavier, Convexificação do problema de distência geométrica através da técnica da suavização hyperbolica, Workshop em Biociências COPPE/UFRJ, in: Workshop em Biociências, COPPE/UFRJ, 2003.
[35] A.E. Xavier, A.A.F. de Oliveira, Optimal covering of plane domains by circles via hyperbolic smoothing, Journal of Global Optimization 31 (2005) 493-504.
[36] O.C. Zienkiewicz, G.N. Pande, Some useful forms of isotropic yield surfaces for soil and rock mechanics, Finite elements in geomechanics (1977) 179-198.
[37] Z. Zou, R.H. Bird, R.B. Schnabel, A stochastic/perturbation global optimization algorithm for distance geometry problems, Journal of Global Optimization 11 (1997) 91-105.
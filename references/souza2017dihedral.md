# Modeling the Molecular Distance Geometry Problem Using Dihedral Angles

Michael Souza
Federal University of Ceará/Brazil
Department of Applied Mathematics

Carlile Lavor
Univerisity of Campinas
Department of Applied Mathematics

![img-0.jpeg](img-0.jpeg)
UNICAMP

![img-1.jpeg](img-1.jpeg)
UNIVERSIDADE FEDERAL DO CEARÁ

---

# Definitions

---

Distance Geometry Problem

Distance Geometry (DG) is the study of geometry based on distance.

Menger, Dimension theorie. Teubner. Berlim, 1928.
Blumenthal, Theory and applic. of distance geometry. Oxford, 1953.
Yemini, The positioning problem, Proc. Dist. Sensor Networks, 1978.

The Distance Geometry Problem (DGP) is to find points

$$
x _ {1}, x _ {2}, \ldots , x _ {n} \in V,
$$

such that

$$
l _ {i j} \leq \left\| x _ {i} - x _ {j} \right\| \leq u _ {i j}, (i, j) \in E \subseteq \{1, 2, \ldots , n \} ^ {2}.
$$

---

# Distance Geometry Problem

The Distance Geometry Problem (DGP) is to find points

$$
x _ {1}, x _ {2}, \ldots , x _ {n} \in V,
$$

such that

$$
l _ {i j} \leq \left\| x _ {i} - x _ {j} \right\| \leq u _ {i j}, (i, j) \in E \subseteq \{1, 2, \ldots , n \} ^ {2}.
$$

![img-2.jpeg](img-2.jpeg)

![img-3.jpeg](img-3.jpeg)

![img-4.jpeg](img-4.jpeg)

---

# Complexity

---

# Distance Geometry Problem

Saxe, Embeddability of weighted graphs in $k$-space is strongly NP-hard. Proc. 17th Allerton Conf. in Comm., Control, and Computing, 1979.

![img-5.jpeg](img-5.jpeg)
Distance Geometry Problem

![img-6.jpeg](img-6.jpeg)

![img-7.jpeg](img-7.jpeg)
Partition Problem

$$
S = \{1, 2, 5, 6, 7, 8, 9\}
$$

$$
S_1 = \{1, 7, 8\} \quad S_2 = \{2, 5, 9\}
$$

DGP is a NP-Hard problem

---

# Applications

---

DGP - Applications

Wireless sensor networks (without GPS)
- Localization of network sensors

Savvides et al., Dynamic fine-grained loc. in ad-hoc nets of sensors, 2001.
Eren et. al, Rigidity, computation, and randomization in net loc., 2004.
Liu et al., Survey of wireless indoor positioning techniques and sys., 2007.

![img-8.jpeg](img-8.jpeg)

![img-9.jpeg](img-9.jpeg)

![img-10.jpeg](img-10.jpeg)
Images from Savvides2001

---

DGP - Applications

## Statics

- Equilibrium under the action of external forces.

- Alfakih, Graph rigidity via Euclidean distance matrices, 2000.
- Alfakih, On dimensional rigidity of bar-and-joint frameworks, 2007.

![img-11.jpeg](img-11.jpeg)
Rigidity and uniqueness of solution

![img-12.jpeg](img-12.jpeg)

---

DGP - Applications

# Dimensionality Reduction

- Multidimensional scaling (MDS)

Borg et al., Modern multidimensional scaling: Theory and appl., 2010.
Everitt et al., The analysis of proximity data, 1997.
Dzemyda et al., Multidimensional data visualization, 2012.

|  Crime | No. | 1 | 2 | 3 | 4 | 5 | 6 | 7  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  Murder | 1 | 1.00 | 0.52 | 0.34 | 0.81 | 0.28 | 0.06 | 0.11  |
|  Rape | 2 | 0.52 | 1.00 | 0.55 | 0.70 | 0.68 | 0.60 | 0.44  |
|  Robbery | 3 | 0.34 | 0.55 | 1.00 | 0.56 | 0.62 | 0.44 | 0.62  |
|  Assault | 4 | 0.81 | 0.70 | 0.56 | 1.00 | 0.52 | 0.32 | 0.33  |
|  Burglary | 5 | 0.28 | 0.68 | 0.62 | 0.52 | 1.00 | 0.80 | 0.70  |
|  Larceny | 6 | 0.06 | 0.60 | 0.44 | 0.32 | 0.80 | 1.00 | 0.55  |
|  Auto theft | 7 | 0.11 | 0.44 | 0.62 | 0.33 | 0.70 | 0.55 | 1.00  |

Images and data from Borg2010.

![img-13.jpeg](img-13.jpeg)

---

DGP - Applications

# Nanostructures

- Unlabeled distances

Juhás et.al, Ab initio determin. of solid-state nanostruct., Nature, 2006.
Billinge et. al, Assigned and unassigned distance geometry, 4OR, 2016.

![img-14.jpeg](img-14.jpeg)

Image from Billinge2016.

---

DGP - Applications

# NMR

- Protein structure determination

|  Method | #Proteins  |
| --- | --- |
|  X-Rray | 105,656  |
|  NMR | 10,257  |
|  Electron Microscopy | 993  |
|  Hybrid | 97  |
|  Other | 181  |
|  Total | 117,184  |

(From: PDB - Accessed in 01/30/2017).

![img-15.jpeg](img-15.jpeg)

---

DGP - Applications

## NMR

- Protein structure determination

NMR doesn’t give positions, but only info about the distance between some pairs of atoms

$$
l_{ij} \leq \left\| x_i - x_j \right\| \leq u_{ij},
$$

$$
(i, j) \in E \subseteq \{1, 2, \dots, n\}^2
$$

![img-16.jpeg](img-16.jpeg)

---

# Algorithms and the underlying hypothesis

---

DGP - Algorithms

Global Optimization – General Methods (No additional hypothesis)

![img-17.jpeg](img-17.jpeg)

- The number of local minima grows exponentially with the problem size;
- BISwas and Ye, SDP for ad hoc wireless sensor network localization, 2004.
- Souza et. al. Solving the MDGP with inaccurate distance data. BMC, 2013.

---

DGP - Algorithms

Geometric build-up – Additional hypothesis

$$
H_1: (i, j) \in E, \forall i, j;
$$

$$
H_2: l_{ij} = u_{ij}, \forall (i, j) \in E;
$$

![img-18.jpeg](img-18.jpeg)

![img-19.jpeg](img-19.jpeg)

DGP can be solved by SVD in $O(n^3)$ operations.

© Crippen and Havel, DG and Molecular Conformation, 1988.

---

DGP - Algorithms

# Geometric build-up - Additional hypothesis

$H_{1}\colon \forall i &gt; 4,$  exist at least four atoms  $i_1,i_2,i_3,i_4 &lt;   i,$  such that

$(i,i_k)\in E,k = 1,2,3,4;$

$H_{2}\colon l_{ij} = u_{ij},\forall (i,j)\in E;$

![img-20.jpeg](img-20.jpeg)

![img-21.jpeg](img-21.jpeg)

DGP can be solved in  $O(n)$  operations.

Dong and Wu, A linear-time algorithm for solving the MDGP with exact inter-atomic distances, 2002.

---

DGP - Algorithms

Geometric build-up – Additional hypothesis

$H_{1}\colon \forall i &gt; 3,$ exist at least three atoms $i_1,i_2,i_3 &lt; i$, such that

$$
(i, i _ {k}) \in E, k = 1, 2, 3;
$$

$$
H _ {2} \colon l _ {i j} = u _ {i j}, \forall (i, j) \in E;
$$

![img-22.jpeg](img-22.jpeg)

![img-23.jpeg](img-23.jpeg)

Find all solutions of DGP.

In practice, $O(n)$ operations, but $O(n^{2})$ in worst case.

Liberti et. al, A branch-and-prune algorithm for the MDGP. ITORS, 2008.

---

DGP - Algorithms

Geometric build-up – Additional hypothesis

$H_{1}\colon \forall i &gt; 3,\text{exist at least three}$

atoms $i_1, i_2, i_3 &lt; i$, such that

$$
(i, i _ {k}) \in E, k = 1, 2, 3;
$$

$$
H _ {2} \colon l _ {i j} = u _ {i j}, \forall (i, j) \in E;
$$

Are these hypothesis still too strong?

For protein backbones, what are the guaranteed valid hypothesis?

---

DGP - Algorithms

## Geometric build-up – Additional hypothesis

## Assumptions

- Distances between close hydrogen (&lt; 5A) atoms will be estimated by NMR;
- Length of covalent bonds are known;
- Planar angles are known;

![img-24.jpeg](img-24.jpeg)

Lavor et. al, On the computation of protein backbones by using artificial backbones of hydrogens. Journal of Global Optimization, 2011.

---

DGP - Algorithms

Geometric build-up – Additional hypothesis

$$
H _ {1} \colon (i, j) \in E, \forall | i - j | &lt;   4;
$$

$$
H _ {2} \colon l _ {i j} = u _ {i j}, \forall | i - j | &lt;   3;
$$

![img-25.jpeg](img-25.jpeg)

Lavor et. al, On the computation of protein backbones by using artificial backbones of hydrogens. Journal of Global Optimization, 2011.

---

DGP - Algorithms

# Geometric build-up - Additional hypothesis

$$
H _ {1} \colon (i, j) \in E, \forall | i - j | &lt;   4;
$$

$$
H _ {2} \colon l _ {i j} = u _ {i j}, \forall | i - j | &lt;   3;
$$

The only degree of freedom is given by the dihedral (torsion) angles;

$$
x _ {k} = R \left(\omega_ {k}, x _ {k} ^ {0}, x _ {k - 1}, x _ {k - 2}\right)
$$

$$
x _ {k} \colon [ 0, 2 \pi ] ^ {k} \to \mathbb {R} ^ {3}
$$

![img-26.jpeg](img-26.jpeg)

This parametrization reduces the number of variables from 3n to n.

The search space is reduced to a hypercube.

---

3

## Torsion Angles

![img-27.jpeg](img-27.jpeg)

![img-28.jpeg](img-28.jpeg)

![img-29.jpeg](img-29.jpeg)

---

6

## Torsion Angles

![img-30.jpeg](img-30.jpeg)

![img-31.jpeg](img-31.jpeg)

![img-32.jpeg](img-32.jpeg)

---

2

## Torsion Angles

![img-33.jpeg](img-33.jpeg)

![img-34.jpeg](img-34.jpeg)

![img-35.jpeg](img-35.jpeg)

---

2

# Torsion Angles

![img-36.jpeg](img-36.jpeg)

![img-37.jpeg](img-37.jpeg)

![img-38.jpeg](img-38.jpeg)

---

3

## Torsion Angles

![img-39.jpeg](img-39.jpeg)

![img-40.jpeg](img-40.jpeg)

![img-41.jpeg](img-41.jpeg)

---

2

## Torsion Angles

![img-42.jpeg](img-42.jpeg)

![img-43.jpeg](img-43.jpeg)

![img-44.jpeg](img-44.jpeg)

---

69

## Torsion Angles

![img-45.jpeg](img-45.jpeg)

![img-46.jpeg](img-46.jpeg)

![img-47.jpeg](img-47.jpeg)

---

DGP - Algorithms

# Geometric build-up - Additional hypothesis

$$
R (\theta , p, u, v) = (1 - \cos (\theta)) \left[ \begin{array}{l} u _ {1} (w _ {2} ^ {2} + w _ {3} ^ {2}) - w _ {1} (\kappa - u _ {1} w _ {1}) \\ u _ {2} (w _ {1} ^ {2} + w _ {3} ^ {2}) - w _ {2} (\kappa - u _ {2} w _ {2}) \\ u _ {3} (w _ {1} ^ {2} + w _ {2} ^ {2}) - w _ {3} (\kappa - u _ {3} w _ {3}) \end{array} \right]
$$

$$
\begin{array}{l} \text {+ \sin (\theta)} \left[ \begin{array}{l} - u _ {3} w _ {2} + u _ {2} w _ {3} - w _ {3} p _ {2} + w _ {2} p _ {3} \\ u _ {3} w _ {1} - u _ {1} w _ {3} + w _ {3} p _ {1} - w _ {1} p _ {3} \\ - u _ {2} w _ {1} + u _ {1} w _ {2} - w _ {2} p _ {2} + w _ {1} p _ {2} \end{array} \right] + \cos (\theta) p, \end{array}
$$

where  $\omega = \frac{u - v}{||u - v||}$  is the normalized direction and  $\kappa = \omega \cdot (u - v)$ .

---

DGP - Algorithms

# Geometric build-up - Additional hypothesis

$$
R (\theta , p, u, v) = (1 - \cos (\theta)) \left[ \begin{array}{l} u _ {1} (w _ {2} ^ {2} + w _ {3} ^ {2}) - w _ {1} (\kappa - u _ {1} w _ {1}) \\ u _ {2} (w _ {1} ^ {2} + w _ {3} ^ {2}) - w _ {2} (\kappa - u _ {2} w _ {2}) \\ u _ {3} (w _ {1} ^ {2} + w _ {2} ^ {2}) - w _ {3} (\kappa - u _ {3} w _ {3}) \end{array} \right]
$$

$$
\begin{array}{l} \left. \begin{array}{c} v \\ \theta \\ p \end{array} \right\} + \sin (\theta) \left[ \begin{array}{c} - u _ {3} w _ {2} + u _ {2} w _ {3} - w _ {3} p _ {2} + w _ {2} p _ {3} \\ u _ {3} w _ {1} - u _ {1} w _ {3} + w _ {3} p _ {1} - w _ {1} p _ {3} \\ - u _ {2} w _ {1} + u _ {1} w _ {2} - w _ {2} p _ {2} + w _ {1} p _ {2} \end{array} \right] + \cos (\theta) p, \\ \sigma_ {j k} ^ {i} = \frac {\partial x _ {j} ^ {i} (\omega)}{\partial \omega_ {k}} = \nabla^ {t} R (\omega_ {i}, x _ {j} ^ {i - 1}, x _ {i + 1} ^ {i - 1}, x _ {i + 2} ^ {i - 1}) (\delta_ {i k}, \sigma_ {j k} ^ {i - 1}, \sigma_ {i + 1, k} ^ {i - 1}, \sigma_ {i + 2, k} ^ {i - 1}) \\ \end{array}
$$

---

DGPω

$[DGP\omega]$ Find $\omega \in [0,2\pi]^n$ such that

$$
l_{ij} \leq \left\| x_i(\omega) - x_j(\omega) \right\| \leq u_{ij}, \quad (i,j) \in E \subseteq \{1,2,\dots,n\}^2,
$$

$$
x_k = R\left(\omega_k, x_k^0, x_{k-1}, x_{k-2}\right),
$$

$$
\omega_k \in [0,2\pi].
$$

Fewer variables and bounded search space.

---

# Proof of concept

---

DGPω

# Numerical Experiments

- Prototype in Matlab;
- 100 random instances with only 8 atoms (toy model);

|   | ns | nf | time(sec)  |
| --- | --- | --- | --- |
|  QN | 67 | 328.18 | 0.35  |
|  SPHω | 90 | 101.56 | 0.98  |

- More solutions (ns), with fewer func evals in average (nf);
- However, the time was bigger.

---

# Conclusions

We presented a parametrization of search space of DGP that

1. Is based on (more) reasonable hypothesis;
2. Reduces the number of variables in 2/3;
3. Reduces the search space to a hypercube;
4. Derivatives are explicitly calculated.

## Future work:

5. Develop an efficient implementation;
6. Explore how other algorithms can benefit by using this model.

Lack of NMR especialists and real instances.

---

Thank you.
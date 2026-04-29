# Proximal Cyclic BCD for Graph Realization
## Set-Up
Fix the anchor

$$x_1=a_1=(0,0),$$

and optimize only over $x_2,\ldots,x_n$. Let

$$z=(x_2,\ldots,x_n)\in \mathbb R^{2(n-1)}.$$

Define the anchored objective

$$F(z)=\frac12\sum_{\{i,j\}\in E}\left(\|x_i-x_j\|^2-d_{ij}^2\right)^2.$$

A point $z^*$ is called stationary if

$$\nabla_{x_i}F(z^*)=0,\qquad i=2,\ldots,n.$$

Equivalently, $z^*$ satisfies the first-order condition for the anchored problem.

## Proximal Cyclic Block Coordinate Descent (BCD)

Let $p=n-1$, and cyclically update $i=2,\ldots,n$. Given $z^r$, define intermediate points $z^{r,0}=z^r$, then for $\kappa=1,\ldots,p$, with $i_\kappa=\kappa+1$,

$$x_{i_\kappa}^{r,\kappa} \in \arg\min_{y\in\mathbb R^2} 
\left\{
F(z^{r,\kappa-1}_{i_\kappa},y) + \frac{\mu}{2}\|y-x_{i_\kappa}^{r,\kappa-1}\|^2
\right\},$$

where $F(z^{r,\kappa-1}_{i_\kappa},y)$ denotes the objective with block $i_\kappa$ set to $y$ and all other free blocks fixed at their values in $z^{r,\kappa-1}$, and where $\mu>0$. All other blocks remain fixed. Set

$$z^{r+1}=z^{r,p}.$$

The subproblem has a solution because the objective is continuous and coercive in $y$.

## Convergence Analysis
### Lemma 1: Bounded Level Sets
For every $\alpha$, the anchored sublevel set

$$\{z:F(z)\le \alpha\}$$

is bounded.

**Proof.**
If $F(z)\le \alpha$, then for every edge $\{i,j\}$,

$$\frac12\left(\|x_i-x_j\|^2-d_{ij}^2\right)^2\le \alpha,$$

hence

$$\|x_i-x_j\|^2\le d_{ij}^2+\sqrt{2\alpha}.$$

Thus every edge length is uniformly bounded on the sublevel set. Since $G$ is connected, every vertex is connected to the fixed anchor $x_1$ by a path of bounded edge lengths. Therefore every $x_i$ is bounded. ∎

### Lemma 2: Sufficient Descent
For every sweep,

$$F(z^r)-F(z^{r+1})
\ge
\frac{\mu}{2}\|z^{r+1}-z^r\|^2.$$

Consequently, $F(z^r)$ converges and

$$\sum_{r=0}^\infty \|z^{r+1}-z^r\|^2<\infty,
\qquad
\|z^{r+1}-z^r\|\to0.$$

**Proof.**
At block $\kappa$, compare the minimizer with the old value $x_{i_\kappa}^{r,\kappa-1}$. This gives

$$F(z^{r,\kappa})
+
\frac{\mu}{2}
\|x_{i_\kappa}^{r,\kappa}-x_{i_\kappa}^{r,\kappa-1}\|^2
\le
F(z^{r,\kappa-1}).$$

Summing over $\kappa=1,\ldots,p$,

$$F(z^r)-F(z^{r+1})
\ge
\frac{\mu}{2}
\sum_{\kappa=1}^p
\|x_{i_\kappa}^{r,\kappa}-x_{i_\kappa}^{r,\kappa-1}\|^2=
\frac{\mu}{2}\|z^{r+1}-z^r\|^2$$

since each block is updated once per sweep.

Since $F\ge0$, the sequence $\{F(z^r)\}$ is monotonically non-increasing and bounded below, hence convergent to some $\bar F\ge0$. Summing the descent inequality from $r=0$ to $m$ gives

$$\frac{\mu}{2}\sum_{r=0}^m \|z^{r+1}-z^r\|^2
\le
F(z^0)-F(z^{m+1})
\le
F(z^0).$$

Letting $m\to\infty$ yields

$$\sum_{r=0}^\infty \|z^{r+1}-z^r\|^2<\infty,$$

and therefore $\|z^{r+1}-z^r\|\to0$. ∎

### Lemma 3: Accumulation Points Are Stationary
Every accumulation point of $\{z^r\}$ is stationary.

**Proof.**
By Lemma 1 and Lemma 2, $\{z^r\}$ is bounded and $\|z^{r+1}-z^r\|\to0$. Take a convergent subsequence $z^{r_j}\to z^*$.

The first-order condition for block $\kappa$ is

$$\nabla_{x_{i_\kappa}}F(z^{r,\kappa})
+
\mu
\left(x_{i_\kappa}^{r,\kappa}-x_{i_\kappa}^{r,\kappa-1}\right)
=0.$$

The block increments tend to zero, and $z^{r_j,\kappa}\to z^*$ for every fixed $\kappa$. Passing to the limit gives

$$\nabla_{x_{i_\kappa}}F(z^*)=0.$$

Since this holds for every free block $i_\kappa=2,\ldots,n$, $z^*$ is stationary. ∎

### Lemma 4: Relative Error Bound
There exists $C>0$ such that

$$\|\nabla F(z^{r+1})\|
\le
C\|z^{r+1}-z^r\|.$$

**Proof.**
All iterates lie in the bounded sublevel set $\{F\le F(z^0)\}$. Since $F$ is polynomial, $\nabla F$ is Lipschitz on the convex hull of this set; say with constant $L$.

For block $i_\kappa$,

$$\nabla_{x_{i_\kappa}}F(z^{r,\kappa})
=
-\mu
\left(x_{i_\kappa}^{r,\kappa}-x_{i_\kappa}^{r,\kappa-1}\right).$$

Therefore,

$$\begin{aligned}
\|\nabla_{x_{i_\kappa}}F(z^{r+1})\|
&\le
\|\nabla_{x_{i_\kappa}}F(z^{r,\kappa})\|
+
\|\nabla_{x_{i_\kappa}}F(z^{r+1})
-\nabla_{x_{i_\kappa}}F(z^{r,\kappa})\| \\
&\le
\mu\|x_{i_\kappa}^{r,\kappa}-x_{i_\kappa}^{r,\kappa-1}\|
+
L\|z^{r+1}-z^{r,\kappa}\| \\
&\le
\mu\|z^{r+1}-z^r\|
+
L\|z^{r+1}-z^r\|.
\end{aligned}$$

Combining the finitely many block estimates gives the claim. ∎

### KL Property
Since $F$ is a polynomial, it satisfies the Kurdyka-Łojasiewicz property at every point. That is, near every stationary point $z^*$, there exist $\eta>0$, a neighborhood $U$, and a function $\varphi:[0,\eta)\to\mathbb R_+$ that is continuous on $[0,\eta)$, concave and $C^1$ on $(0,\eta)$, with $\varphi(0)=0$ and $\varphi'>0$ on $(0,\eta)$, such that

$$\varphi'(F(z)-F(z^*))\|\nabla F(z)\|\ge1$$

whenever

$$z\in U,\qquad F(z^*)<F(z)<F(z^*)+\eta.$$

This is stated without proof (see [[1]](#ref1), [[2]](#ref2)).

### Theorem: Whole-Sequence Convergence
The proximal cyclic BCD sequence $\{z^r\}$ converges to a stationary point $z^*$. Moreover,

$$\sum_{r=0}^\infty \|z^{r+1}-z^r\|<\infty.$$

Thus the flattened block-update sequence also converges to the same $z^*$.

**Proof.**
By Lemma 2, $F(z^r)\downarrow \bar F$. By boundedness, choose an accumulation point $z^*$. By continuity,

$$F(z^*)=\bar F.$$

By Lemma 3, $z^*$ is stationary, so the KL property applies at $z^*$.

Let

$$s_r=F(z^r)-\bar F.$$

If $s_r=0$ for some $r$, then Lemma 2 forces $z^{r+1}=z^r$, and the sequence is constant afterward. Otherwise $s_r>0$.

Let

$$\Delta_r=\|z^{r+1}-z^r\|.$$

By the KL property at $z^*$, there are $\rho>0$, $\eta>0$, and a desingularizing function $\varphi$ such that the KL inequality holds for every $z\in B(z^*,\rho)$ satisfying

$$F(z^*)<F(z)<F(z^*)+\eta.$$

Let $a=\mu/2$ and let $b>0$ be the constant from Lemma 4, so that

$$s_r-s_{r+1}\ge a\Delta_r^2, \qquad \|\nabla F(z^r)\|\le b\Delta_{r-1},\qquad r\ge1.$$

Since $z^*$ is an accumulation point, $s_r\to0$, $\Delta_r\to0$, and $\varphi(s_r)\to\varphi(0)=0$, we may choose $N\ge1$ along a subsequence converging to $z^*$ such that

$$s_N<\eta,\qquad
\|z^N-z^*\|+\Delta_{N-1}+\frac{b}{a}\varphi(s_N)<\rho.$$

We claim that $z^r\in B(z^*,\rho)$ for every $r\ge N$. Suppose first that $z^N,\ldots,z^M\in B(z^*,\rho)$ and that $s_r>0$ for $r=N,\ldots,M$. Then the KL inequality and the relative error bound give

$$\varphi'(s_r)
\ge
\frac{1}{b\Delta_{r-1}}.$$

Since $\varphi$ is concave,

$$\varphi(s_r)-\varphi(s_{r+1})
\ge
\varphi'(s_r)(s_r-s_{r+1})
\ge
\frac{a\Delta_r^2}
{b\Delta_{r-1}}.$$

Hence

$$\Delta_r^2
\le
\frac{b}{a}
\Delta_{r-1}
\left(\varphi(s_r)-\varphi(s_{r+1})\right).$$

Taking square roots and applying $2\sqrt{uv}\le u+v$ (AM-GM) with $u=\Delta_{r-1}$ and $v=\tfrac{b}{a}(\varphi(s_r)-\varphi(s_{r+1}))$ gives

$$2\Delta_r
\le
\Delta_{r-1}
+
\frac{b}{a}
\left(\varphi(s_r)-\varphi(s_{r+1})\right).$$

Summing from $r=N$ to $M$ and rearranging,

$$2\sum_{r=N}^M\Delta_r
\le
\sum_{r=N}^M\Delta_{r-1}
+
\frac{b}{a}\sum_{r=N}^M\bigl(\varphi(s_r)-\varphi(s_{r+1})\bigr).$$

The right-hand sum of $\Delta_{r-1}$ equals $\sum_{r=N-1}^{M-1}\Delta_r$, so subtracting $\sum_{r=N}^{M-1}\Delta_r$ from both sides leaves

$$\sum_{r=N}^M\Delta_r+\Delta_M\le\Delta_{N-1}+\frac{b}{a}\bigl(\varphi(s_N)-\varphi(s_{M+1})\bigr).$$

More usefully, the $\varphi$-sum telescopes to $\varphi(s_N)-\varphi(s_{M+1})\le\varphi(s_N)$ gives

$$\sum_{r=N}^M \Delta_r + \Delta_M
\le
\Delta_{N-1}
+
\frac{b}{a}\varphi(s_N).$$

This bounds the total path length from $z^N$ onward by quantities fixed at step $N$. Applying the triangle inequality,

$$\begin{aligned}
\|z^{M+1} - z^*\|
&\le \|z^N - z^*\| + \|z^{N+1} - z^N\| + \cdots + \|z^{M+1} - z^M\| \\
&= \|z^N - z^*\| + \Delta_N + \cdots + \Delta_M \\
&= \|z^N - z^*\| + \sum_{r=N}^M \Delta_r
\end{aligned}$$

So
$$\|z^{M+1}-z^*\|
\le
\|z^N-z^*\|+\sum_{r=N}^M \Delta_r
\le
\|z^N-z^*\|+\Delta_{N-1}+\frac{b}{a}\varphi(s_N)
<
\rho,$$

where the last inequality uses exactly the choice of $N$ made above. Since $M$ is arbitrary, this proves the claim by induction. If $s_r=0$ for some $r\ge N$, then the sequence is constant afterward, as above; otherwise the preceding estimates hold for all $r\ge N$. Letting $M\to\infty$ gives

$$\sum_{r=N}^\infty \Delta_r<\infty.$$

Thus $\{z^r\}$ is Cauchy and converges to some $\bar z$. Since the subsequence used to choose $N$ still has infinitely many terms after $N$ and converges to $z^*$, we must have $\bar z=z^*$. Finally, Lemma 4 gives

$$\|\nabla F(z^r)\|\to0,$$

so the limit is stationary. ∎

## References

1. <a id="ref1"></a> Attouch, H., Bolte, J., and Svaiter, B. F. (2013). Convergence of descent methods for semi-algebraic and tame problems: proximal algorithms, forward–backward splitting, and regularized Gauss–Seidel methods. *Mathematical Programming*, 137(1), 91–129.

2. <a id="ref2"></a> Bolte, J., Sabach, S., and Teboulle, M. (2014). Proximal alternating linearized minimization for nonconvex and nonsmooth problems. *Mathematical Programming*, 146(1), 459–494.

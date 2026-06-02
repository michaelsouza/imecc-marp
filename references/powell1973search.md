Mathematical Programming 4 (1973) 193–201. North-Holland Publishing Company

# ON SEARCH DIRECTIONS FOR MINIMIZATION ALGORITHMS

M.J.D. POWELL

Atomic Energy Research Establishment, Harwell, Great Britain

Received 11 August 1972

Revised manuscript received 4 January 1973

Some examples are given of differentiable functions of three variables, having the property that if they are treated by the minimization algorithm that searches along the coordinate directions in sequence, then the search path tends to a closed loop. On this loop the gradient of the objective function is bounded away from zero. We discuss the relevance of these examples to the problem of proving general convergence theorems for minimization algorithms that use search directions.

# 1. Introduction

Suppose that we wish to calculate the least value of a differentiable function $F(x_{1},x_{2},\ldots ,x_{n}) = F(\pmb {x})$, say, of $n$ real variables, given a starting point $\pmb{x}^{(1)}$. Further suppose that the set $\{x\colon F(x)\leq F(x^{(1)})\}$ is bounded, and that on this set the gradient $\pmb {g} = \pmb {\nabla}F(\pmb {x})$ is continuous. Then, if the steepest descent algorithm is used to minimize $F(x)$, it will generate a sequence of points $\pmb{x}^{(k)}$ ($k = 1,2,3,\dots$) having the property that the corresponding gradient vectors $\pmb{g}^{(k)} = \pmb {\nabla}F(\pmb{x}^{(k)})$ converge to zero [1].

Because it does not seem to be over-ambitious to find other algorithms for unconstrained minimization that also have this sureness of convergence, we consider in this paper the minimization algorithm that seeks the least value of $F(x)$ by adjusting the variables in sequence. Specifically, for $k = 1, 2, \dots, n$, the point $x^{(k + 1)}$ is calculated from $x^{(k)}$ by changing the $k^{\text{th}}$ component of $x^{(k)}$ so that

$$
F \left(\boldsymbol {x} ^ {(k + 1)}\right) = \min  _ {\boldsymbol {x}} F \left(x _ {1} ^ {(k)}, \dots , x _ {k - 1} ^ {(k)}, x, x _ {k + 1} ^ {(k)}, \dots , x _ {n} ^ {(k)}\right).
$$

When $\pmb{x}^{(n + 1)}$ has been found, another search is made along the first coordinate direction, and so the process continues iteratively changing the components of $x^{(k)}$ in sequence. We call this algorithm “the method that changes one variable at a time”.

I used to believe that this method had convergence properties that are like those of the steepest descent algorithm, because if one is at a point where $\nabla F(x)$ is nonzero, then the value of $F(x)$ can be reduced appreciably by searching along the coordinate direction that is nearest to the steepest descent direction. However, the main purpose of this paper is to present numerical examples showing that the method that changes one variable at a time may cycle without calculating any point where the gradient of $F(x)$ is small. There is no conflict with Polak’s [2] convergence theorem on the method of local variations, because, although the search directions of this method are the coordinate directions, they are not used in a strict cyclic order.

These numerical examples are given in Section 2. The first one is straightforward, and it includes functions of the form

$$
(x - c) _ {+} ^ {2} = \left\{ \begin{array}{l l} 0, &amp; \text{if } x \leq c, \\ (x - c) ^ {2}, &amp; \text{if } x \geq c. \end{array} \right. \tag {1}
$$

However, the remarkable properties of this example can be destroyed by making a small perturbation to the starting vector $x^{(1)}$. Therefore we extend the example so that the limiting behaviour of the iterative method that changes one variable at a time is not sensitive to either small changes in the initial data or to small errors introduced during the iterative process, for example computer rounding errors. Thirdly, there is an example of a function that is infinitely differentiable that also causes an endless loop in the iterative minimization method that changes one variable at a time.

It is pointed out in Section 3 that the examples emphasize the difficulty of proving convergence theorems for algorithms that minimize general functions. In particular they show that iterative algorithms that use linearly independent search directions may fail to find any stationary points of smooth functions.

## 2. Three examples

For all three examples, the objective function depends on only three variables. Therefore, to avoid subscripts, we write $(x,y,z)$ in place of $(x_{1},x_{2},x_{3})$.

First consider applying the method that changes one variable at a time to the function

$$
\varphi_{1}(x, y, z) = -xy - yz - zx + (x - 1)_{+}^{2} + (-x - 1)_{+}^{2} + (y - 1)_{+}^{2} + (-y - 1)_{+}^{2} + (z - 1)_{+}^{2} + (-z - 1)_{+}^{2}. \tag{2}
$$

In this case, a search along the $x$-direction from the point $(x, y, z)$ changes the $x$-coordinate to the value $\{1 + \frac{1}{2}|y + z|\}$ sign $(y + z)$, a search along the $y$-direction from the point $(x, y, z)$ changes the $y$-coordinate to the value $\{1 + \frac{1}{2}|z + x|\}$ sign $(z + x)$, and a search along the $z$-direction from the point $(x, y, z)$ changes the $z$-coordinate to the value $\{1 + \frac{1}{2}|x + y|\}$ sign $(x + y)$. It follows that if the starting value of $(x, y, z)$ is the point $(-1 - \epsilon, 1 + \frac{1}{2}\epsilon, -1 - \frac{1}{4}\epsilon)$, then the first six steps of the algorithm generate the six points $(1 + \frac{1}{6}\epsilon, 1 + \frac{1}{2}\epsilon, -1 - \frac{1}{4}\epsilon)$, $(1 + \frac{1}{6}\epsilon, -1 - \frac{1}{4}\epsilon, -1 - \frac{1}{4}\epsilon)$, $(1 + \frac{1}{8}\epsilon, -1 - \frac{1}{4}\epsilon, 1 + \frac{1}{32}\epsilon)$, $(-1 - \frac{1}{64}\epsilon, -1 - \frac{1}{4}\epsilon, 1 + \frac{1}{32}\epsilon)$, $(-1 - \frac{1}{64}\epsilon, 1 + \frac{1}{128}\epsilon, 1 + \frac{1}{32}\epsilon)$ and $(-1 - \frac{1}{64}\epsilon, 1 + \frac{1}{128}\epsilon, -1 - \frac{1}{256}\epsilon)$. Now this last point is the same as the starting point except that $\epsilon$ has been replaced by $\frac{1}{64}\epsilon$. Therefore the calculated sequence of points does not converge to a single limit point. Instead, the steps of the method tend to cycle round six edges of the cube whose vertices are $(\pm 1, \pm 1, \pm 1)$.

The most remarkable property of this example is that on the limiting path the gradient vector of the objective function is bounded away from zero. Specifically on this path the gradient vector is $(-y - z, -z - x, -x - y) = (g_x, g_y, g_z)$, say, and it satisfies the equation

$$
\left| g_{x} \right| + \left| g_{y} \right| + \left| g_{z} \right| = 2. \tag{3}
$$

In fact this example is unstable with respect to small perturbations. In other words, small changes in the starting vector $(-1 - \epsilon, 1 + \frac{1}{2}\epsilon, -1 - \frac{1}{4}\epsilon)$, or small errors in the numbers that are computed during the calculation, almost always destroy the cyclic behaviour. We verify this statement by letting the initial point be $(-1 - \epsilon_1, 1 + \epsilon_2, -1 - \epsilon_3)$ and by letting the sequence of points that are calculated by the iterative process be $(1 + \epsilon_4, 1 + \epsilon_2, -1 - \epsilon_3), (1 + \epsilon_4, -1 - \epsilon_5, -1 - \epsilon_3), (1 + \epsilon_4, -1 - \epsilon_5, 1 + \epsilon_6), (-1 - \epsilon_7, -1 - \epsilon_5, 1 + \epsilon_6), (-1 - \epsilon_7, 1 + \epsilon_8, 1 + \epsilon_6), (-1 - \epsilon_7, 1 + \epsilon_8, -1 - \epsilon_9)$, etc. Then we find that the successive values of $\epsilon_k$ satisfy the recurrence relation

$$
\epsilon_{k} = \frac{1}{2} \left(\epsilon_{k-2} - \epsilon_{k-1}\right), \tag{4}
$$

provided that $\epsilon_{k-2} &gt; \epsilon_{k-1}$. This condition is necessary to ensure that the step that yields the value of $\epsilon_k$ does correspond to an edge of the cube whose vertices are $(\pm 1, \pm 1, \pm 1)$. Otherwise the cyclic behaviour of the example finishes.

Now the solution of the recurrence relation (4) has the form

$$
\epsilon_k = A \left(\frac{1}{2}\right)^k + B (-1)^k, \tag{5}
$$

where $A$ and $B$ are constants, and where so far we have considered only the case when $B = 0$. However, practically any perturbation to the data causes an effect that is equivalent to $B$ being nonzero, and in this case for sufficiently large $k$ the second term of expression (5) must dominate the first term. It follows that the condition $\epsilon_{k-2} \leq \epsilon_{k-1}$ holds eventually, so we have shown that the cyclic behaviour of our example is unstable with respect to small perturbations.

Therefore, in order to show that the cyclic behaviour may occur in practical computations, we now extend the function (2) to a form that gives stable cycling. Specifically, we consider applying the method that changes one variable at a time to the function

$$
\begin{aligned}
\varphi_2(x, y, z) = &amp; -xy - yz - zx + 3\psi(y)\psi(-z)(x - \frac{2}{3} + \frac{1}{3}z)_+^2 \\
&amp; + 3\psi(-y)\psi(z)(-x - \frac{2}{3} - \frac{1}{3}z)_+^2 \\
&amp; + 3\psi(z)\psi(-x)(y - \frac{2}{3} + \frac{1}{3}x)_+^2 \\
&amp; + 3\psi(-z)\psi(x)(-y - \frac{2}{3} - \frac{1}{3}x)_+^2 \\
&amp; + 3\psi(x)\psi(-y)(z - \frac{2}{3} + \frac{1}{3}y)_+^2 \\
&amp; + 3\psi(-x)\psi(y)(-z - \frac{2}{3} - \frac{1}{3}y)_+^2, \tag{6}
\end{aligned}
$$

where $\psi(\theta)$ is a monotonically increasing differentiable function that satisfies the conditions

$$
\begin{aligned}
\psi(\theta) &amp;\equiv 0, \quad \theta \leq -\frac{1}{2}, \\
\psi(\theta) &amp;\equiv 1, \quad \theta \geq \frac{1}{2}. \tag{7}
\end{aligned}
$$

If a search is made in the $x$-direction from the point $(x, y, z)$ to minimize the value of $\varphi_2(x, y, z)$, and if $y &gt; -z \geq 1$, then the required value of $x$ minimizes the function

$$
\begin{aligned}
-x(y+z) + 3\left(x - \frac{2}{3} + \frac{1}{3}z\right)_+^2 + 3\psi(x)\left(-y - \frac{2}{3} - \frac{1}{3}x\right)_+^2 \\
+ 3\psi(-x)\left(-z - \frac{2}{3} - \frac{2}{3}y\right)_+^2. \tag{8}
\end{aligned}
$$

Now the next to last term of this expression is identically zero when $y &gt; 1$, so it follows that the required value of $x$ is $\frac{1}{6}(4 + y - z)$. Similarly, if $-y &gt; z \geq 1$, then a search from $(x, y, z)$ along the $x$-direction changes the $x$-coordinate to the value $\frac{1}{6}(-4 + y - z)$. Further, by symmetry, if a search is made in the $y$-direction from the point $(x, y, z)$, then the $y$-coordinate is changed to $\frac{1}{6}(4 + z - x)$ if $z &gt; -x \geq 1$, and it is changed to $\frac{1}{6}(-4 + z - x)$ if $-z &gt; x \geq 1$. Also, if a search is made in the $z$-direction from the point $(x, y, z)$ then the $z$-coordinate is changed to $\frac{1}{6}(4 + x - y)$ if $x &gt; -y \geq 1$, and it is changed to $\frac{1}{6}(-4 + x - y)$ if $-x &gt; y \geq 1$. Therefore, if again the starting point for the method that changes one variable at a time is $(-1 - \epsilon, 1 + \frac{1}{2}\epsilon, -1 - \frac{1}{4}\epsilon)$, then again the first six steps of the method provide the sequence of points $(1 + \frac{1}{6}\epsilon, 1 + \frac{1}{2}\epsilon, -1 - \frac{1}{4}\epsilon)$, $(1 + \frac{1}{6}\epsilon, -1 - \frac{1}{16}\epsilon, -1 - \frac{1}{4}\epsilon)$, $(1 + \frac{1}{6}\epsilon, -1 - \frac{1}{16}\epsilon, 1 + \frac{1}{32}\epsilon)$, $(-1 - \frac{1}{64}\epsilon, -1 - \frac{1}{16}\epsilon, 1 + \frac{1}{32}\epsilon)$, $(-1 - \frac{1}{64}\epsilon, 1 + \frac{1}{128}\epsilon, 1 + \frac{1}{32}\epsilon)$ and $(-1 - \frac{1}{64}\epsilon, 1 + \frac{1}{128}\epsilon, -1 - \frac{1}{256}\epsilon)$. Thus again we have the cyclic behaviour that was shown in our first example, and again equation (3) holds on the limiting path.

To consider the stability properties of our second example, we define the numbers $\epsilon_{k}$ ($k = 1,2,3,\ldots$) in the same way as before. Therefore, for instance, the first step of the method that changes one variable at a time is from the point $(-1 - \epsilon_{1}, 1 + \epsilon_{2}, -1 - \epsilon_{3})$ to $(1 + \epsilon_{4}, 1 + \epsilon_{2}, -1 - \epsilon_{3})$. It follows from the previous paragraph that if $\epsilon_{2} &gt; \epsilon_{3} \geq 0$, then $(1 + \epsilon_{4})$ has the value $\frac{1}{6}(6 + \epsilon_{2} + \epsilon_{3})$. In other words, the recurrence relation

$$
\epsilon_{k} = \frac{1}{6} \left(\epsilon_{k-2} + \epsilon_{k-1}\right) \tag{9}
$$

is satisfied for $k = 4$. Now it is straightforward to show that this recurrence relation holds for all subsequent values of $k$ provided that the numbers in the sequence $(\epsilon_{2}, \epsilon_{3}, \epsilon_{4}, \ldots)$ decrease strictly monotonically, and the solution of the recurrence relation has the form

$$
\epsilon_{k} = A \left(\frac{1}{2}\right)^{k} + B \left(-\frac{1}{3}\right)^{k}. \tag{10}
$$

We note in particular that, in contrast to equation (5), the second term of the solution decreases in magnitude faster than the first term. Therefore the cyclic behaviour of the second example is stable with respect to small perturbations, such as computer rounding errors. Indeed the numbers $\epsilon_{k}$ ($k = 2, 3, 4, \ldots$) decrease strictly monotonically provided only that $\epsilon_{2}$ and $\epsilon_{3}$ satisfy the conditions

$$
\epsilon_{2} &gt; \epsilon_{3} &gt; \frac{1}{5} \epsilon_{2} &gt; 0. \tag{11}
$$


Therefore the cyclic behaviour is attained from any starting vector in a wedge-shaped region of the space of the variables.

Because the functions (2) and (6) have discontinuous second derivatives, it is interesting to ask whether the cyclic behaviour can occur for smoother objective functions. Our third example answers this question, for it provides cycling although it is infinitely differentiable, and again the gradient of the objective function is bounded away from zero on the limiting path.

To define this example, we let $\mu(\theta)$ be the infinitely differentiable function

$$
\mu(\theta) = \begin{cases}
0 &amp; \text{if } \theta \leq 0, \\
\frac{1}{2} \theta \exp(1 - \theta^{-1}) &amp; \text{if } \theta &gt; 0,
\end{cases} \tag{12}
$$

and we impose on the monotonically increasing function $\psi(\theta)$, that satisfies condition (7), the extra condition that it too is to be infinitely differentiable. Here I wish to acknowledge that the existence of such a function $\psi(\theta)$ was shown to me by A.R. Curtis, for he suggested constructing it by taking a suitable multiple of the indefinite integral of the function

$$
\lambda(\theta) = \begin{cases}
0 &amp; \text{if } \theta \leq -\frac{1}{2}, \\
\exp\{(4\theta^2 - 1)^{-1}\} &amp; \text{if } -\frac{1}{2} &lt; \theta &lt; \frac{1}{2}, \\
0 &amp; \text{if } \theta \geq \frac{1}{2}.
\end{cases} \tag{13}
$$

In this example, the objective function is the expression

$$
\begin{aligned}
\varphi_3(x, y, z) &amp;= -xy - yz - zx - (1 + z) \psi(y) \psi(-z) \mu(x + \frac{1}{2} (1 + z)) \\
&amp;\quad - (1 - z) \psi(-y) \psi(z) \mu(-x + \frac{1}{2} (1 - z)) \\
&amp;\quad - (1 + x) \psi(z) \psi(-x) \mu(y + \frac{1}{2} (1 + x)) \\
&amp;\quad - (1 - x) \psi(-z) \psi(x) \mu(-y + \frac{1}{2} (1 - x)) \\
&amp;\quad - (1 + y) \psi(x) \psi(-y) \mu(z + \frac{1}{2} (1 + y)) \\
&amp;\quad - (1 - y) \psi(-x) \psi(y) \mu(-z + \frac{1}{2} (1 - y)).
\end{aligned} \tag{14}
$$

Therefore, if a search is made in the $x$-direction from the point $(x, y, z)$ and if $y &gt; -z \geq 1$, then instead of calculating $x$ to minimize the function (8), we calculate it to minimize the expression

$$
\begin{aligned}
-x(y + z) &amp;- (1 + z) \mu(x + \frac{1}{2} (1 + z)) - (1 - x) \psi(x) \mu(-y + \frac{1}{2} (1 - x)) \\
&amp;\quad - (1 - y) \psi(-x) \mu(-z + \frac{1}{2} (1 - y)).
\end{aligned} \tag{15}
$$

Note that again the next to last term is identically zero, and that the last term is also unimportant if the required value of $x$ is greater than $\frac{1}{2}$. We let $y = 1 + \frac{1}{2}\epsilon$ and $z = -1 - \frac{1}{4}\epsilon$, and in this case the required value of $x$ minimizes the expression

$$
\frac {1}{4} \epsilon \left\{\mu \left(x - \frac {1}{8} \epsilon\right) - x \right\} = \sigma (x), \tag {16}
$$

say. Now from the definition (12) we obtain the derivatives

$$
\sigma^ {\prime} (x) = \left\{ \begin{array}{l l} - \frac {1}{4} \epsilon &amp; \text {if} x \leq \frac {1}{8} \epsilon , \\ \frac {1}{4} \epsilon \left\{\frac {1}{2} \left(1 + \theta^ {- 1}\right) \exp \left(1 - \theta^ {- 1}\right) - 1 \right\} &amp; \text {if} x &gt; \frac {1}{8} \epsilon , \end{array} \right. \tag {17}
$$

$$
\sigma^ {\prime \prime} (x) = \left\{ \begin{array}{l l} 0 &amp; \text {if} x \leq \frac {1}{8} \epsilon , \\ \frac {1}{8} \epsilon \theta^ {- 3} \exp (1 - \theta^ {- 1}) &amp; \text {if} x &gt; \frac {1}{8} \epsilon , \end{array} \right. \tag {18}
$$

where $\theta = x - \frac{1}{8}\epsilon$. It follows that $\sigma(x)$ is a convex function whose least value is obtained when $\theta = 1$. Thus, given the starting point $(-1 - \epsilon, 1 + \frac{1}{2}\epsilon, -1 - \frac{1}{4}\epsilon)$, the first step of the method that changes one variable at a time reaches the point $(1 + \frac{1}{8}\epsilon, 1 + \frac{1}{2}\epsilon, -1 - \frac{1}{4}\epsilon)$. By continuing this analysis it may be shown that the sequence of points that is calculated for the objective function (14) is the same as the sequence of points for the objective function (2). In particular, the limiting behaviour is again that the cycle $(-1, 1, -1) \to (1, 1, -1) \to (1, -1, -1) \to (1, -1, 1) \to (-1, -1, 1) \to (-1, 1, 1) \to (-1, 1, -1)$ is repeated endlessly.

Equation (3) is not satisfied by the function (14) on the limiting path. For example on the line segment that joins $(-1, 1, -1)$ to $(1, 1, -1)$ we have the derivative values

$$
g _ {x} = 0, \quad g _ {y} = (1 - x) + \frac {1}{2} \psi (- x), \quad g _ {z} = - (1 + x) - \mu (x), \tag {19}
$$

and on the line segment that joins $(1, -1, 1)$ to $(-1, -1, 1)$ we have the derivative values

$$
g _ {x} = 0, \quad g _ {y} = - (1 + x) - \frac {1}{2} \psi (x), \quad g _ {z} = (1 - x) + \mu (- x). \tag {20}
$$

Thus, because of the symmetry of the function (14), the inequality

$$
\left| g _ {x} \right| + \left| g _ {y} \right| + \left| g _ {z} \right| \geq 2 \tag {21}
$$

is satisfied on the limiting path.

Moreover, this example is not stable with respect to perturbations. To support this statement, we note that if $y$ is close to one and $z$ is close to minus one, then the value of $x$ that minimizes expression (15) is close to one only if the ratio $(y + z)/(1 + z)$ is close to minus one. The instability occurs because of the cancellation in the numerator and denominator of this ratio.

I do not know if it is possible to find objective functions with continuous second derivatives, having the property that, if they are treated by the method that changes one variable at a time, then cycling occurs, the gradients are bounded away from zero on the limiting path, and the cycling is stable with respect to small perturbations.

## 3. Discussion

Although the examples show that the method of changing one variable at a time may fail to find the least value of the objective function, this fact alone is rather unimportant. What is important is that there are many very useful algorithms that use search directions, and we have shown that one of these algorithms may fail. Therefore, to prove that a different algorithm cannot fail in the way shown by our examples, one has to use properties of the algorithm that are not shared by the method that changes one variable at a time.

In particular, because the method that changes one variable at a time uses the coordinate directions as search directions, one cannot prove that an algorithm will not cycle just by showing that it uses search directions with good linear independence properties. Instead, one probably has to consider the angles between search directions and gradients. For example, the first theorem given by Zoutendijk [4] shows that if an iterative method calculates the point $x^{(k+1)}$ from $x^{(k)}$ by searching in the direction $d^{(k)}$ ($k = 1, 2, 3, \ldots$), and if we let $c^{(k)}$ be the cosine of the angle between $d^{(k)}$ and $\nabla F(x^{(k)})$, then the cyclic behaviour of the examples of Section 2 cannot occur if the sum $\Sigma(c^{(k)})^2$ is divergent.

Finally we recall some remarks made by Rosenbrock [3]. He suggests that instead of using the method that changes one variable at a time, it is better to rotate the set of orthogonal search directions after each direction has been used once. He makes this suggestion in order to reduce the number of line searches that are required to achieve a given accuracy, but now we can add to his argument that the rotation of search directions may be necessary to prevent the cyclic behaviour shown by the examples of Section 2.

# Acknowledgments

This work was begun by some discussions with R.W.H. Sargent at the stimulating NATO meeting on Mathematical Programming, held at Figueira da Foz in June 1972. Also I am grateful to A.R. Curtis for his comments on the manuscript, and in particular for his help with identifying infinitely differentiable functions for the third example of Section 2.

# References

[1] H.B. Curry, “The method of steepest descent for nonlinear minimization problems”, Quarterly of Applied Mathematics 2 (1944) 258–261.
[2] E. Polak, Computational methods in optimization: a unified approach (Academic Press, New York, 1971).
[3] H.H. Rosenbrock, "An automatic method for finding the greatest or least value of a function", Computer Journal 3 (1960) 175-184.
[4] G. Zoutendijk, “Nonlinear programming, computational methods”, in: Integer and nonlinear programming, Ed. J. Abadie (North-Holland, Amsterdam, 1970) 37–86.
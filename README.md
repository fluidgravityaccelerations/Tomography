# Tomography

## Preconditioning for CT Reconstruction (`preconditioning_CAT.m`)

Role of **preconditioning** when solving linear systems arising in **X-ray computed tomography (CT)**. The algorithm reconstructs a Shepp–Logan phantom from simulated Radon projections, forms the normal equations  

$$(A^\top A)\,x = A^\top b$$,

and solves it using **Conjugate Gradient** (CG) with two different preconditioning strategies:
  - **Jacobi-preconditioned CG** 
  - **Incomplete Cholesky (IC) preconditioned CG** - it uses auto-adaptive strategies to stabilize `ichol()` for ill-conditioned CT matrices  

---

### Jacobi Preconditioning 
A simple diagonal scaling based on the diagonal of \(A^\top A\):

$$M_{\text{Jac}}^{-1} = \mathrm{diag}(A^\top A)^{-1}$$.

Usually provides a mild accelaration, but it is very cheap.

---

### Incomplete Cholesky (IC) Preconditioning  
Builds a sparse approximation of the Cholesky factorization:

$$A^\top A + \tau I \approx L L^\top$$.

Typically leads to much faster convergence.

The script approximate the Cholesky factor but drop entries during factorization to preserve sparsity and:

- tests multiple regularization values $\tau$  
- switches between  
  - `type = 'nofill'` - When one computes a full Cholesky factorization, new nonzeros appear all over the matrix (“fill-in”). Even if $A^\top A$ is sparse, the Cholesky factor $L$ becomes dense. With `type = 'nofill'`, only entries that are already nonzero in the matrix $A^\top A$ are allowed to appear in the Cholesky factor $L$. Fast (almost no overhead), memory-efficient, produces a very sparse $𝐿$, works well when the matrix $A^\top A$ is diagonally dominant or nearly banded. However, it is a weaker preconditioner, often fails on ill-conditioned matrices (like those in tomography), it is more likely to produce nonpositive pivots, causing 1) error using `ichol`; 2) encountered nonpositive pivot.
  - `type = 'ict'` with `droptol` - threshold-based dropping, stronger - This is why in tomography one often needs a diagonal shift $\tau I$ or a different type ('ict'). `type = 'ict'` with `droptol` allows some fill-in, but drops small elements depending on droptol. It is a much stronger preconditioner.
  - falls back to Jacobi if `ichol()` fails.

Requires tuning of `tau`, `type`, and `droptol`.

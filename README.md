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

The script:

- tests multiple regularization values $\tau$  
- switches between  
  - `type = 'nofill'` (very sparse, weaker)  
  - `type = 'ict'` with `droptol` (threshold-based dropping, stronger)  
- falls back to Jacobi if `ichol()` fails

**Effects:**

- Stronger reduction of condition number  
- Often yields noticeably fewer iterations  
- Requires tuning of `tau`, `type`, and `droptol`

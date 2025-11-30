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


# Algebraic Reconstruction Technique (ART) – Theoretical Notes

This repository implements a tomographic reconstruction algorithm based on the **Algebraic Reconstruction Technique (ART)**.  
This README provides a concise but rigorous explanation of:

- the **standard ART update formula**,  
- the role of \( A^T A \),  
- why \( A^T A \) is **never explicitly formed in practice**,  
- how **matrix-free implementations** using `radon` and `iradon` handle the update implicitly,  
- how the **relaxation parameter** \( \lambda \) compensates for missing row normalization.

---

## 1. The Standard ART Update Formula

Let

- \( A \in \mathbb{R}^{M \times N} \) be the system matrix,
- \( a_i^T \) be the \( i \)-th row of \( A \),
- \( x \in \mathbb{R}^N \) be the image to reconstruct,
- \( b_i \in \mathbb{R} \) be the measured projection for ray \( i \).

The classical ART / Kaczmarz row-relaxation update is:

\[
x^{k+1}
= x^k +
\lambda_k \,
\frac{\,b_i - a_i^T x^k\,}{\|a_i\|^2} \, a_i.
\]

This comes from projecting the current iterate onto the solution hyperplane of equation:

\[
a_i^T x = b_i.
\]

### Meaning of the terms

- **Forward projection**: \( a_i^T x^k = (Ax^k)_i \)  
- **Residual**: \( r_i = b_i - a_i^T x^k \)  
- **Backprojection of a single ray**: \( a_i \)  
- **Normalization**: division by \( \|a_i\|^2 = \sum_j a_{ij}^2 \)  
- **Relaxation parameter** \( \lambda_k \): controls convergence speed and stability.

---

## 2. Why We *Do Not* Form \( A^T A \)

A naïve interpretation of ART suggests constructing the normal equations:

\[
A^T A x = A^T b.
\]

This is **never done** in tomography because:

1. **A is enormous** (often \(10^6 - 10^8\) nonzeros).  
2. \(A^TA\) is **dense**, even if \(A\) is sparse.  
3. Explicitly computing or storing \(A^TA\) is **computationally infeasible**.  
4. Forming \(A^TA\) destroys the ray geometry and introduces smoothing that is **not physically meaningful**.

Thus, modern tomography always uses **matrix-free operators**:

- **\(Ax\)** → forward projection (Radon transform)  
- **\(A^T y\)** → backprojection (adjoint Radon transform)

without ever assembling \(A\) or \(A^T A\).

---

## 3. How `radon` and `iradon` Replace Explicit Matrix Multiplication

MATLAB's:

- `radon(x, theta)` computes **\( Ax \)**,
- `iradon(p, theta, 'linear', 'none')` computes **\( A^T p \)**.

Therefore, the ART update becomes:

\[
x^{k+1}
= x^k +
\lambda_k \, A^T (b - Ax).
\]

### Where did the division by \( \|a_i\|^2 \) go?

In matrix-free ART, the exact row normalization is **absorbed into** the relaxation parameter:

\[
\lambda_k \quad \longrightarrow \quad
\frac{\lambda_k}{\|a_i\|^2}.
\]

Why?

- We do not have explicit access to each \( \|a_i\|^2 \).  
- Forward and backward projectors are designed so their **relative magnitudes** approximate the correct scaling.  
- A global tuning parameter \( \lambda \) stabilizes the iteration.

This is standard practice in:

- SART (Simultaneous ART),  
- OS-ART (Ordered Subsets ART),  
- RAMLA (Row Action Maximum Likelihood Algorithm),  
- ASTRA, TIGRE, ODL, and MATLAB iterative tomography.

---

## 4. Practical ART Update Used in This Repository

The implemented update is:

\[
x^{k+1}
= x^k + \lambda_k \, A^T(b - Ax).
\]

This is a *matrix-free approximation* of the classical ART rule.  
Key points:

### ✔ Multiplication by \( a_i \)
Handled by the **backprojection operator** \( A^T \).

### ✔ Multiplication by \( a_i^T \)
Handled by the **forward projector** \( A \).

### ✔ Division by \( \|a_i\|^2 \)
Absorbed into the choice of \( \lambda_k \), which plays the role of a **combined relaxation and normalization factor**.

---

## 5. Choice and Decay of the Relaxation Parameter \( \lambda_k \)

A typical schedule is:

\[
\lambda_k = \frac{\lambda_0}{1 + \alpha k},
\]

or exponential decay:

\[
\lambda_k = \lambda_0 \, \rho^k,
\quad 0 < \rho < 1.
\]

The decay serves to:

- prevent oscillations caused by inconsistent data,
- compensate for approximate normalization,
- improve convergence stability,
- handle noise in the sinogram.

---

## 6. Summary

| Concept | Classical ART | Matrix-Free ART (radon / iradon) |
|--------|----------------|----------------------------------|
| \( Ax \) | explicit multiplication | `radon` |
| \( A^T y \) | explicit multiplication | `iradon` |
| multiply by \( a_i \) | direct | backprojection |
| normalize by \( \|a_i\|^2 \) | explicit | absorbed into \( \lambda \) |
| assemble \( A^TA \) | required | **never used** |
| memory use | huge | minimal |

This makes iterative tomography feasible for large-scale 2D and 3D problems.

---

## 7. References

- Kak & Slaney, *Principles of Computerized Tomographic Imaging*  
- Herman, *Fundamentals of Computerized Tomography*  
- Andersen & Kak, *SART: A new reconstruction algorithm*  

---

If you use this README, please cite the relevant references and adapt the equations to the exact version of ART implemented in your code.



# Tomography


# MATLAB Demo: Preconditioning for CT Reconstruction

This repository contains a MATLAB demo illustrating the role of **preconditioning** when solving linear systems arising in **X-ray computed tomography (CT)**. The example reconstructs a Shepp–Logan phantom from simulated Radon projections, forms the normal equations  
\[
(A^\top A)\,x = A^\top b,
\]  
and solves them using **Conjugate Gradient** (CG) with three different preconditioning strategies.

The goal is to show students how preconditioning affects convergence speed, stability, and reconstruction quality.

---

## Features

- Builds an explicit system matrix \(A\) using MATLAB’s `radon` operator  
- Solves the CT inverse problem using:  
  - **Unpreconditioned CG**  
  - **Jacobi-preconditioned CG**  
  - **Incomplete Cholesky (IC) preconditioned CG**  
- Uses auto-adaptive strategies to stabilize `ichol()` for ill-conditioned CT matrices  
- Compares reconstructions visually  
- Plots **residual histories** to highlight convergence improvements

---

## Preconditioning Strategies

### 1. Jacobi Preconditioning  
A simple diagonal scaling based on the diagonal of \(A^\top A\):

\[
M_{\text{Jac}}^{-1} = \mathrm{diag}(A^\top A)^{-1}.
\]

**Effects:**

- Normalizes the magnitude of the unknowns  
- Usually provides mild acceleration  
- Very cheap and robust

---

### 2. Incomplete Cholesky (IC) Preconditioning  
Builds a sparse approximation of the Cholesky factorization:

\[
A^\top A + \tau I \approx L L^\top.
\]

This preconditioner significantly improves eigenvalue clustering, which typically leads to much faster CG convergence.

The script automatically:

- tests multiple regularization values \(\tau\)  
- switches between  
  - `type = 'nofill'` (very sparse, weaker)  
  - `type = 'ict'` with `droptol` (threshold-based dropping, stronger)  
- falls back to Jacobi if `ichol()` fails

**Effects:**

- Stronger reduction of condition number  
- Often yields noticeably fewer iterations  
- Requires tuning of `tau`, `type`, and `droptol`

---

## Output

The script displays:

- The true phantom  
- Reconstructions from:
  - **CG**
  - **Jacobi PCG**
  - **IC-PCG**  
- A semilog plot of **relative residual history**, clearly showing how preconditioning changes convergence behavior.

---


# Computing the leading RV eigenpair without storing RV

The CLI now defaults to `--rv-method matrix-free`. This changes how the same
RV eigenproblem is solved; it does not change the DISTATIS definition, axis
selection, outlier thresholds, imputation or filtering loop.

For centered gene matrices S_i, define

```text
<A,B> = sum_r A_rr B_rr + 2 sum_{r<c} A_rc B_rc
d_i   = sqrt(<S_i,S_i>)
RV_ij = <S_i,S_j> / (d_i d_j)
```

Let B have one row per gene, containing its diagonal and sqrt(2)-weighted upper
triangle, divided by d_i. Then RV = B B^T, so RV x = B (B^T x).
Neither B nor RV needs to be stored as an additional matrix.

The implementation uses the equivalent operations

```text
T = sum_j (x_j / d_j) S_j            # only diagonal and upper triangle needed
y_i = <S_i,T> / d_i                 # y = RV x
```

This preserves the R code's upper-triangle convention even for asymmetric
matrices produced by its legacy nearest-neighbor imputation fallback.

## Solver and accuracy

A fully reorthogonalized Krylov subspace with Rayleigh–Ritz extraction finds
the leading eigenpair. It starts with both a constant and a deterministic
pseudorandom direction; a constant-only power iteration can miss a leading
eigenvector orthogonal to that start. The subspace is capped at 48 vectors,
with explicit restarts retaining leading Ritz vectors.

Convergence is checked against a fresh application of the original operator:

```text
norm(RV v - lambda v) / max(1, abs(lambda)) <= 2e-14
```

Failure to converge within 2,048 products raises an explicit error suggesting
`--rv-method dense`. There is no silent fallback that allocates a large RV
matrix. As with other iterative eigensolvers, the residual measures convergence
of the computed pair; it is not a universal proof that an arbitrary starting
subspace captures every eigenspace.

With tied leading eigenvalues, weights are not uniquely defined. Nearly tied
eigenvalues and values on outlier/acceptance boundaries can make results
sensitive to floating-point details, including in the R implementation.
The dense reference method remains available for comparisons.

## Cost and retained data

For K genes and N taxa:

* Each implicit product costs O(K N^2).
* Extra solver memory is O(K m + N^2), with m <= 48.
* Dense RV construction costs O(K^2 N^2), followed by a full O(K^3)
  eigendecomposition in the previous C++ implementation.

The centered/input matrices still require O(K N^2) storage. The compromise
eigendecomposition and automatic factor selection are unchanged. This change
does not yet implement packed distance matrices or a reduced-axis fast mode.

Independent aggregate rows and output vector entries can run in parallel.
Each entry retains the same summation order across thread counts.

## Usage and diagnostics

```sh
phylter run --trees genes.nwk --out results
phylter run --trees genes.nwk --out dense-reference --rv-method dense
```

The normal CLI path does not allocate RV. `--diagnostics DIRECTORY` explicitly
requests RV output and therefore materializes it **after accepting a state**,
while still obtaining the eigenpair with the selected solver. Each diagnostic
state also includes a `.rv-solver.tsv` file with operator-product count and
relative residual (zero products and an unevaluated `NA` residual for the dense solver).

For C++ API compatibility, direct `distatis()` calls retain RV by default.
Pass `RvOptions{RvMethod::matrix_free, false}` to suppress that output.
`analyze()` uses the compact setting by default; observers that need RV must
explicitly set `options.rv.materialize = true`.

Validation covers all 18 R-reference configurations, comparisons with the
previous dense C++ output, explicit operator-versus-matrix products, negative
entries, rank deficiency, restarts, near/tied leading eigenvalues, an
eigenvector orthogonal to the constant start, non-convergence, and a 5,000-gene
regression that does not allocate RV. See the current benchmark report for
measured speed and RAM rather than extrapolating the asymptotic costs.

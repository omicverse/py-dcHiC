# MATH.md — py-dcHiC derivations

## 1. Acceleration: full SVD → symmetric eigendecomposition (exact identity)

**Rewrite.** Replace `np.linalg.svd(C2, full_matrices=False)` with `np.linalg.eigh(C2)`,
keeping the eigenvectors of the `n_pcs` largest eigenvalues.

**Claim.** For the matrix `C2` that dcHiC feeds to its decomposition, this is an *exact*
identity on the gated output (compartment PC1).

**Proof.** `C2` is a Pearson correlation matrix of the (rounded) correlation matrix `C1`:

```
Z2 = column-zscore(C1)         (each column mean 0, sample-std 1)
C2 = Z2ᵀ Z2 / (n-1)
```

`C2` is symmetric (`C2ᵀ = C2`) and positive semi-definite (`xᵀC2x = ‖Z2 x‖² /(n-1) ≥ 0`).
For a symmetric PSD matrix the singular value decomposition `C2 = U Σ Vᵀ` and the
eigendecomposition `C2 = Q Λ Qᵀ` coincide: `U = V = Q` and `Σ = Λ`. Hence the right
singular vectors returned by `svd` and the eigenvectors returned by `eigh` span the same
per-eigenvalue invariant subspaces. The PC scores are `PC = Z2 · V`.

Two residual degrees of freedom:

1. **Per-vector sign.** `eigh` and `svd` may return `±qₖ`. dcHiC's orientation step
   `PC[:,k] ← sign(corr(PC[:,k], gcc)) · PC[:,k]` (dchicf.r:526) removes the sign, so the
   final oriented PC is invariant to the choice.
2. **Degenerate eigenvalues.** Within an exactly repeated eigenvalue the eigenbasis is only
   defined up to rotation, so `svd` and `eigh` may disagree on vectors *inside* such a block.
   PC1 corresponds to the **dominant, well-separated** eigenvalue (the compartment axis),
   which is non-degenerate, so PC1 is uniquely determined and identical under both routines.

**Empirical confirmation.** On the canonical fixture, compartment PC1 from the `eigh` path
matches the dcHiC C++/`svd` reference to `max |Δ| = 1.6e-14` (ordinal `|r| = 1.0`), i.e. bit
equivalence to f64 rounding. Speedup: 1.83× end-to-end Stage-A on the fixture (n=200);
3.0–3.2× on the isolated decomposition at n=1000–2000 (`eigh` exploits symmetry via
tridiagonal reduction and skips the discarded singular vectors). See `ITERATION_LOG.md`.

No bounded-ε approximation is introduced, so there is no `Σ bound` term to track and the
manifest stays at the `ordinal` threshold (0.99) for this output.

## 2. Documented cross-implementation deviations (NOT acceleration rewrites)

These are unavoidable library substitutions, not algebraic rewrites. They are why the Stage-B
gate is the rank-based **inference** class (Spearman on −log10 p + top-K Jaccard), per the
omicverse-rebuildr FAQ on `mgcv≠pygam`-style divergence — never widen the gate to force
element-wise agreement.

### 2.1 `robust::covRob` → `sklearn.covariance.MinCovDet`

dcHiC estimates the condition×condition covariance for the Mahalanobis quadratic form with
`robust::covRob`, whose default selects the **Donoho–Stahel** estimator for `n<1000, p<10`.
No faithful Python port of Donoho–Stahel exists; `MinCovDet` (FastMCD) is the closest robust
estimator. The condition-level statistic is

```
sample_maha[k] = max_min_diff[k] · (g[k] − c[k])ᵀ Σ⁻¹ (g[k] − c[k])
```

where `max_min_diff[k] = (maxₙ gₖₙ − minₙ gₖₙ)²` is the squared cross-condition range. This
scalar effect-size factor is **identical** under both estimators and dominates the ranking, so
differing `Σ` perturbs the *magnitude* but largely preserves the *order* of bins by
significance. Measured: Spearman on −log10 p = 0.98, top-50 Jaccard = 0.96, and all 89 of R's
`padj<0.1` bins are recovered.

### 2.2 `IHW::ihw` → Benjamini–Hochberg

For multiple-testing correction with replicates dcHiC uses Independent Hypothesis Weighting
(`IHW`, Bioconductor), with `replicate_wt` as the covariate. IHW learns data-dependent group
weights via convex optimization; it has no faithful Python implementation. The port uses
Benjamini–Hochberg (`p.adjust(p,'BH')` is also dcHiC's no-replicate path). Both are monotone
in the per-bin p-value, so the −log10 p ranking that the inference gate measures is preserved.
The raw chi-square p-value (pre-adjustment) is bit-faithful to dcHiC modulo §2.1.

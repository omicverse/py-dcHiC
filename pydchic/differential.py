"""Stage B — differential compartment analysis (dcHiC ``analyze``).

Pure-Python re-implementation of dcHiC's differential statistic. Mirrors
``pcanalyze`` / ``calcen`` / ``mdweight`` in ``dchicf.r``.

Per genomic bin, across conditions (each with replicates):

1. Per-condition quantile-normalise replicate PC vectors (``limma::normalizeQuantiles``,
   ties=TRUE), then scale each column by its mean absolute value -> ``pcQnm`` matrix.
2. Replicate-level Mahalanobis distance from a shrunk per-replicate centre using a
   diagonal covariance -> ``replicate_wt = -log P(chi2_{nrep-1} <= d)`` (covariate for FDR).
3. Condition means -> ``intra_grp`` (bins x conditions).
4. Condition-level Mahalanobis distance from a shrunk centre (``calcen``), with the
   inverse robust covariance (``robust::covRob``; here ``sklearn.MinCovDet``) scaled by the
   squared cross-condition range ``max_min_diff``; p-value = ``chi2_{ncond-1}`` survival.
5. Optional refinement: re-estimate the robust covariance on the null bins (p > 1-rconf)
   and recompute.
6. Multiple-testing adjustment. dcHiC uses ``IHW`` (covariate = ``replicate_wt``) when
   replicates exist; for cross-language determinism this port uses Benjamini-Hochberg
   (documented deviation — IHW has no faithful Python equivalent; see MATH.md / README).

Default parameters match the dcHiC CLI: ``rzscore=2``, ``szscore=0``, ``refine=True``,
``rconf=0.90``.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import chi2, norm, rankdata

__all__ = [
    "normalize_quantiles",
    "quantile_normalize_by_condition",
    "calcen",
    "bh_adjust",
    "differential_compartments",
]


def normalize_quantiles(a: np.ndarray) -> np.ndarray:
    """``limma::normalizeQuantiles(A, ties=TRUE)`` (bins x samples).

    Sort each column, average the sorted values across columns at each rank, then map
    each observation back through its average rank by linear interpolation.
    """
    a = np.asarray(a, dtype=np.float64)
    n = a.shape[0]
    s = np.sort(a, axis=0)
    m = s.mean(axis=1)                      # mean of sorted values per rank
    grid = np.arange(1, n + 1, dtype=np.float64)
    out = np.empty_like(a)
    for j in range(a.shape[1]):
        r = rankdata(a[:, j])               # average ranks (ties -> mean rank)
        out[:, j] = np.interp(r, grid, m)
    return out


def quantile_normalize_by_condition(pc_raw: np.ndarray, conditions) -> np.ndarray:
    """Per-condition quantile normalisation + mean-abs scaling (dchicf.r:880-882)."""
    conditions = np.asarray(conditions)
    out = np.empty_like(pc_raw, dtype=np.float64)
    seen = []
    for c in conditions:
        if c not in seen:
            seen.append(c)
    for c in seen:
        cols = np.where(conditions == c)[0]
        qn = normalize_quantiles(pc_raw[:, cols])
        # dcHiC scales by column mean-abs via `qn / apply(abs(qn),2,mean)` (dchicf.r:881).
        # R divides a matrix by a length-p vector with column-MAJOR recycling, not per
        # column. Replicate that exactly (it matters only when p does not divide n_bins;
        # after quantile normalisation the column means are ~equal so the effect is tiny,
        # but we match R's literal semantics to be safe).
        colmeans = np.abs(qn).mean(axis=0)
        flat = qn.flatten(order="F")
        div = np.resize(colmeans, flat.size)
        out[:, cols] = (flat / div).reshape(qn.shape, order="F")
    return out


def _rowwise_spread(df: np.ndarray) -> np.ndarray:
    """``sqrt(colSums(as.matrix(dist(x)))/(p-1))`` per row (dchicf.r:753)."""
    p = df.shape[1]
    absdiff = np.abs(df[:, :, None] - df[:, None, :])   # n x p x p
    colsums = absdiff.sum(axis=2)                        # n x p
    return np.sqrt(colsums / (p - 1))


def calcen(df: np.ndarray, cls: str, rzscore: float = 2.0, szscore: float = 0.0) -> np.ndarray:
    """Shrunk centre estimate (dchicf.r:751).

    ``cls='rep'`` z-scores the spread over the whole matrix; ``cls='sample'`` z-scores
    per column. ``cen = round(df * (1 - rowmax(pnorm(z, mean=zscore))), 5)``.
    """
    df = np.asarray(df, dtype=np.float64)
    df_dist = _rowwise_spread(df)
    if cls == "rep":
        val = df_dist.ravel()
        z = (df_dist - val.mean()) / val.std(ddof=1)
        pvl = np.round(norm.cdf(z - rzscore), 5)
    elif cls == "sample":
        z = np.empty_like(df_dist)
        for n in range(df_dist.shape[1]):
            col = df_dist[:, n]
            z[:, n] = (col - col.mean()) / col.std(ddof=1)
        pvl = np.round(norm.cdf(z - szscore), 5)
    else:
        raise ValueError(f"cls must be 'rep' or 'sample', got {cls!r}")
    pvl_max = pvl.max(axis=1)
    return np.round(df * (1.0 - pvl_max[:, None]), 5)


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (``p.adjust(p, 'BH')``)."""
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1.0)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=np.float64)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def _robust_cov(x: np.ndarray, seed: int = 123) -> np.ndarray:
    """Robust covariance, mirroring ``robust::covRob`` (here sklearn FastMCD)."""
    from sklearn.covariance import MinCovDet
    return MinCovDet(random_state=seed).fit(x).covariance_


def differential_compartments(
    pc_qnm: np.ndarray,
    conditions,
    rzscore: float = 2.0,
    szscore: float = 0.0,
    refine: bool = True,
    rconf: float = 0.90,
    seed: int = 123,
) -> dict:
    """Differential compartment statistic. Columns of ``pc_qnm`` must be grouped by
    condition (matching ``conditions``). Returns ``sample_maha``, ``pval``, ``padj`` and
    (when replicates exist) ``replicate_wt``.
    """
    pc_qnm = np.asarray(pc_qnm, dtype=np.float64)
    conditions = np.asarray(conditions)
    nbins, ntot = pc_qnm.shape

    order_conditions = []
    for c in conditions:
        if c not in order_conditions:
            order_conditions.append(c)

    cen_cols, grp_cols = [], []
    rep_count = 0
    for c in order_conditions:
        cols = np.where(conditions == c)[0]
        sub = pc_qnm[:, cols]
        if sub.shape[1] > 1:
            cen_cols.append(calcen(sub, "rep", rzscore, szscore))
            grp_cols.append(sub.mean(axis=1))
            rep_count += 1
        else:
            cen_cols.append(np.full((nbins, 1), sub.mean()))
            grp_cols.append(sub[:, 0])

    intra_cen = np.hstack(cen_cols)                     # bins x ntot
    intra_grp = np.column_stack(grp_cols)               # bins x ncond
    ncond = intra_grp.shape[1]

    out = {}
    if rep_count > 0:
        cov_all = np.cov(pc_qnm, rowvar=False, ddof=1)
        inv_var = 1.0 / np.diag(cov_all)                # diagonal covariance
        diff = pc_qnm - intra_cen
        replicate_maha = (diff ** 2 * inv_var).sum(axis=1)
        out["replicate_wt"] = -chi2.logcdf(replicate_maha, df=ntot - 1)

    intra_grp_cen = calcen(intra_grp, "sample", rzscore, szscore)
    max_min_diff = (intra_grp.max(axis=1) - intra_grp.min(axis=1)) ** 2
    diffg = intra_grp - intra_grp_cen

    def _maha_pval(cov_rows: np.ndarray):
        cov_rob = _robust_cov(cov_rows, seed=seed)
        inv_cov = np.linalg.inv(cov_rob)
        quad = np.einsum("ij,jk,ik->i", diffg, inv_cov, diffg)
        smaha = max_min_diff * quad
        return smaha, chi2.sf(smaha, df=ncond - 1)

    sample_maha, pval = _maha_pval(intra_grp)
    if refine:
        thr = 1.0 - rconf                               # = chi2.sf(chi2.ppf(rconf,df),df)
        null_mask = pval > thr
        sample_maha, pval = _maha_pval(intra_grp[null_mask])

    finite = np.isfinite(pval)
    padj = np.full(nbins, np.nan)
    padj[finite] = bh_adjust(pval[finite])

    out["sample_maha"] = sample_maha
    out["pval"] = pval
    out["padj"] = padj
    return out

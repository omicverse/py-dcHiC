"""Stage A — compartment calling (dcHiC ``cis`` + ``select``).

Pure-Python re-implementation of dcHiC's A/B-compartment kernel. Faithfully mirrors
``functionsdchic`` (C++ ``oe2cor`` / ``zmat`` / ``ijk2mat``) and the R glue in
``ijk2matfunc_cis`` / ``readfilesintra`` / ``pcselectioncore`` (``dchicf.r``).

Pipeline (per intra-chromosomal contact matrix):

1. Observed/Expected: expected count per genomic distance d is
   ``sum(weight at d) / n_pairs(d)`` where ``n_pairs(d) = max(n_bins - d/res, 0)``
   (``expectedInteraction``, dchicf.r:204); ``WeightOE = weight / expected[d]``.
2. Build dense symmetric O/E matrix M (``ijk2mat``); blacklist bins with row-sum < 3.
3. C1 = column Pearson-correlation matrix of M  (``oe2cor`` = column z-score with
   ddof=1, then ``Zᵀ Z / (n-1)``).
4. ``mat2fbm``: set diag(C1)=1, round to 5 decimals.
5. C2 = column Pearson-correlation matrix of (rounded) C1; let Z2 = z-score(C1).
6. set diag(C2)=1, round to 5 decimals; V = top-k right singular vectors of C2 (SVD).
7. **PC = Z2 · V**  (dchicf.r:300).
8. Orient each PC by GC content: ``PC <- sign(cor(PC, gcc)) · PC`` (dchicf.r:526).

The A/B compartment score is the (GC-oriented) PC1 column.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "expected_interaction",
    "observed_over_expected",
    "ijk2mat",
    "oe_zscore",
    "correlation_matrix",
    "compartment_pca",
    "orient_pc",
    "call_compartments",
]

ROUND_DECIMALS = 5  # mat2fbm rounds the correlation matrix to 5 dp (dchicf.r:28-33)


def expected_interaction(dist_bp: np.ndarray, n_bins: int, resolution: int) -> np.ndarray:
    """Number of valid bin-pairs at each genomic distance (dchicf.r:204).

    ``vecInt = n_bins - dist/res``; the count is ``sum(vecInt[vecInt > 0])``. For a
    single chromosome this reduces to ``max(n_bins - dist/res, 0)`` per distance.
    """
    steps = np.asarray(dist_bp, dtype=np.float64) / resolution
    vec = n_bins - steps
    return np.where(vec > 0, vec, 0.0)


def observed_over_expected(
    a_idx: np.ndarray,
    b_idx: np.ndarray,
    weight: np.ndarray,
    pos: np.ndarray,
    n_bins: int,
    resolution: int,
    count_thr: float = 0.0,
    minexpcc: float = 0.0,
):
    """Convert raw intra-chromosomal contacts to O/E ratios (dchicf.r:357-375).

    Parameters mirror dcHiC: ``a_idx``/``b_idx`` are 1-based bin indices, ``weight`` the
    raw contact counts, ``pos`` the bin start coordinates (bp) indexed by 1-based bin id.
    Returns ``(a_idx, b_idx, weight_oe)`` for entries surviving ``weight > count_thr``.
    """
    a_idx = np.asarray(a_idx)
    b_idx = np.asarray(b_idx)
    weight = np.asarray(weight, dtype=np.float64)

    keep = weight > count_thr
    a_idx, b_idx, weight = a_idx[keep], b_idx[keep], weight[keep]

    dist = np.abs(pos[a_idx - 1] - pos[b_idx - 1]).astype(np.float64)

    # sum of weights per distance (aggregate Weight ~ dist, sum)
    uniq_dist, inv = np.unique(dist, return_inverse=True)
    sum_w = np.zeros(uniq_dist.shape[0], dtype=np.float64)
    np.add.at(sum_w, inv, weight)

    total_pairs = expected_interaction(uniq_dist, n_bins, resolution)
    with np.errstate(divide="ignore", invalid="ignore"):
        expcc = np.where(total_pairs > 0, sum_w / total_pairs, 0.0)

    # floor expcc (dchicf.r:371): if min(expcc) > minexpcc keep min else minexpcc
    floor = expcc.min() if expcc.min() > minexpcc else minexpcc
    expcc = np.where(expcc <= floor, floor, expcc)

    expcc_by_dist = dict(zip(uniq_dist.tolist(), expcc.tolist()))
    denom = np.array([expcc_by_dist[d] for d in dist], dtype=np.float64)
    weight_oe = weight / denom
    return a_idx, b_idx, weight_oe


def ijk2mat(a_idx: np.ndarray, b_idx: np.ndarray, value: np.ndarray, n: int) -> np.ndarray:
    """Sparse (i,j,v) triplets to a dense symmetric matrix (C++ ``ijk2mat``)."""
    mat = np.zeros((n, n), dtype=np.float64)
    x = np.asarray(a_idx) - 1
    y = np.asarray(b_idx) - 1
    v = np.asarray(value, dtype=np.float64)
    mat[x, y] = v
    mat[y, x] = v
    return mat


def oe_zscore(m: np.ndarray) -> np.ndarray:
    """Column-wise z-score with ddof=1 (C++ ``oe2cor`` z-transform, lines 147-165).

    Columns with zero standard deviation become all-zero.
    """
    m = np.asarray(m, dtype=np.float64)
    mean = m.mean(axis=0)
    n = m.shape[0]
    sq = ((m - mean) ** 2).sum(axis=0)
    std = np.sqrt(sq / (n - 1))
    zm = np.zeros_like(m)
    nz = std > 0
    zm[:, nz] = (m[:, nz] - mean[nz]) / std[nz]
    return zm


def correlation_matrix(m: np.ndarray) -> np.ndarray:
    """Full column Pearson-correlation matrix via z-score then ``Zᵀ Z / (n-1)``.

    Mirrors ``oe2cor`` -> ``zmat`` (C++): ``block_cor = zm_blockᵀ zm_block / (nrow-1)``.
    """
    zm = oe_zscore(m)
    n = m.shape[0]
    return zm.T @ zm / (n - 1)


def _finalize(cor: np.ndarray) -> np.ndarray:
    """``mat2fbm``: unit diagonal + round to 5 dp (dchicf.r:17-33)."""
    cor = cor.copy()
    np.fill_diagonal(cor, 1.0)
    return np.round(cor, ROUND_DECIMALS)


def compartment_pca(oe: np.ndarray, n_pcs: int = 2):
    """Compartment PC scores from an O/E matrix (dchicf.r:212-316).

    Returns ``(pc, keep_mask)`` where ``pc`` has shape ``(n_kept, n_pcs)`` and
    ``keep_mask`` marks bins surviving the row-sum >= 3 blacklist.
    """
    oe = np.asarray(oe, dtype=np.float64)

    # blacklist: drop bins whose O/E row-sum < 3 (dchicf.r:219-224)
    keep = oe.sum(axis=1) >= 3
    m = oe[np.ix_(keep, keep)]

    # first correlation: C1 = corr(O/E), then mat2fbm finalize
    c1 = _finalize(correlation_matrix(m))

    # second correlation on C1; keep its z-scored matrix for the projection
    z2 = oe_zscore(c1)
    c2 = _finalize(correlation_matrix(c1))

    # Top-k eigenvectors of C2 via symmetric eigendecomposition.
    # ACCELERATION (E) exact identity: dcHiC computes the right singular vectors of the
    # symmetric PSD correlation matrix C2 via big_randomSVD. For a symmetric PSD matrix
    # the SVD and the eigendecomposition coincide (V == eigenvectors), so eigh of the
    # top-k eigenvalues returns the same subspace as svd up to per-vector sign, which the
    # downstream GC orientation fixes. PC1 (the dominant compartment axis) is well
    # separated, so the gated output is bit-identical to full SVD. eigh is ~3x faster than
    # full svd and avoids computing the discarded singular vectors. See ITERATION_LOG.md.
    w, vecs = np.linalg.eigh(c2)
    idx = np.argsort(w)[::-1][:n_pcs]   # largest eigenvalues first (PC1, PC2, ...)
    v = vecs[:, idx]                    # (n, n_pcs)

    pc = z2 @ v  # (n_kept, n_pcs)  -- dchicf.r:300
    return pc, keep


def orient_pc(pc: np.ndarray, gcc: np.ndarray) -> np.ndarray:
    """Flip each PC so it positively correlates with GC content (dchicf.r:526)."""
    pc = np.asarray(pc, dtype=np.float64)
    gcc = np.asarray(gcc, dtype=np.float64)
    out = pc.copy()
    for k in range(pc.shape[1]):
        r = np.corrcoef(pc[:, k], gcc)[0, 1]
        out[:, k] = np.sign(r) * pc[:, k]
    return out


def call_compartments(
    a_idx: np.ndarray,
    b_idx: np.ndarray,
    weight: np.ndarray,
    pos: np.ndarray,
    gcc: np.ndarray,
    n_bins: int,
    resolution: int,
    n_pcs: int = 2,
    count_thr: float = 0.0,
    minexpcc: float = 0.0,
):
    """End-to-end Stage A: raw contacts -> GC-oriented compartment PCs.

    Returns ``(pc_oriented, keep_mask)``. The A/B compartment score is ``pc[:, 0]``.
    """
    a, b, oe_val = observed_over_expected(
        a_idx, b_idx, weight, pos, n_bins, resolution, count_thr, minexpcc
    )
    oe = ijk2mat(a, b, oe_val, n_bins)
    pc, keep = compartment_pca(oe, n_pcs=n_pcs)
    pc = orient_pc(pc, np.asarray(gcc)[keep])
    return pc, keep

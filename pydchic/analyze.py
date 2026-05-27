"""High-level ``analyze`` driver — the Python equivalent of

    Rscript dchicf.r --file <input.txt> --pcatype analyze --dirovwt T --diffdir <name>

Assembles the per-replicate selected-PC bedGraphs produced by a prior ``cis``/``select`` run,
reproduces dcHiC's differential-compartment statistic, and writes the combined result.

Faithful to ``pcanalyze`` (dchicf.r:808): quantile-normalise across ALL replicates per
chromosome (dchicf.r:968), compute the robust-covariance Mahalanobis p-value per chromosome,
then adjust genome-wide. The per-bin statistic core is parity-validated against the R
reference; multiple-testing uses Benjamini-Hochberg (dcHiC uses IHW when replicates exist —
see MATH.md; both preserve the −log10 p ranking).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from .differential import bh_adjust, differential_compartments, normalize_quantiles

__all__ = ["read_input_file", "analyze"]


def read_input_file(input_file: str) -> pd.DataFrame:
    """Parse a dcHiC input.txt (matrix, bed, prefix, prefix.master[, ...])."""
    df = pd.read_csv(input_file, sep=r"\s+", header=None, comment="#")
    cols = ["mat", "bed", "prefix", "prefix_master"]
    df = df.iloc[:, : len(cols)]
    df.columns = cols[: df.shape[1]]
    return df


def _pc_dir(workdir: str, prefix: str, pc_type: str) -> Path:
    return Path(workdir) / f"{prefix}_pca" / f"{pc_type}_pca" / f"{prefix}_mat"


def _global_quantile_normalize(mat: np.ndarray) -> np.ndarray:
    """QN across all replicate columns + mean-abs scaling (dchicf.r:968-969)."""
    qn = normalize_quantiles(mat)
    colmeans = np.abs(qn).mean(axis=0)
    flat = qn.flatten(order="F")                       # R column-major recycling
    return (flat / np.resize(colmeans, flat.size)).reshape(qn.shape, order="F")


def analyze(
    input_file: str,
    workdir: str = ".",
    diffdir: str = "sample",
    pc_type: str = "intra",
    rzscore: float = 2.0,
    szscore: float = 0.0,
    refine: bool = True,
    rconf: float = 0.90,
    fdr_thr: float = 0.1,
    seed: int = 123,
    write: bool = True,
) -> pd.DataFrame:
    """Differential compartment analysis across all samples in ``input_file``.

    Equivalent to ``dchicf.r --pcatype analyze --diffdir <diffdir>``. Requires that a prior
    ``cis``/``select`` run produced ``<prefix>_pca/<pc_type>_pca/<prefix>_mat/<chr>.pc.bedGraph``
    for every replicate. Returns a genome-wide bins table with ``sample_maha``, ``pval``,
    ``padj`` (and the per-condition mean PC). When ``write`` is True, writes
    ``DifferentialResult/<diffdir>/fdr_result/differential.<pc_type>_sample_combined.pcQnm.bedGraph``
    and the ``.Filtered.`` subset (``padj < fdr_thr``) under ``workdir``.
    """
    inp = read_input_file(input_file)
    cond_of = dict(zip(inp["prefix"], inp["prefix_master"]))
    cond_order = list(dict.fromkeys(inp["prefix_master"]))               # first-occurrence order
    reps = [p for c in cond_order for p in inp["prefix"] if cond_of[p] == c]
    cond_labels = [cond_of[p] for p in reps]

    first_dir = _pc_dir(workdir, reps[0], pc_type)
    chroms = sorted(
        f[: -len(".pc.bedGraph")] for f in os.listdir(first_dir) if f.endswith(".pc.bedGraph")
    )

    per_chrom = []
    for chrom in chroms:
        frames = {}
        for p in reps:
            f = _pc_dir(workdir, p, pc_type) / f"{chrom}.pc.bedGraph"
            d = pd.read_csv(f, sep="\t", header=None, names=["chr", "start", "end", p])
            d.index = d["chr"].astype(str) + "_" + d["start"].astype(str) + "_" + d["end"].astype(str)
            frames[p] = d
        common = set.intersection(*[set(frames[p].index) for p in reps])
        common = [k for k in frames[reps[0]].index if k in common]       # preserve genomic order
        if len(common) <= len(cond_order):
            continue
        bed = frames[reps[0]].loc[common, ["chr", "start", "end"]].reset_index(drop=True)
        mat = np.column_stack([frames[p].loc[common, p].to_numpy(float) for p in reps])

        qn = _global_quantile_normalize(mat)
        res = differential_compartments(
            qn, cond_labels, rzscore=rzscore, szscore=szscore,
            refine=refine, rconf=rconf, seed=seed,
        )
        out = bed.copy()
        for c in cond_order:
            cols = [i for i, lab in enumerate(cond_labels) if lab == c]
            out[c] = qn[:, cols].mean(axis=1)
        out["sample_maha"] = res["sample_maha"]
        out["pval"] = res["pval"]
        per_chrom.append(out)

    combined = pd.concat(per_chrom, ignore_index=True).dropna(subset=["pval"])
    combined["padj"] = bh_adjust(combined["pval"].to_numpy())             # genome-wide BH

    if write:
        outdir = Path(workdir) / "DifferentialResult" / diffdir / "fdr_result"
        outdir.mkdir(parents=True, exist_ok=True)
        base = outdir / f"differential.{pc_type}_sample_combined.pcQnm.bedGraph"
        combined.to_csv(base, sep="\t", index=False)
        filt = combined[combined["padj"] < fdr_thr]
        filt.to_csv(str(base).replace(".bedGraph", ".Filtered.bedGraph"), sep="\t", index=False)
        print(f"[analyze] wrote {base} ({len(combined)} bins, {len(filt)} at padj<{fdr_thr})")
    return combined

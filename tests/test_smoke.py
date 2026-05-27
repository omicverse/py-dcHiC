"""Fast 'does it run at all' checks — no R reference required."""

from pathlib import Path

import numpy as np

import pydchic
from pydchic import (
    analyze,
    call_compartments,
    compartment_pca,
    differential_compartments,
    normalize_quantiles,
    quantile_normalize_by_condition,
)


def test_imports_and_version():
    assert pydchic.__version__
    assert hasattr(pydchic, "dcHiC")


def test_normalize_quantiles_equalizes_distributions():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(100, 3)) + np.array([0.0, 5.0, -3.0])
    qn = normalize_quantiles(a)
    # after QN every column shares the same sorted distribution
    s = np.sort(qn, axis=0)
    assert np.allclose(s[:, 0], s[:, 1], atol=1e-6)
    assert np.allclose(s[:, 0], s[:, 2], atol=1e-6)


def test_compartment_pca_recovers_checkerboard():
    # block-structured O/E -> PC1 separates the two compartments
    n = 40
    comp = np.r_[np.ones(20), -np.ones(20)]
    oe = 1.0 + 0.5 * np.outer(comp, comp) + 5.0
    pc, keep = compartment_pca(oe, n_pcs=2)
    assert pc.shape[0] == keep.sum()
    # PC1 should correlate (in magnitude) with the compartment label
    r = abs(np.corrcoef(pc[:, 0], comp[keep])[0, 1])
    assert r > 0.9


def test_differential_runs():
    rng = np.random.default_rng(1)
    bins, reps = 80, 6
    pc = rng.normal(size=(bins, reps))
    conditions = ["A", "A", "A", "B", "B", "B"]
    # inject a differential signal in some bins
    pc[:10, 3:] += 4.0
    qn = quantile_normalize_by_condition(pc, conditions)
    res = differential_compartments(qn, conditions, seed=123)
    assert res["pval"].shape == (bins,)
    assert np.all((res["padj"] >= 0) & (res["padj"] <= 1))
    # the injected bins should be among the most significant
    assert res["pval"][:10].mean() < res["pval"][10:].mean()


def test_analyze_driver_end_to_end(tmp_path):
    """`analyze` reads a dcHiC-style _pca tree and writes a combined differential bedGraph."""
    rng = np.random.default_rng(3)
    samples = [("A_1", "A"), ("A_2", "A"), ("B_1", "B"), ("B_2", "B")]
    chroms = {"chr1": 60, "chr2": 50}
    res_bp = 100_000
    for prefix, _ in samples:
        d = tmp_path / f"{prefix}_pca" / "intra_pca" / f"{prefix}_mat"
        d.mkdir(parents=True)
        for chrom, n in chroms.items():
            base = rng.normal(size=n)
            starts = np.arange(n) * res_bp
            with open(d / f"{chrom}.pc.bedGraph", "w") as fh:
                for i in range(n):
                    fh.write(f"{chrom}\t{starts[i]}\t{starts[i]+res_bp}\t{base[i]:.6f}\n")
    inp = tmp_path / "input.txt"
    inp.write_text("\n".join(f"m\tb\t{p}\t{c}" for p, c in samples) + "\n")

    out = analyze(str(inp), workdir=str(tmp_path), diffdir="test", write=True)
    assert set(out["chr"]) == set(chroms)
    assert {"sample_maha", "pval", "padj", "A", "B"} <= set(out.columns)
    assert np.all((out["padj"] >= 0) & (out["padj"] <= 1))
    assert (tmp_path / "DifferentialResult" / "test" / "fdr_result"
            / "differential.intra_sample_combined.pcQnm.bedGraph").exists()

# py-dcHiC

**Pure-Python port of [dcHiC](https://github.com/ay-lab/dcHiC)** — differential compartment
analysis of Hi-C datasets. Reproduces dcHiC's A/B-compartment calling and differential-compartment
statistic in NumPy/SciPy/scikit-learn, validated against the R reference under a pre-registered
parity gate (omicverse-rebuildr protocol).

> **Scope.** Compartment calling (`cis` + `select`) and differential analysis (`analyze`).
> Out of scope: subcompartment HMM (`depmixS4`), fithic differential loops, IGV/HTML
> visualization, gene enrichment, trans/inter compartments. See `AUDIT.md`.

## Parity with R dcHiC

| Output | dcHiC step | Parity class | Threshold | **Measured** |
|---|---|---|---|---|
| `compartment_pc1` | compartment calling | ordinal (Pearson \|r\|) | ≥ 0.99 | **1.0000** (max \|Δ\| = 1.6e-14) |
| `differential_padj` | differential analysis | inference (Spearman −log10p; top-50 Jaccard) | ≥ 0.90; ≥ 0.7 | **0.981; 0.961** |

Compartment PC1 is **bit-equivalent** to dcHiC's C++ kernel (`functionsdchic::oe2cor` + SVD).
The differential ranking matches dcHiC at Spearman 0.98 — all of R's `padj<0.1` bins are
recovered — despite the unavoidable robust-covariance / IHW library substitutions (see `MATH.md`).

## Install

```bash
pip install pydchic
```

## Quickstart

```python
import numpy as np
import pydchic

# --- Compartment calling: raw intra-chromosomal contacts -> A/B score (PC1) ---
pc, keep = pydchic.call_compartments(
    a_idx, b_idx, weight,   # sparse (bin_i, bin_j, count); 1-based bin ids
    pos, gcc,               # bin start coords (bp); GC-content track (for orientation)
    n_bins, resolution, n_pcs=2,
)
compartment_score = pc[:, 0]    # >0 = A compartment, <0 = B compartment

# --- Differential analysis across conditions (each with replicates) ----------
pc_raw = ...                    # bins x replicates, columns grouped by condition
conditions = ["ES","ES","NPC","NPC","CN","CN"]
qn  = pydchic.quantile_normalize_by_condition(pc_raw, conditions)
res = pydchic.differential_compartments(qn, conditions)   # -> pval, padj, sample_maha
```

Class API:

```python
from pydchic import dcHiC
model = dcHiC(resolution=100_000).differential(pc_raw, conditions)
model.differential_.sort_values("padj").head()
```

## What's included (Python ⇄ R map)

| Python | dcHiC R | Stage |
|---|---|---|
| `call_compartments` | `ijk2matfunc_cis` + `pcselect` orientation | compartment calling |
| `observed_over_expected` | `readfilesintra` O/E model | compartment calling |
| `compartment_pca` | `oe2cor` ×2 + SVD + projection | compartment calling |
| `orient_pc` | `pcselectioncore` GC sign-flip | select |
| `normalize_quantiles` | `limma::normalizeQuantiles` | differential |
| `calcen` | `calcen` | differential |
| `differential_compartments` | `pcanalyze` (Mahalanobis + chi-sq) | differential |

## Reproducing the R parity check

```bash
pip install -e ".[dev]"
export R_TEST_ENV=/path/to/conda/env/with/dcHiC   # R + functionsdchic, limma, robust, jsonlite
export PYTHON_TEST_ENV=/path/to/this/python/env
pytest -q            # builds the fixture, runs the R reference + the port, asserts both gates
```

The pre-registered gate lives in [`data/manifest.yaml`](data/manifest.yaml); the R reference
driver is [`tests/r_reference_driver.R`](tests/r_reference_driver.R).

## Notebooks (pre-executed under `examples/`)

- `compare_R_vs_Python.ipynb` — pipeline-level parity vs R, one plot per output.
- `tutorial_dchic.ipynb` — Python-only tour of the public API.
- `function_by_function_R_parity.ipynb` — R⇄Python per-function dictionary with live comparisons.
- `evolution.ipynb` — per-iteration acceleration narrative (`svd → eigh`, exact identity).

## Relationship to omicverse

Produced under the [omicverse-rebuildr](https://github.com/omicverse) reference-driven porting
protocol. Reconstruction's goal is **identical** output to R, not "better" — see
`RECONSTRUCTION_REPORT.md`.

## Citation

If you use this port, please cite the original dcHiC paper:

> Chakraborty A., Wang J.G., Ay F. *dcHiC detects differential compartments across multiple
> Hi-C datasets.* **Nature Communications** 13, 6827 (2022).

## License

GPL-3.0. The dcHiC tool is MIT, but its `functionsdchic` C++ kernels (re-implemented here from
scratch) are GPL (≥2); the port matches the more restrictive license. See `LICENSE`.

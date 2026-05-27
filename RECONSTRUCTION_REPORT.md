# Reconstruction Report — py-dcHiC

## 1. Identity

| Field | Value |
|---|---|
| Package | `py-dcHiC` (PyPI `pydchic`, import `pydchic`) |
| Upstream | [ay-lab/dcHiC](https://github.com/ay-lab/dcHiC) @ `144b4df` (2024-01-10); `functionsdchic` 1.1 |
| Paper | Chakraborty, Wang, Ay. *Nat. Commun.* 13, 6827 (2022) |
| Scope | compartment calling (`cis` + `select`) + differential analysis (`analyze`) |
| Algorithm classes | (1) ordinal — compartment PC1; (2) inference — differential p-values |
| Thresholds (read-only) | PC1 Pearson \|r\| ≥ 0.99; differential Spearman(−log10p) ≥ 0.90 & top-50 Jaccard ≥ 0.7 |
| Final parity | PC1 \|r\| = **1.0** (max \|Δ\| = 1.6e-14); differential Spearman = **0.981**, Jaccard = **0.961** |
| Audit class | **A** (translation-only parity; one exact-identity acceleration) |
| LOC (port) | ~330 (`pydchic/`) |
| Speedup | Stage-A decomposition 1.83× on fixture (3.0–3.2× at n=1000–2000) via `svd→eigh` |

## 2. R function coverage audit

Full table in [`AUDIT.md`](AUDIT.md). Summary: **4/7** `functionsdchic` C++ helpers and **8/25**
`dchicf.r` functions ported — i.e. the complete numerical path for the two in-scope stages. The
remaining 17 driver functions are out of scope (HMM subcompartments, fithic loops, HTML/IGV viz,
gene enrichment, trans) per the registered scope and the protocol's "when NOT to use" guidance
(S4-heavy / GUI / plotting).

### Dependencies reused from omicverse

| R dep | omicverse port | Reused as | Saved work |
|---|---|---|---|
| (none) | — | — | greenfield; no upstream dep had an omicverse mirror (see `DISCOVERY.md`) |

## 3. Parity evidence

| Output | Fixture | Metric | Value | Threshold | Pass |
|---|---|---|---|---|---|
| `compartment_pc1` | synthetic 200-bin Hi-C (compartment + decay) | Pearson \|r\| | 1.0000 (max\|Δ\|=1.6e-14) | ≥ 0.99 | ✅ |
| `differential_padj` | dcHiC demo ESC/NPC/CN chr19, 583 bins × 12 reps | Spearman(−log10p) | 0.981 | ≥ 0.90 | ✅ |
| `differential_padj` | (same) | top-50 Jaccard | 0.961 | ≥ 0.7 | ✅ |

Significant-bin recovery: R reports 89 bins at `padj<0.1`; the port recovers **all 89** (Jaccard 0.967).

Reproducible reference command:
```bash
Rscript tests/r_reference_driver.R data/fixture_dchic.json data/reference_output.json
python  tests/_run_candidate.py    data/fixture_dchic.json data/candidate_output.json
pytest -q   # asserts both gates from data/manifest.yaml
```
The R reference invokes the **real** dcHiC primitives: `functionsdchic::oe2cor` (C++),
`limma::normalizeQuantiles`, `robust::covRob`, `stats::mahalanobis`.

## 4. Acceleration evidence

Two-panel evolution figure: [`examples/evolution.png`](examples/evolution.png); per-iteration
narrative: [`examples/evolution.ipynb`](examples/evolution.ipynb); log: [`ITERATION_LOG.md`](ITERATION_LOG.md).

| iter | action | admissibility | mean time | speedup | parity | status |
|---|---|---|---|---|---|---|
| 0 | baseline (full `svd`) | — | 5.91 ms | 1× | \|r\|=1.0 | — |
| 1 | `svd → eigh` (top-k) | **(E) exact identity** | 3.23 ms | 1.83× | \|r\|=1.0 | ACCEPT |

(E) proof in [`MATH.md §1`](MATH.md): `C2` is symmetric PSD ⇒ SVD = eigendecomposition; sign fixed
by GC orientation; PC1 dominant/non-degenerate ⇒ gated output bit-identical. No parity dip.

## 5. Code quality audit

| Check | Status |
|---|---|
| `pip install -e .` in clean env | ✅ |
| `pytest -q` | ✅ 5 passed (4 smoke + 1 parity gate) |
| 4 mandatory notebooks pre-executed (0 errors) | ✅ compare / tutorial / function-by-function / evolution |
| `evolution.png` rendered from `ITERATION_LOG.md` | ✅ |
| License compatible (GPL-3 ≥ upstream functionsdchic GPL≥2) | ✅ |
| Version pinned `0.1.0` | ✅ |

## 6. Known limitations

- **Robust covariance**: dcHiC's `robust::covRob` (Donoho–Stahel default) is approximated by
  `sklearn.MinCovDet` (FastMCD). Differential **p-values match in ranking, not element-wise**
  (inference-class gate); the gate is never widened to force element-wise agreement.
- **Multiple testing**: dcHiC's `IHW` (with replicate covariate) is replaced by Benjamini–Hochberg
  — no faithful Python IHW exists. Both monotone in the p-value, so −log10p ranking is preserved.
- **PC selection**: GC-oriented **PC1** is gated; dcHiC's multi-sample hierarchical-clustering PC
  *selection* (choosing PC1 vs PC2 per chromosome across replicates) is not part of the gate.
- **Out of scope**: trans/inter compartments, subcompartment HMM, fithic loops, viz, enrichment.

## 7. Integration into omicverse

- Vendor location: `omicverse/py-dcHiC` (proposed); import surface `pydchic`.
- Public API: `dcHiC` class + functional mirrors (`call_compartments`, `compartment_pca`,
  `differential_compartments`, `normalize_quantiles`, `calcen`, …).
- Tutorial slot: `examples/tutorial_dchic.ipynb`.

## 8. Sign-off

| Field | Value |
|---|---|
| Author | omicverse-rebuildr agent (Claude) |
| Date | 2026-05-26 |
| Final audit class | **A** |
| Both pre-registered gates | ✅ PASS |

# R function coverage audit — py-dcHiC

> `engine.r_function_audit` parses a package `NAMESPACE` + `R/`. dcHiC ships its algorithm
> as a **standalone driver script** (`dchicf.r`, 2976 lines) plus a thin Rcpp helper package
> (`functionsdchic`), so the automatic exported-symbol audit returns 0 (no `export()` lines;
> the helper uses `exportPattern`). This file is therefore a manual, line-referenced audit of
> both surfaces, with ported / out-of-scope status per the registered scope
> (compartment calling + differential analysis).

## Coverage summary

| Surface | Ported | In scope | Out of scope | Total |
|---|---|---|---|---|
| `functionsdchic` C++ helpers | 4 | 4 | 3 | 7 |
| `dchicf.r` functions | 8 | 8 | 17 | 25 |

"Ported" = re-implemented in pure Python with parity coverage. "Out of scope" = excluded by
the registered scope (viz/HTML, fithic loops, HMM subcompartments, gene enrichment, trans).

## `functionsdchic` C++ helpers (src/functionsdchic.cpp)

| R/C++ function | Python equivalent | Status | Note |
|---|---|---|---|
| `ijk2mat` | `compartment.ijk2mat` | ✅ ported | sparse (i,j,v) → dense symmetric matrix |
| `oe2cor` | `compartment.oe_zscore` + `correlation_matrix` | ✅ ported | column z-score (ddof=1) then `Zᵀ Z/(n-1)` |
| `zmat` | `compartment.correlation_matrix` | ✅ ported | block correlation tiling collapses to one matmul |
| `eigenMapMatMult` | `numpy` `@` | ✅ ported | `Aᵀ B` |
| `fbm2mat` | n/a | ⬜ out of scope | FBM↔matrix glue (bigstatsr file-backed matrices) |
| `transmat` | n/a | ⬜ out of scope | trans (inter-chromosomal) matrix assembly |
| `createtransijk` | n/a | ⬜ out of scope | trans triplet extraction |

## `dchicf.r` driver functions

| R function | dchicf.r | Python equivalent | Status |
|---|---|---|---|
| `expectedInteraction` | 204 | `compartment.expected_interaction` | ✅ ported |
| `readfilesintra` (O/E model) | 319 | `compartment.observed_over_expected` | ✅ ported (numeric core) |
| `ijk2matfunc_cis` | 212 | `compartment.compartment_pca` | ✅ ported |
| `mat2fbm` | 17 | `compartment._finalize` | ✅ ported (unit diag + round-5) |
| `pcselectioncore` (orientation) | 506 | `compartment.orient_pc` | ✅ ported (GC sign-flip) |
| `mdweight` / `calcen` | 737 / 751 | `differential.calcen` | ✅ ported |
| `calchclust` | 789 | `differential` (distclust=-1 path) | ✅ ported (default no-op) |
| `pcanalyze` (differential stat) | 808 | `differential.differential_compartments` | ✅ ported (numeric core) |
| `pcselectioncore` (hclust PC selection) | 559 | — | ⬜ partial: GC-oriented PC1 ported; multi-sample hclust PC *selection* not gated |
| `pcselect` (genome goldenpath I/O) | 662 | — | ⬜ out of scope (bedtools/curl GC+TSS track generation) |
| `readfilesinter` | 166 | — | ⬜ out of scope (trans/inter) |
| `extractTrans` | 37 | — | ⬜ out of scope (trans) |
| `subcompartment_level` | 1267 | — | ⬜ out of scope (HMM) |
| `hmmsegment` | 1316 | — | ⬜ out of scope (depmixS4 HMM) |
| `subcompartment` | 1394 | — | ⬜ out of scope (HMM) |
| `callfithiC` | 1488 | — | ⬜ out of scope (fithic loops) |
| `formatconversion` | 1523 | — | ⬜ out of scope (fithic loops) |
| `fithicformat` | 1595 | — | ⬜ out of scope (fithic loops) |
| `compartmentLoop` | 1676 | — | ⬜ out of scope (differential loops) |
| `htmlheader` / `scriptbody` / `htmlbody` | 1784 / 1807 / 1832 | — | ⬜ out of scope (HTML viz) |
| `generateTrackfiles` | 2050 | — | ⬜ out of scope (IGV tracks) |
| `geneEnrichment` | 2365 | — | ⬜ out of scope (functional enrichment) |
| `datadownload` / `.findExecutable` | 478 / 468 | — | ⬜ out of scope (CLI/env glue) |

## Dependencies reused from omicverse

| R dep | omicverse port | Reused as | Saved work |
|---|---|---|---|
| (none) | — | — | greenfield port; no upstream dep had an omicverse mirror |

## Deviations from the R reference

| Upstream | Port | Reason |
|---|---|---|
| `bigstatsr::big_randomSVD` | `numpy.linalg.eigh` (top-k) | exact identity for symmetric PSD C2; ~3× faster (see ITERATION_LOG.md iter 1) |
| `robust::covRob` (Donoho-Stahel) | `sklearn.covariance.MinCovDet` (FastMCD) | no faithful Python port of Donoho-Stahel; both robust; ranking preserved (Spearman 0.98) |
| `IHW::ihw` | Benjamini-Hochberg | no faithful Python IHW; gate is rank-based on −log10 p (see MATH.md) |

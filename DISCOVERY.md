# Discovery — py-dcHiC

> Committed BEFORE any algorithmic code (Phase 0.5). Records the decision to start
> this port and which omicverse-org Python ports get reused as dependencies.

## 1. Is this package already ported?

Output of `python -m engine.discover_omicverse_deps --check dcHiC` (and `--check dchic`):

```
## Discovery — `dcHiC`

**No existing omicverse port found.** Safe to start a new port.
```

**Decision**: START_PORT — no `omicverse/py-dcHiC` (or `py-dchic`) exists.

## 2. Dependency audit

The upstream `functionsdchic` `DESCRIPTION` only declares the Rcpp/RcppEigen build
deps (it is a thin C++ helper package). The *algorithmic* dependencies live in the
`dchicf.r` driver via `library()` / `::` calls. Audited from `dcHiC-ref/dchicf.r`:

```
library(): Rcpp bench bigstatsr data.table depmixS4 functionsdchic hashmap limma optparse parallel
::       : data.table(29) parallel(16) bigstatsr(16) hashmap(10) R.utils(10)
           functionsdchic(9) rjson(6) robust(5) networkD3(4) limma(4)
           htmlwidgets(4) depmixS4(4) colorspace(4) IHW(4) bench(3)
```

`engine.discover_omicverse_deps --description` reports **0 of 4** declared deps have a
omicverse mirror (the DESCRIPTION only lists Rcpp/RcppEigen). None of the algorithmic
deps below have an omicverse Python mirror either — all map to native scientific Python.

## 3. Decisions per R dep (scoped to compartment-calling + differential stages)

| R dep | omicverse match | Decision | Python equivalent | Note |
|---|---|---|---|---|
| `functionsdchic` (C++: `oe2cor`, `zmat`, `eigenMapMatMult`, `ijk2mat`) | — | native-python (reimplement) | `numpy` | GPL≥2 kernels — reimplemented from scratch, not linked |
| `bigstatsr` (`big_SVD` PCA) | — | native-python | `numpy.linalg.eigh` / `scipy.sparse.linalg` | correlation-matrix eigendecomposition |
| `limma` (`normalizeQuantiles`) | — | native-python | reimplement quantile normalization | small, well-specified |
| `robust` (`covRob`) | — | native-python | `sklearn.covariance.MinCovDet` | robust covariance for Mahalanobis |
| `IHW` | — | native-python | `statsmodels` BH + weighting | independent hypothesis weighting |
| `data.table` | — | native-python | `pandas` | I/O + group ops |
| `parallel` | — | native-python | `joblib` / `multiprocessing` | not algorithmic |
| `depmixS4` (HMM subcompartments) | — | out-of-scope | — | S4-heavy; excluded per protocol |
| `networkD3` / `htmlwidgets` / `colorspace` | — | out-of-scope | — | HTML/IGV visualization |
| `hashmap` / `R.utils` / `bench` / `rjson` / `optparse` | — | out-of-scope | — | utility / CLI / benchmark glue |

## 4. Reusable work saved

| Reused omicverse port | Approx. LOC saved | Notes |
|---|---|---|
| (none) | 0 | No upstream dep has an existing omicverse mirror. |

**Total**: 0 LOC reused — this is a greenfield port. All numerical kernels reimplemented
on numpy/scipy/sklearn/pandas.

## 5. New ports surfaced

- `IHW` (Independent Hypothesis Weighting, Bioconductor) — used by several Bioc DE tools;
  no Python equivalent. Worth a future standalone port. Suggested rank: T3.
- `depmixS4` (HMM) — general dependency-mixture HMM; out of scope here but broadly useful.

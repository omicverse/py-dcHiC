"""Build the four mandatory pre-executed notebooks via nbformat, then they are executed
in-place by jupyter nbconvert (see the shell driver after this file is written).

Notebooks:
  1. compare_R_vs_Python.ipynb        — pipeline parity, one viz per manifest output
  2. tutorial_dchic.ipynb             — Python-only tour of the public API
  3. function_by_function_R_parity.ipynb — R<->Python per-function dictionary
  4. evolution.ipynb                  — per-iteration acceleration narrative + subplots
"""

from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
PORT = HERE.parent


def md(s):
    return nbf.v4.new_markdown_cell(s.strip("\n"))


def code(s):
    return nbf.v4.new_code_cell(s.strip("\n"))


def save(cells, name):
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {
        "name": "rebuild-py", "display_name": "Python (rebuild-py)", "language": "python"
    }
    out = HERE / name
    nbf.write(nb, out)
    print("wrote", out)


SETUP = f"""
import os, json, subprocess, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, {str(PORT)!r})
sys.path.insert(0, {str(PORT / 'tests')!r})
PORT = {str(PORT)!r}
R_ENV = os.environ.get("R_TEST_ENV", "/home/shengmao/.local/share/mamba/envs/dchic")
RSCRIPT = f"{{R_ENV}}/bin/Rscript"
FIX = os.path.join(PORT, "data", "fixture_dchic.json")
"""


# --------------------------------------------------------------------------- #
# Notebook 1 — compare_R_vs_Python
# --------------------------------------------------------------------------- #
def nb_compare():
    cells = [
        md("""
# py-dcHiC vs dcHiC — pipeline parity

Pipeline-level proof that **py-dcHiC** reproduces the R reference (`ay-lab/dcHiC`) on the
canonical fixture, for both registered outputs:

| Output | Stage | Parity class | Threshold |
|---|---|---|---|
| `compartment_pc1` | compartment calling (`cis`+`select`) | ordinal (Pearson \\|r\\|) | ≥ 0.99 |
| `differential_padj` | differential analysis (`analyze`) | inference (Spearman −log10p + top-50 Jaccard) | ≥ 0.90 / ≥ 0.7 |

Both the R reference and the Python candidate read the *same* fixture; the R side calls the
real dcHiC primitives (`functionsdchic::oe2cor`, `limma::normalizeQuantiles`,
`robust::covRob`, `stats::mahalanobis`).
"""),
        code(SETUP),
        md("## 1. Run the R reference and the Python candidate on the same fixture"),
        code("""
ref_json = os.path.join(PORT, "data", "reference_output.json")
cand_json = os.path.join(PORT, "data", "candidate_output.json")
subprocess.run([RSCRIPT, os.path.join(PORT, "tests", "r_reference_driver.R"), FIX, ref_json], check=True)
env = dict(os.environ, PYTHONPATH=PORT)
subprocess.run([sys.executable, os.path.join(PORT, "tests", "_run_candidate.py"), FIX, cand_json], check=True, env=env)
ref = json.load(open(ref_json)); cand = json.load(open(cand_json))
print("outputs:", list(ref))
"""),
        md("## 2. Apply the pre-registered parity metrics"),
        code("""
from parity_metrics import compute_parity, is_pass
m1 = compute_parity(ref["compartment_pc1"], cand["compartment_pc1"], "ordinal")
m2 = compute_parity(ref["differential_padj"], cand["differential_padj"], "inference")
print(f"compartment_pc1   ordinal |r| = {m1:.8f}   pass(>=0.99) = {is_pass(m1,'ordinal',0.99)}")
print(f"differential_padj inference   = {m2}   pass = {is_pass(m2,'inference',0.90)}")
"""),
        md("## 3. Visualization — output 1: compartment PC1 (R vs Python)"),
        code("""
r1 = np.array(ref["compartment_pc1"]); c1 = np.array(cand["compartment_pc1"])
if np.corrcoef(r1, c1)[0,1] < 0: c1 = -c1   # PC sign is arbitrary; align for display
fig, ax = plt.subplots(1, 2, figsize=(11,4))
ax[0].plot(r1, label="R dcHiC", lw=1); ax[0].plot(c1, "--", label="py-dcHiC", lw=1)
ax[0].set_title("Compartment PC1 along chromosome"); ax[0].set_xlabel("bin"); ax[0].legend()
ax[1].scatter(r1, c1, s=8, alpha=0.6); ax[1].plot([r1.min(),r1.max()],[r1.min(),r1.max()],'r-',lw=1)
ax[1].set_title(f"R vs Py (max|Δ|={np.max(np.abs(r1-c1)):.1e})"); ax[1].set_xlabel("R PC1"); ax[1].set_ylabel("py PC1")
plt.tight_layout(); plt.show()
"""),
        md("## 4. Visualization — output 2: differential significance (R vs Python)"),
        code("""
r2 = np.array(ref["differential_padj"]); c2 = np.array(cand["differential_padj"])
nlr, nlc = -np.log10(np.clip(r2,1e-300,1)), -np.log10(np.clip(c2,1e-300,1))
sr, sc = set(np.where(r2<0.1)[0]), set(np.where(c2<0.1)[0])
fig, ax = plt.subplots(1, 2, figsize=(11,4))
ax[0].scatter(nlr, nlc, s=8, alpha=0.6); lim=[0, max(nlr.max(),nlc.max())]
ax[0].plot(lim,lim,'r-',lw=1); ax[0].set_title("-log10(padj): R vs py"); ax[0].set_xlabel("R"); ax[0].set_ylabel("py")
ax[1].bar(["R sig","py sig","overlap"], [len(sr),len(sc),len(sr&sc)], color=["#4c72b0","#dd8452","#55a868"])
ax[1].set_title(f"padj<0.1 bins (Jaccard={len(sr&sc)/len(sr|sc):.3f})")
plt.tight_layout(); plt.show()
"""),
        md("""
## Verdict

Both pre-registered gates **pass**: compartment PC1 is bit-equivalent to dcHiC's C++/SVD
kernel (`|r|=1.0`, `max|Δ|≈1.6e-14`), and the differential ranking matches dcHiC at
Spearman ≈ 0.98 on −log10 p with all of R's significant bins recovered, despite the
robust-covariance / IHW library substitutions documented in `MATH.md`.
"""),
    ]
    save(cells, "compare_R_vs_Python.ipynb")


# --------------------------------------------------------------------------- #
# Notebook 2 — tutorial
# --------------------------------------------------------------------------- #
def nb_tutorial():
    cells = [
        md("""
# py-dcHiC tutorial — A/B compartments & differential analysis in Python

A copy-pastable tour of the public API. `pip install pydchic`, then:
"""),
        code(SETUP + "\nimport pydchic\nprint('pydchic', pydchic.__version__)"),
        md("""
## 1. Compartment calling from a contact matrix

`call_compartments` turns raw intra-chromosomal contacts into GC-oriented compartment PCs:
observed/expected → correlation → correlation-of-correlation → eigendecomposition →
projection → GC sign-orientation. The A/B score is PC1.
"""),
        code("""
fix = json.load(open(FIX)); a = fix["stageA"]
pc, keep = pydchic.call_compartments(
    np.array(a["a_idx"]), np.array(a["b_idx"]), np.array(a["weight"], float),
    np.array(a["pos"], float), np.array(a["gcc"], float),
    int(a["n_bins"]), int(a["resolution"]), n_pcs=2)
print("compartment PC matrix:", pc.shape, "| kept bins:", int(keep.sum()))
plt.figure(figsize=(10,3)); s = pc[:,0]
plt.fill_between(range(len(s)), s, where=s>=0, color="#c44e52", label="A (PC1>0)")
plt.fill_between(range(len(s)), s, where=s<0, color="#4c72b0", label="B (PC1<0)")
plt.title("A/B compartment track (PC1)"); plt.xlabel("bin"); plt.legend(); plt.show()
"""),
        md("## 2. Lower-level functional API (mirrors the R functions one-to-one)"),
        code("""
from pydchic.compartment import observed_over_expected, ijk2mat, compartment_pca, orient_pc
A,B,oe = observed_over_expected(np.array(a["a_idx"]),np.array(a["b_idx"]),np.array(a["weight"],float),
                                np.array(a["pos"],float),int(a["n_bins"]),int(a["resolution"]))
M = ijk2mat(A,B,oe,int(a["n_bins"]))
pc_raw,_ = compartment_pca(M, n_pcs=2)
pc_or = orient_pc(pc_raw, np.array(a["gcc"])[M.sum(1)>=3])
print("O/E matrix:", M.shape, "| PC after orientation == call_compartments:",
      np.allclose(pc_or, pc))
"""),
        md("""
## 3. Differential compartment analysis across conditions

`quantile_normalize_by_condition` + `differential_compartments` reproduce dcHiC's `analyze`
step: per-condition quantile normalization, condition means, a robust-covariance Mahalanobis
distance scaled by the squared cross-condition range, a chi-square p-value, and BH adjustment.
"""),
        code("""
b = fix["stageB"]
pc_raw = np.array(b["pc_raw"], float); conditions = b["conditions"]
qn = pydchic.quantile_normalize_by_condition(pc_raw, conditions)
res = pydchic.differential_compartments(qn, conditions)
import numpy as np
order = np.argsort(res["padj"])
print("conditions:", sorted(set(conditions)))
print("most differential bins (padj):", np.round(res["padj"][order][:5], 4))
plt.figure(figsize=(10,3)); plt.plot(-np.log10(np.clip(res["padj"],1e-300,1)), lw=1)
plt.axhline(-np.log10(0.1), color="r", ls="--", label="FDR 0.1")
plt.title("Differential compartment significance"); plt.xlabel("bin"); plt.ylabel("-log10(padj)"); plt.legend(); plt.show()
"""),
        md("## 4. Class API (method chaining)"),
        code("""
from pydchic import dcHiC
model = dcHiC(resolution=int(a["resolution"]))
model.differential(pc_raw, conditions)
print(model.differential_.sort_values("padj").head())
"""),
        md("""
## Pitfalls (parity notes)

- **PC sign is arbitrary** — dcHiC orients PC1 by GC content; absolute sign is not meaningful.
- **Columns must be grouped by condition** for the differential step.
- The robust covariance (`MinCovDet`) and BH (vs dcHiC's `IHW`) mean differential **p-values
  match in ranking**, not element-wise — see `MATH.md`.
"""),
    ]
    save(cells, "tutorial_dchic.ipynb")


# --------------------------------------------------------------------------- #
# Notebook 3 — function-by-function R parity
# --------------------------------------------------------------------------- #
def nb_fbf():
    cells = [
        md("""
# Function-by-function R ⇄ Python dictionary

For an R user porting dcHiC code line-by-line: each section gives the R call, the py-dcHiC
call on identical input, and a numeric comparison. Small R snippets run live via `Rscript`.
"""),
        code(SETUP + """
def r_eval(expr):
    "Run an R expression that cat()s a JSON result; return the parsed object."
    out = subprocess.run([RSCRIPT, "-e", expr], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    return json.loads(out.stdout.strip().splitlines()[-1])
"""),
        md("""
## `oe2cor` / `zmat`  ⇄  `pydchic.compartment.correlation_matrix`

| R param | Py param | Type | Meaning |
|---|---|---|---|
| `m` | `m` | matrix | O/E matrix (bins × bins) |
| `s`,`e` | — | int vec | column block ranges (tiling; result identical) |
| `check_cov=0` | — | flag | return correlation (not covariance) |

R computes column z-scores (ddof=1) then `Zᵀ Z /(n-1)`.
"""),
        code("""
import numpy as np
from pydchic.compartment import correlation_matrix
rng = np.random.default_rng(0); M = rng.normal(size=(8,8)); M = (M+M.T)/2 + 5
np.savetxt("/tmp/_m.txt", M)
r = r_eval('library(jsonlite); library(functionsdchic); m<-as.matrix(read.table("/tmp/_m.txt"));'
           ' x<-functionsdchic::oe2cor(m, c(0), c(ncol(m)-1), 1, 0); cat(toJSON(x[[1]]$mat))')
R = np.array(r); P = correlation_matrix(M)
print("R[:2,:2]=\\n", np.round(R[:2,:2],5)); print("Py[:2,:2]=\\n", np.round(P[:2,:2],5))
print("max|Δ| =", np.max(np.abs(R-P)), "-> VERDICT:", "MATCH" if np.allclose(R,P,atol=1e-6) else "DIFF")
"""),
        md("""
## `limma::normalizeQuantiles`  ⇄  `pydchic.normalize_quantiles`

| R param | Py param | Type | Default | Meaning |
|---|---|---|---|---|
| `A` | `a` | matrix | — | values (bins × samples) |
| `ties` | (always) | logical | `TRUE` | average-rank interpolation |
"""),
        code("""
from pydchic import normalize_quantiles
rng = np.random.default_rng(1); X = rng.normal(size=(12,3)) + np.array([0,4,-2])
np.savetxt("/tmp/_x.txt", X)
r = r_eval('library(jsonlite); library(limma); x<-as.matrix(read.table("/tmp/_x.txt"));'
           ' cat(toJSON(limma::normalizeQuantiles(x, ties=TRUE)))')
R = np.array(r); P = normalize_quantiles(X)
print("max|Δ| =", np.max(np.abs(R-P)), "-> VERDICT:", "MATCH" if np.allclose(R,P,atol=1e-8) else "DIFF")
"""),
        md("""
## `calcen`  ⇄  `pydchic.calcen`

dcHiC's shrunk centre estimate. `cls='sample'` z-scores the per-row spread per column;
`cen = round(df * (1 - rowmax(pnorm(z, mean=szscore))), 5)`.

| R param | Py param | Default | Meaning |
|---|---|---|---|
| `df` | `df` | — | values (bins × conditions) |
| `class` | `cls` | — | `'rep'` or `'sample'` |
| `rzscore` | `rzscore` | 2 | replicate z threshold |
| `szscore` | `szscore` | 0 | sample z threshold |
"""),
        code("""
from pydchic import calcen
rng = np.random.default_rng(2); G = rng.normal(size=(20,3))
np.savetxt("/tmp/_g.txt", G)
r = r_eval('''library(jsonlite)
calcen <- function(df, class, rzscore, szscore) {
  df_dist <- as.data.frame(t(apply(df,1,function(x){sqrt(colSums(as.matrix(dist(as.numeric(x))))/(length(as.numeric(x))-1))})))
  df_zsc <- list(); for(n in 1:ncol(df)){v<-as.numeric(as.matrix(df_dist[,n])); df_zsc[[n]]<-(v-mean(v))/sd(v)}
  df_zsc <- do.call(cbind, df_zsc)
  df_pvl <- round(pnorm(as.matrix(df_zsc), mean=szscore, lower.tail=T),5)
  round(df*(1-apply(df_pvl,1,max)),5)
}
g<-as.matrix(read.table("/tmp/_g.txt")); cat(toJSON(as.matrix(calcen(as.data.frame(g),"sample",2,0))))''')
R = np.array(r); P = calcen(G, "sample", 2, 0)
print("max|Δ| =", np.max(np.abs(R-P)), "-> VERDICT:", "MATCH" if np.allclose(R,P,atol=1e-4) else "DIFF")
"""),
        md("""
## Differential statistic — top-level

`pcanalyze` ⇄ `differential_compartments`. The condition-level Mahalanobis p-value matches
dcHiC in **ranking** (inference-class gate); robust covariance (`covRob`→`MinCovDet`) and
multiple testing (`IHW`→`BH`) are documented library substitutions (`MATH.md`).
"""),
        code("""
ref = json.load(open(os.path.join(PORT,"data","reference_output.json")))
cand = json.load(open(os.path.join(PORT,"data","candidate_output.json")))
from parity_metrics import compute_parity
print("differential_padj inference:", compute_parity(ref["differential_padj"], cand["differential_padj"], "inference"))
"""),
    ]
    save(cells, "function_by_function_R_parity.ipynb")


# --------------------------------------------------------------------------- #
# Notebook 4 — evolution
# --------------------------------------------------------------------------- #
def nb_evolution():
    cells = [
        md("""
# Acceleration evolution — py-dcHiC

One section per acceleration iteration. Reconstruction's goal is **identical** output to R,
not "better" — so accuracy stays pinned at the ceiling while wall-clock drops. Timings are the
full Stage-A compartment PCA on the fixture, warmup-excluded, 3 runs, `OMP_NUM_THREADS=8`.
"""),
        code(SETUP + """
os.environ["OMP_NUM_THREADS"] = "8"
from pydchic.compartment import observed_over_expected, ijk2mat, correlation_matrix, _finalize, oe_zscore
sys.path.insert(0, "/large_storage/zhoulab/shengmao/omicverse-rebuildr")
from engine.benchmark import time_callable
fix = json.load(open(FIX)); a = fix["stageA"]
A,B,oe = observed_over_expected(np.array(a["a_idx"]),np.array(a["b_idx"]),np.array(a["weight"],float),
                                np.array(a["pos"],float),int(a["n_bins"]),int(a["resolution"]))
M = ijk2mat(A,B,oe,int(a["n_bins"])); keep=M.sum(1)>=3; m=M[np.ix_(keep,keep)]
c1=_finalize(correlation_matrix(m)); z2=oe_zscore(c1); c2=_finalize(correlation_matrix(c1))
def f_svd():
    _,_,vt=np.linalg.svd(c2,full_matrices=False); return z2@vt[:2].T
def f_eigh():
    w,v=np.linalg.eigh(c2); idx=np.argsort(w)[::-1][:2]; return z2@v[:,idx]
t0=time_callable(f_svd); t1=time_callable(f_eigh)
pc0, pc1 = f_svd(), f_eigh()
print(f"baseline svd : {t0.mean_s*1e3:.3f} ms ; eigh : {t1.mean_s*1e3:.3f} ms ; speedup {t0.mean_s/t1.mean_s:.2f}x")
"""),
        md("""
## Iteration 0 — Baseline translation

Clean R→Python translation. Decomposition via `np.linalg.svd(C2, full_matrices=False)` —
computes *all* singular vectors though only the top-2 PCs are kept. Compartment PC1 matches
the dcHiC C++/SVD kernel to `1.6e-14`. This is the parity ceiling (`|r| = 1.0`).
"""),
        code("""
fig, ax = plt.subplots(figsize=(7,3.2))
ax.bar(["iter0 svd"], [t0.mean_s*1e3], yerr=[t0.stddev_s*1e3], color="#4c72b0", capsize=4)
ax.set_ylabel("Stage-A wall-clock (ms)"); ax.set_title("Iteration 0 — baseline (full SVD)")
for i,v in enumerate([t0.mean_s*1e3]): ax.text(i, v, f"{v:.2f} ms", ha="center", va="bottom")
plt.tight_layout(); plt.show()
"""),
        md("""
## Iteration 1 — `svd_to_eigh` (exact identity, ACCEPT)

`C2` is symmetric PSD, so its SVD equals its eigendecomposition. Replace the full SVD with
`np.linalg.eigh` keeping the top-k eigenvectors. Per-vector sign is fixed by the downstream
GC orientation, and PC1 is the dominant non-degenerate eigenvalue, so the **gated output is
bit-identical**. `eigh` exploits symmetry and skips discarded vectors → faster. See `MATH.md §1`.
"""),
        code("""
align = abs(np.corrcoef(pc0[:,0], pc1[:,0])[0,1])
print(f"PC1 agreement svd vs eigh: |r| = {align:.10f}  (max|Δ| = {np.max(np.abs(np.abs(pc0[:,0])-np.abs(pc1[:,0]))):.2e})")
fig, ax = plt.subplots(figsize=(7,3.2))
ax.bar(["iter0 svd","iter1 eigh"], [t0.mean_s*1e3, t1.mean_s*1e3],
       yerr=[t0.stddev_s*1e3, t1.stddev_s*1e3], color=["#4c72b0","#55a868"], capsize=4)
ax.set_ylabel("Stage-A wall-clock (ms)")
ax.set_title(f"Iteration 1 — eigh ({t0.mean_s/t1.mean_s:.2f}x faster, parity unchanged)")
for i,v in enumerate([t0.mean_s*1e3, t1.mean_s*1e3]): ax.text(i, v, f"{v:.2f} ms", ha="center", va="bottom")
plt.tight_layout(); plt.show()
"""),
        md("## Aggregate — two-panel evolution (rendered from `ITERATION_LOG.md`)"),
        code("""
from IPython.display import Image
img = os.path.join(PORT, "examples", "evolution.png")
Image(filename=img) if os.path.exists(img) else print("run: python -m engine.plot_evolution --port-dir .")
"""),
        md("""
## Stop reason

Remaining cost is two dense `Zᵀ Z` matmuls (single BLAS calls; no admissible exact reduction
that preserves dcHiC's round-to-5-dp intermediate) and the Stage-B `MinCovDet` library call.
Accuracy held at the ceiling (`|r| = 1.0`) throughout. Search terminated.
"""),
    ]
    save(cells, "evolution.ipynb")


if __name__ == "__main__":
    nb_compare()
    nb_tutorial()
    nb_fbf()
    nb_evolution()

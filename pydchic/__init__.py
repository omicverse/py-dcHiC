"""py-dcHiC — pure-Python port of dcHiC (differential compartment analysis of Hi-C).

Scope: compartment calling (``cis`` + ``select``) and differential analysis (``analyze``).
Out of scope: subcompartment HMM, fithic differential loops, IGV/HTML viz, gene enrichment.

Provenance: re-implemented under reference-driven parity against ay-lab/dcHiC
(Nature Communications 2022) via the omicverse-rebuildr protocol.
"""

from .core import dcHiC
from .compartment import (
    call_compartments,
    compartment_pca,
    observed_over_expected,
    orient_pc,
)
from .differential import (
    differential_compartments,
    normalize_quantiles,
    quantile_normalize_by_condition,
    calcen,
    bh_adjust,
)

__all__ = [
    "dcHiC",
    "call_compartments",
    "compartment_pca",
    "observed_over_expected",
    "orient_pc",
    "differential_compartments",
    "normalize_quantiles",
    "quantile_normalize_by_condition",
    "calcen",
    "bh_adjust",
]
__version__ = "0.1.0"

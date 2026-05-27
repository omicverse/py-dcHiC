"""Class API orchestrating the two dcHiC stages."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .compartment import call_compartments
from .differential import differential_compartments, quantile_normalize_by_condition

__all__ = ["dcHiC"]


class dcHiC:
    """Differential compartment analysis of Hi-C data (compartment calling + differential).

    Mirrors the dcHiC ``cis``/``select`` (compartment calling) and ``analyze``
    (differential) steps. Method-chaining writes results onto the instance.
    """

    def __init__(self, resolution: int):
        self.resolution = int(resolution)
        self.pc1_: dict[str, np.ndarray] = {}
        self.keep_: dict[str, np.ndarray] = {}
        self.differential_: pd.DataFrame | None = None

    # -- Stage A ---------------------------------------------------------------
    def call_compartment(
        self,
        sample: str,
        a_idx,
        b_idx,
        weight,
        pos,
        gcc,
        n_bins: int,
        n_pcs: int = 2,
        count_thr: float = 0.0,
        minexpcc: float = 0.0,
    ) -> "dcHiC":
        """Compute GC-oriented compartment PCs for one sample; store PC1."""
        pc, keep = call_compartments(
            a_idx, b_idx, weight, pos, gcc, n_bins, self.resolution,
            n_pcs=n_pcs, count_thr=count_thr, minexpcc=minexpcc,
        )
        self.pc1_[sample] = pc[:, 0]
        self.keep_[sample] = keep
        return self

    # -- Stage B ---------------------------------------------------------------
    def differential(
        self,
        pc_raw: np.ndarray,
        conditions,
        bed: pd.DataFrame | None = None,
        rzscore: float = 2.0,
        szscore: float = 0.0,
        refine: bool = True,
        rconf: float = 0.90,
        seed: int = 123,
        already_normalized: bool = False,
    ) -> "dcHiC":
        """Differential compartment statistic across conditions.

        ``pc_raw`` is bins x replicates (columns grouped by condition). Set
        ``already_normalized=True`` to skip the per-condition quantile normalisation.
        """
        pc_qnm = pc_raw if already_normalized else quantile_normalize_by_condition(pc_raw, conditions)
        res = differential_compartments(
            pc_qnm, conditions, rzscore=rzscore, szscore=szscore,
            refine=refine, rconf=rconf, seed=seed,
        )
        df = pd.DataFrame({"pval": res["pval"], "padj": res["padj"], "sample_maha": res["sample_maha"]})
        if bed is not None:
            df = pd.concat([bed.reset_index(drop=True)[["chr", "start", "end"]], df], axis=1)
        self.differential_ = df
        return self

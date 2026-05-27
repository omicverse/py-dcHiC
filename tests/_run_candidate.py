"""Candidate runner — invoked under conda env $PYTHON_TEST_ENV.

Usage: python _run_candidate.py <fixture.json> <output.json>

Runs the pydchic port on the canonical fixture; dumps outputs keyed to match
manifest.yaml::outputs[*].location_candidate.
"""

import json
import sys

import numpy as np

from pydchic import call_compartments
from pydchic.differential import differential_compartments, quantile_normalize_by_condition


def main():
    if len(sys.argv) != 3:
        print("Usage: python _run_candidate.py <fixture.json> <output.json>")
        sys.exit(2)
    fixture_path, output_path = sys.argv[1], sys.argv[2]

    with open(fixture_path) as f:
        fix = json.load(f)
    np.random.seed(123)

    # ---- Stage A: compartment calling --------------------------------------
    a = fix["stageA"]
    pc, keep = call_compartments(
        np.array(a["a_idx"]), np.array(a["b_idx"]), np.array(a["weight"], dtype=float),
        np.array(a["pos"], dtype=float), np.array(a["gcc"], dtype=float),
        int(a["n_bins"]), int(a["resolution"]), n_pcs=2, count_thr=0.0, minexpcc=0.0,
    )
    compartment_pc1 = pc[:, 0]

    # ---- Stage B: differential ---------------------------------------------
    b = fix["stageB"]
    pc_raw = np.asarray(b["pc_raw"], dtype=float)        # bins x replicates
    conditions = b["conditions"]
    pc_qnm = quantile_normalize_by_condition(pc_raw, conditions)
    res = differential_compartments(
        pc_qnm, conditions, rzscore=2.0, szscore=0.0, refine=True, rconf=0.90, seed=123,
    )

    out = {
        "compartment_pc1": np.asarray(compartment_pc1, dtype=float).tolist(),
        "differential_padj": np.asarray(res["padj"], dtype=float).tolist(),
    }
    with open(output_path, "w") as f:
        json.dump(out, f)
    print(f"[cand] wrote: {output_path}")


if __name__ == "__main__":
    main()

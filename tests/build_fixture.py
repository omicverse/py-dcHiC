"""Build the canonical parity fixture -> data/fixture_dchic.json.

Stage A (compartment calling): a deterministic synthetic single-chromosome Hi-C contact
matrix with clear A/B compartment (checkerboard) structure on top of a power-law
distance decay, plus a GC-content track correlated with the compartments.

Stage B (differential): real dcHiC demo PC1 values for chr19 across the ESC/NPC/CN
mouse dataset (4 + 3 + 5 replicates), as shipped in dcHiC_demo.zip.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PORT = HERE.parent
# Vendored demo PC files (chr19, ESC/NPC/CN) so the fixture is reproducible on a clone.
VENDORED_DEMO = HERE / "demo_data"

RESOLUTION = 100_000
N_BINS = 200
SEED = 123


def build_stage_a():
    rng = np.random.default_rng(SEED)
    # compartment label per bin: 4 blocks of alternating A(+1)/B(-1)
    block = N_BINS // 4
    comp = np.concatenate([
        np.full(block, 1.0), np.full(block, -1.0),
        np.full(block, 1.0), np.full(N_BINS - 3 * block, -1.0),
    ])
    # GC track correlated with compartment (A-rich = high GC)
    gcc = 0.45 + 0.10 * comp + rng.normal(0, 0.01, N_BINS)

    a_idx, b_idx, weight = [], [], []
    for i in range(N_BINS):
        for j in range(i, N_BINS):
            d = abs(i - j)
            decay = (d + 1.0) ** (-0.9)
            compartment = 1.0 + 0.6 * comp[i] * comp[j]
            noise = rng.lognormal(0.0, 0.15)
            count = 800.0 * decay * compartment * noise
            c = int(round(count))
            if c > 0:
                a_idx.append(i + 1)        # 1-based bin ids
                b_idx.append(j + 1)
                weight.append(c)

    pos = (np.arange(N_BINS) * RESOLUTION).tolist()  # bin start coordinates (bp)
    return {
        "a_idx": a_idx,
        "b_idx": b_idx,
        "weight": weight,
        "pos": pos,
        "gcc": gcc.tolist(),
        "n_bins": N_BINS,
        "resolution": RESOLUTION,
    }


def build_stage_b():
    samples = [
        ("ES_1", "ES"), ("ES_2", "ES"), ("ES_3", "ES"), ("ES_4", "ES"),
        ("NPC_2", "NPC"), ("NPC_3", "NPC"), ("NPC_4", "NPC"),
        ("CN_1", "CN"), ("CN_2", "CN"), ("CN_3", "CN"), ("CN_4", "CN"), ("CN_5", "CN"),
    ]
    cols, conditions, coords = [], [], None
    import pandas as pd
    for sample, cond in samples:
        f = VENDORED_DEMO / f"{sample}_chr19.pc.txt"
        df = pd.read_csv(f, sep="\t", header=0)
        if coords is None:
            coords = df[["chr", "start", "end"]].copy()
        cols.append(df["PC1"].to_numpy())
        conditions.append(cond)

    pc_raw = np.column_stack(cols)  # bins x replicates (grouped by condition order)
    # orient each replicate's PC1 to positively correlate with the first column
    ref = pc_raw[:, 0]
    for k in range(pc_raw.shape[1]):
        if np.corrcoef(pc_raw[:, k], ref)[0, 1] < 0:
            pc_raw[:, k] = -pc_raw[:, k]

    return {
        "pc_raw": pc_raw.tolist(),         # rows = bins, cols = replicates
        "conditions": conditions,
        "chr": coords["chr"].tolist(),
        "start": coords["start"].tolist(),
        "end": coords["end"].tolist(),
    }


def main():
    out = PORT / "data" / "fixture_dchic.json"
    fixture = {"stageA": build_stage_a(), "stageB": build_stage_b()}
    out.write_text(json.dumps(fixture))
    a, b = fixture["stageA"], fixture["stageB"]
    print(f"[fixture] Stage A: {a['n_bins']} bins, {len(a['weight'])} nonzero contacts")
    print(f"[fixture] Stage B: {len(b['start'])} bins x {len(b['conditions'])} replicates "
          f"({sorted(set(b['conditions']))})")
    print(f"[fixture] wrote {out} ({out.stat().st_size/1e6:.2f} MB)")


if __name__ == "__main__":
    sys.exit(main())

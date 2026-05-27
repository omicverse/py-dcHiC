"""I/O helpers for dcHiC file formats (sparse matrix, bed, pc.txt)."""

from __future__ import annotations

import pandas as pd

__all__ = ["read_contact_matrix", "read_bed", "read_pc_txt"]


def read_contact_matrix(path: str) -> pd.DataFrame:
    """Read a 3-column sparse contact matrix (bin_i, bin_j, count); 1-based indices."""
    df = pd.read_csv(path, sep="\t", header=None, names=["A", "B", "Weight"])
    return df


def read_bed(path: str) -> pd.DataFrame:
    """Read a dcHiC bed file (4 or 5 columns)."""
    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] == 4:
        df.columns = ["chr", "start", "end", "index"]
    elif df.shape[1] == 5:
        df.columns = ["chr", "start", "end", "index", "blacklist"]
    else:
        raise ValueError(f"bed file must have 4 or 5 columns, got {df.shape[1]}")
    return df


def read_pc_txt(path: str) -> pd.DataFrame:
    """Read a dcHiC ``*.pc.txt`` file (chr, start, end, index, PC1, PC2, ...)."""
    return pd.read_csv(path, sep="\t", header=0)

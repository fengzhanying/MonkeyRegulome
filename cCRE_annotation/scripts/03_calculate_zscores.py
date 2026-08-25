#!/usr/bin/env python3
"""Compute mark-wise maximum z-scores across macaque experiments."""

from __future__ import annotations

import os
import sys
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd


COL_NAMES = ["name", "size", "covered", "sum", "mean0", "mean"]


def require_env(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"ERROR: set {name} in the config file.")
    return value


def main() -> None:
    workdir = Path(require_env("CCRE_WORKDIR"))
    sigdir = workdir / "02.signals"
    zsdir = workdir / "03.zscores"
    zsdir.mkdir(parents=True, exist_ok=True)

    mark_patterns = {
        "ATAC": str(sigdir / "ATAC-*.tab"),
        "H3K4me3": str(sigdir / "CUT-Tag-*-H3K4me3.tab"),
        "H3K27ac": str(sigdir / "CUT-Tag-*-H3K27ac.tab"),
        "CTCF": str(sigdir / "CUT-Tag-*-CTCF.tab"),
    }

    for mark, pattern in mark_patterns.items():
        files = sorted(glob(pattern))
        if not files:
            print(f"[WARN] {mark}: no files matched pattern '{pattern}'")
            continue

        print(f"[{mark}] {len(files)} experiment(s) found")
        per_exp_zscores = {}

        for fpath in files:
            sname = Path(fpath).name.replace(".tab", "")
            try:
                df = pd.read_csv(
                    fpath,
                    sep="\t",
                    header=None,
                    names=COL_NAMES,
                    usecols=["name", "mean0"],
                )
            except Exception as exc:
                print(f"  [SKIP] {sname}: read error -- {exc}")
                continue

            log_signal = np.log1p(df["mean0"].clip(lower=0).to_numpy())
            mean = log_signal.mean()
            std = log_signal.std(ddof=1)
            z = (log_signal - mean) / std if std > 0 else np.zeros_like(log_signal)

            per_exp_zscores[sname] = pd.Series(
                z, index=df["name"].to_numpy(), dtype=np.float32
            )

        if not per_exp_zscores:
            print(f"  [SKIP] {mark}: no valid experiments")
            continue

        zmat = pd.DataFrame(per_exp_zscores)
        print(f"  Z-score matrix: {zmat.shape[0]} rDHSs x {zmat.shape[1]} experiments")
        zmat["max_zscore"] = zmat.max(axis=1)

        zmat.to_csv(
            zsdir / f"{mark}_all_samples_zscore.txt.gz",
            sep="\t",
            compression="gzip",
            index=True,
        )
        (
            zmat[["max_zscore"]]
            .reset_index()
            .rename(columns={"index": "name"})
            .to_csv(zsdir / f"{mark}_max_zscore.txt", sep="\t", index=False)
        )

        print(f"  Saved: {zsdir / f'{mark}_max_zscore.txt'}")

    print("Z-score computation complete.")


if __name__ == "__main__":
    main()

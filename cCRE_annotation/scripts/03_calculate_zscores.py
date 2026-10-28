#!/usr/bin/env python3
"""Compute per-experiment log1p(mean0) z-scores and mark-wise maxima."""
from pathlib import Path
import numpy as np
import pandas as pd
from common import MARKS, cli, parser, read_anchors, read_manifest, read_signal, record, write_tsv

def standardize(values):
    values = np.asarray(values, dtype=np.float64)
    if len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("Z-scores require at least two finite values")
    # Negative input signals were clipped in the original analysis.
    logged = np.log1p(np.clip(values, 0, None))
    std = logged.std(ddof=1)
    z = (logged - logged.mean()) / std if std > 0 else np.zeros_like(logged)
    return z.astype(np.float32)

def main():
    p = parser(__doc__)
    p.add_argument("--manifest", type=Path, help="Default: WORKDIR/logs/experiments.tsv")
    p.add_argument("--save-matrix", action="store_true", help="Also write the optional per-experiment matrix")
    args = p.parse_args()
    manifest_path = args.manifest or args.workdir / "logs/experiments.tsv"
    manifest = read_manifest(manifest_path)
    anchor_path = args.workdir / "01.masterlist/monkey_rDHS.bed"
    anchors = read_anchors(anchor_path)
    sigdir, zsdir = args.workdir / "02.signals", args.workdir / "03.zscores"
    expected = {sample + ".tab" for sample in manifest["sample"]}
    if {path.name for path in sigdir.glob("*.tab")} != expected:
        raise ValueError("Signal files do not exactly match the experiment manifest; use a fresh workdir")
    inputs = [anchor_path, manifest_path]
    outputs = []
    for mark in MARKS:
        maximum = np.full(len(anchors), -np.inf, dtype=np.float32)
        matrix = {}
        samples = manifest.loc[manifest.mark == mark, "sample"]
        for sample in samples:
            path = sigdir / f"{sample}.tab"
            frame = read_signal(path, anchors.id)
            if not np.array_equal(frame["size"].to_numpy(), (anchors.end - anchors.start).to_numpy()):
                raise ValueError(f"{path}: signal region sizes disagree with anchors")
            z = standardize(frame.mean0)
            maximum = np.maximum(maximum, z)
            if args.save_matrix:
                matrix[sample] = z
            inputs.append(path)
        target = zsdir / f"{mark}_max_zscore.txt"
        write_tsv(pd.DataFrame({"name": anchors.id, "max_zscore": maximum}), target)
        outputs.append(target)
        if args.save_matrix:
            frame = pd.DataFrame(matrix, index=anchors.id)
            frame.index.name = "name"
            frame["max_zscore"] = maximum
            target = zsdir / f"{mark}_all_samples_zscore.txt.gz"
            write_tsv(frame.reset_index(), target, compression={"method": "gzip", "mtime": 0})
            outputs.append(target)
        print(f"{mark}: {len(samples)} experiments, {len(anchors)} anchors")
    record(args.workdir, "03_zscores", inputs, outputs,
           {"transform": "log1p(max(mean0, 0))", "ddof": 1, "per_experiment_dtype": "float32"})

if __name__ == "__main__":
    cli(main)

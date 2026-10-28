#!/usr/bin/env python3
"""Validate production inputs and freeze the selected experiments in a fresh workdir."""
import json
import os
import shutil
from pathlib import Path
import pandas as pd
from common import MARKS, atomic_path, cli, parser, read_manifest, record, write_tsv

def main():
    args = parser(__doc__).parse_args()
    if os.environ.get("CCRE_DISTANCE_MODE", "legacy") not in {"legacy", "corrected"}:
        raise ValueError("CCRE_DISTANCE_MODE must be legacy or corrected")
    workdir = args.workdir.resolve()
    if workdir.exists() and any(workdir.iterdir()):
        raise ValueError("Use a new or empty workdir; automatic reuse of previous intermediates is disabled")
    selected = os.environ.get("CCRE_EXPERIMENTS")
    if selected:
        manifest_path = Path(selected).resolve()
        frame = read_manifest(manifest_path)
        frame["bigwig"] = [str((manifest_path.parent / p).resolve()) for p in frame.bigwig]
    else:
        rows = []
        for mark in MARKS:
            var = "CCRE_ATAC_BW_DIR" if mark == "ATAC" else "CCRE_CUTTAG_BW_DIR"
            value = os.environ.get(var)
            if not value or not Path(value).is_dir():
                raise ValueError(f"Set {var} to an existing directory or supply CCRE_EXPERIMENTS")
            pattern = "ATAC-*.fc.signal.bw" if mark == "ATAC" else f"CUT-Tag-*-{mark}.fc.signal.bw"
            for path in sorted(Path(value).resolve().glob(pattern)):
                rows.append((path.name.removesuffix(".fc.signal.bw"), mark, str(path)))
        frame = pd.DataFrame(rows, columns=["sample", "mark", "bigwig"])
    for value in frame.bigwig:
        if not Path(value).is_file():
            raise ValueError(f"Missing bigWig: {value}")
    sources = {}
    for var in ["CCRE_MASTERLIST", "CCRE_GTF"]:
        value = os.environ.get(var)
        if not value or not Path(value).is_file():
            raise ValueError(f"Set {var} to an existing input file")
        sources[var] = str(Path(value).resolve())
    for executable in ["bedtools", "bigWigAverageOverBed"]:
        if not shutil.which(executable):
            raise ValueError(f"Missing executable: {executable}")
    target = workdir / "logs/experiments.tsv"
    if frame.empty or set(frame.mark) != set(MARKS) or frame["sample"].duplicated().any():
        raise ValueError("Provide all four marks and unique sample names")
    if not frame["sample"].str.fullmatch(r"[A-Za-z0-9_.-]+").all():
        raise ValueError("Invalid sample filename")
    write_tsv(frame, target)
    read_manifest(target)
    with atomic_path(workdir / "logs/run_config.json") as tmp:
        tmp.write_text(json.dumps({**sources, "workdir": str(workdir),
            "distance_mode": os.environ.get("CCRE_DISTANCE_MODE", "legacy")}, indent=2) + "\n")
    record(workdir, "prepare_run", [Path(p) for p in sources.values()], [target])
    print(f"Validated {len(frame)} experiments in {workdir}")

if __name__ == "__main__":
    cli(main)

#!/usr/bin/env python3
"""Extract one manifest-selected bigWig to an atomically written, validated signal table."""
import os
import shutil
import subprocess
from common import atomic_path, cli, parser, read_anchors, read_manifest, read_signal, record

def main():
    p = parser(__doc__)
    p.add_argument("--index", type=int, default=os.environ.get("SLURM_ARRAY_TASK_ID"),
                   required=not os.environ.get("SLURM_ARRAY_TASK_ID"))
    args = p.parse_args()
    manifest_path = args.workdir / "logs/experiments.tsv"
    manifest = read_manifest(manifest_path)
    if not 1 <= args.index <= len(manifest):
        raise ValueError("Experiment index is outside the manifest (one-based)")
    row = manifest.iloc[args.index - 1]
    anchor_path = args.workdir / "01.masterlist/monkey_rDHS.bed"
    anchors = read_anchors(anchor_path)
    target = args.workdir / "02.signals" / (row["sample"] + ".tab")
    executable = shutil.which("bigWigAverageOverBed")
    if not executable:
        raise ValueError("bigWigAverageOverBed is not on PATH")
    with atomic_path(target) as tmp:
        subprocess.run([executable, row.bigwig, str(anchor_path), str(tmp)], check=True)
        read_signal(tmp, anchors.id)
    record(args.workdir, "02_signal_" + row["sample"],
           [anchor_path, manifest_path, row.bigwig, executable], [target],
           {"sample": row["sample"], "mark": row.mark})
    print(f"Extracted {row['sample']}")

if __name__ == "__main__":
    cli(main)

#!/usr/bin/env python3
"""Normalize a supplied merged ATAC BED; this step does not merge raw peaks."""
import os
from pathlib import Path
import pandas as pd
from common import cli, parser, read_anchors, record, write_tsv

def main():
    p = parser(__doc__)
    p.add_argument("--masterlist", type=Path, default=os.environ.get("CCRE_MASTERLIST"),
                   required=not os.environ.get("CCRE_MASTERLIST"))
    args = p.parse_args()
    rows = []
    with args.masterlist.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"BED line {line_number}: fewer than three fields")
            chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"BED line {line_number}: invalid interval")
            # Historical rule, not a full canonical-contig whitelist.
            if "_" in chrom or chrom in {"chrM", "chrMT"}:
                continue
            rows.append((chrom, start, end, f"{chrom}_{start}_{end}"))
    frame = pd.DataFrame(rows, columns=["chr", "start", "end", "id"])
    if frame.empty or frame.id.duplicated().any():
        raise ValueError("No retained anchors or duplicate intervals in master list")
    frame = frame.sort_values(["chr", "start", "end"], kind="stable")
    target = args.workdir / "01.masterlist/monkey_rDHS.bed"
    write_tsv(frame, target, header=False)
    read_anchors(target)
    record(args.workdir, "00_masterlist", [args.masterlist], [target])
    print(f"Retained {len(frame)} ATAC anchors: {target}")

if __name__ == "__main__":
    cli(main)

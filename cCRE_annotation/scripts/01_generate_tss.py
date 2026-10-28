#!/usr/bin/env python3
"""Generate deduplicated transcript TSS positions from a plain or gzipped GTF."""
import gzip
import os
import re
from pathlib import Path
import pandas as pd
from common import cli, parser, read_anchors, record, write_tsv

def main():
    p = parser(__doc__)
    p.add_argument("--gtf", type=Path, default=os.environ.get("CCRE_GTF"),
                   required=not os.environ.get("CCRE_GTF"))
    args = p.parse_args()
    anchors_path = args.workdir / "01.masterlist/monkey_rDHS.bed"
    anchors = read_anchors(anchors_path)
    has_prefix = anchors.chr.iloc[0].startswith("chr")
    rows, seen = [], set()
    opener = gzip.open if args.gtf.suffix == ".gz" else open
    with opener(args.gtf, "rt", encoding="utf-8") as stream:
        for n, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {n}: expected nine fields")
            if fields[2] != "transcript":
                continue
            chrom, strand = fields[0], fields[6]
            if strand not in {"+", "-"}:
                raise ValueError(f"GTF line {n}: invalid transcript strand")
            if has_prefix and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            elif not has_prefix and chrom.startswith("chr"):
                chrom = chrom[3:]
            if "_" in chrom or chrom in {"chrM", "chrMT", "M", "MT", "chrEBV"}:
                continue
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"GTF line {n}: invalid coordinates")
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not match:
                match = re.search(r'gene_id "([^"]+)"', fields[8])
            gene = match.group(1) if match else "NA"
            pos = (start if strand == "+" else end) - 1
            key = chrom, pos, strand
            if key not in seen:
                seen.add(key)
                rows.append((chrom, pos, pos + 1, gene, 0, strand))
    if not rows or not set(anchors.chr).intersection(r[0] for r in rows):
        raise ValueError("No TSS on any anchor chromosome; check assembly/naming")
    rows.sort(key=lambda row: (row[0], row[1]))
    target = args.workdir / "01.masterlist/monkey_TSS.bed"
    write_tsv(pd.DataFrame(rows), target, header=False)
    record(args.workdir, "01_tss", [args.gtf, anchors_path], [target])
    print(f"Wrote {len(rows)} unique TSS positions: {target}")

if __name__ == "__main__":
    cli(main)

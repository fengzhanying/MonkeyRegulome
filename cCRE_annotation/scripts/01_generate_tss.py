#!/usr/bin/env python3
"""Generate a non-redundant macaque TSS BED file from a GTF annotation."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path


def require_env(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"ERROR: set {name} in the config file.")
    return value


def main() -> None:
    workdir = Path(require_env("CCRE_WORKDIR"))
    gtf = Path(require_env("CCRE_GTF"))
    outdir = workdir / "01.masterlist"
    rdhs = outdir / "monkey_rDHS.bed"
    tss_out = outdir / "monkey_TSS.bed"

    if not gtf.exists():
        sys.exit(f"ERROR: GTF file not found: {gtf}")
    if not rdhs.exists():
        sys.exit(f"ERROR: rDHS BED file not found: {rdhs}")

    with rdhs.open() as handle:
        first_chrom = handle.readline().split("\t", 1)[0]
    has_chr_prefix = first_chrom.startswith("chr")

    tss_seen = set()
    records = []

    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            chrom = fields[0]
            strand = fields[6]
            attr = fields[8]

            if has_chr_prefix and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            elif not has_chr_prefix and chrom.startswith("chr"):
                chrom = chrom[3:]

            if "_" in chrom or chrom in {"chrM", "chrMT", "M", "MT", "chrEBV"}:
                continue

            match = re.search(r'gene_name "([^"]+)"', attr)
            gene_name = match.group(1) if match else None
            if not gene_name:
                match = re.search(r'gene_id "([^"]+)"', attr)
                gene_name = match.group(1) if match else "NA"

            start_gtf = int(fields[3])
            end_gtf = int(fields[4])
            tss_start = start_gtf - 1 if strand == "+" else end_gtf - 1
            tss_end = tss_start + 1

            key = (chrom, tss_start, strand)
            if key in tss_seen:
                continue
            tss_seen.add(key)
            records.append((chrom, tss_start, tss_end, gene_name, 0, strand))

    records.sort(key=lambda record: (record[0], record[1]))
    with tss_out.open("w", newline="\n") as out:
        for record in records:
            out.write("\t".join(map(str, record)) + "\n")

    print(f"Total unique TSS positions: {len(records)}")
    print(f"Output: {tss_out}")


if __name__ == "__main__":
    main()

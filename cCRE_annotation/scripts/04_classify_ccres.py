#!/usr/bin/env python3
"""Classify global accessible anchors, with explicit historical distance compatibility."""
import os
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from common import CLASSES, MARKS, atomic_path, cli, parser, read_anchors, record, write_tsv

Z_THRESH = 1.64
TSS_CORE_BP = 200
TSS_PROX_BP = 2000

def classify(row):
    k4 = row["z_h3k4me3"] > Z_THRESH
    k27 = row["z_h3k27ac"] > Z_THRESH
    ctcf = row["z_ctcf"] > Z_THRESH
    if row["tss_core"] and k4:
        return "PLS"
    if row["tss_proximal"] and k27:
        return "pELS"
    if k27:
        return "dELS"
    if k4:
        return "CA-H3K4me3"
    if ctcf:
        return "CA-CTCF"
    return "CA"

def parse_nearest(path, anchors, mode):
    near = pd.read_csv(path, sep="\t", header=None, keep_default_na=False)
    if near.shape[1] != 11:
        raise ValueError("Nearest-TSS input must be BED4 + BED6 + signed distance (11 columns)")
    near.columns = ["chr", "start", "end", "id", "tss_chr", "tss_start",
                    "tss_end", "gene", "score", "strand", "dist"]
    if near.id.duplicated().any() or set(near.id) != set(anchors.id):
        raise ValueError("Nearest-TSS records must match anchor IDs exactly, one per anchor")
    near = near.set_index("id").loc[anchors.id].reset_index()
    for col in ["start", "end", "tss_start", "tss_end", "dist"]:
        near[col] = pd.to_numeric(near[col], errors="raise")
        if not np.isfinite(near[col]).all() or (near[col] % 1 != 0).any():
            raise ValueError("Nearest-TSS coordinates/distances must be finite integers")
    if not near[["chr", "start", "end", "id"]].equals(anchors.reset_index(drop=True)):
        raise ValueError("Nearest-TSS coordinates disagree with anchor BED")
    missing = near.tss_chr.isin(["none", "."]) | (near.tss_start < 0) | (near.tss_end < 0)
    valid = ~missing
    if (near.loc[valid, "tss_chr"] != near.loc[valid, "chr"]).any() or (
        near.loc[valid, "tss_end"] - near.loc[valid, "tss_start"] != 1
    ).any():
        raise ValueError("Nearest TSS must be a one-base interval on the anchor chromosome")
    if mode == "legacy":
        # Deliberate historical compatibility: an upstream distance of -1 was
        # treated as missing. Corrected mode uses only the missing-B fields.
        missing = missing | (near.dist == -1)
    near["center_dist"] = abs((near.start + near.end - near.tss_start - near.tss_end) / 2)
    near.loc[missing, "center_dist"] = 1e12
    return near

def main():
    p = parser(__doc__)
    p.add_argument("--distance-mode", choices=["legacy", "corrected"],
                   default=os.environ.get("CCRE_DISTANCE_MODE", "legacy"))
    p.add_argument("--nearest-tss", type=Path,
                   help="Explicit precomputed bedtools result (used by the core demo)")
    args = p.parse_args()
    rdhs_path = args.workdir / "01.masterlist/monkey_rDHS.bed"
    tss_path = args.workdir / "01.masterlist/monkey_TSS.bed"
    anchors = read_anchors(rdhs_path)
    near_path = args.nearest_tss or args.workdir / "01.masterlist/monkey_rDHS_nearTSS.bed"
    external = "precomputed; bedtools not invoked"
    if args.nearest_tss is None:
        external = subprocess.check_output(["bedtools", "--version"], text=True).strip()
        with atomic_path(near_path) as temporary, temporary.open("w", encoding="utf-8", newline="\n") as stream:
            subprocess.run(["bedtools", "closest", "-a", str(rdhs_path), "-b", str(tss_path),
                            "-D", "ref", "-t", "first"], stdout=stream, check=True)
    near = parse_nearest(near_path, anchors, args.distance_mode)
    frame = near[["chr", "start", "end", "id", "center_dist", "gene"]].copy()
    frame["tss_core"] = frame.center_dist <= TSS_CORE_BP
    frame["tss_proximal"] = frame.center_dist <= TSS_PROX_BP
    inputs = [rdhs_path, tss_path, near_path]
    for mark in MARKS:
        path = args.workdir / "03.zscores" / f"{mark}_max_zscore.txt"
        scores = pd.read_csv(path, sep="\t")
        if list(scores.columns) != ["name", "max_zscore"] or scores.name.duplicated().any() or set(scores.name) != set(anchors.id):
            raise ValueError(f"{path}: expected one maximum z-score per anchor")
        scores = scores.set_index("name").loc[anchors.id, "max_zscore"]
        if not np.isfinite(scores).all():
            raise ValueError(f"{path}: non-finite maximum z-scores")
        frame["z_" + mark.lower()] = scores.to_numpy()
        inputs.append(path)
    frame["cCRE_class"] = frame.apply(classify, axis=1)
    # Sorted coordinates simplify browser and interval-tool reuse; labels are unchanged.
    frame = frame.sort_values(["chr", "start", "end"], kind="stable")
    frame = frame[["id"] + [c for c in frame.columns if c != "id"]]
    out = args.workdir / "04.classification"
    annotated, bed, stats = (out / name for name in
        ["monkey_cCREs_annotated.txt", "monkey_cCREs.bed", "cCRE_class_stats.txt"])
    write_tsv(frame, annotated)
    write_tsv(frame[["chr", "start", "end", "id", "cCRE_class"]], bed, header=False)
    counts = frame.cCRE_class.value_counts().reindex(CLASSES, fill_value=0)
    write_tsv(counts.rename_axis("cCRE_class").reset_index(name="count"), stats)
    record(args.workdir, "04_classification", inputs, [annotated, bed, stats],
           {"distance_mode": args.distance_mode, "threshold_strict_gt": Z_THRESH,
            "core_bp_inclusive": TSS_CORE_BP, "proximal_bp_inclusive": TSS_PROX_BP,
            "bedtools": external, "ATAC_filter": False})
    print(counts.to_string())
    print(f"Total: {len(frame)}; distance mode: {args.distance_mode}")

if __name__ == "__main__":
    cli(main)

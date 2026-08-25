#!/usr/bin/env python3
"""Classify macaque accessible anchors into global cCRE categories."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


Z_THRESH = 1.64
TSS_CORE_BP = 200
TSS_PROX_BP = 2000


def require_env(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"ERROR: set {name} in the config file.")
    return value


def main() -> None:
    workdir = Path(require_env("CCRE_WORKDIR"))
    outdir = workdir / "04.classification"
    outdir.mkdir(parents=True, exist_ok=True)

    rdhs_bed = workdir / "01.masterlist" / "monkey_rDHS.bed"
    tss_bed = workdir / "01.masterlist" / "monkey_TSS.bed"
    near_tss = workdir / "01.masterlist" / "monkey_rDHS_nearTSS.bed"
    zsdir = workdir / "03.zscores"

    if not rdhs_bed.exists():
        sys.exit(f"ERROR: rDHS BED file not found: {rdhs_bed}")
    if not tss_bed.exists():
        sys.exit(f"ERROR: TSS BED file not found: {tss_bed}")

    if not near_tss.exists():
        print("Computing distance from each anchor to the nearest TSS...")
        with near_tss.open("w", newline="\n") as out:
            ret = subprocess.run(
                [
                    "bedtools",
                    "closest",
                    "-a",
                    str(rdhs_bed),
                    "-b",
                    str(tss_bed),
                    "-D",
                    "ref",
                    "-t",
                    "first",
                ],
                stdout=out,
                check=False,
            )
        if ret.returncode != 0:
            sys.exit(
                "ERROR: bedtools closest failed. Check genome assembly and "
                "chromosome naming consistency."
            )
    else:
        print(f"Reusing existing nearest-TSS file: {near_tss}")

    near = pd.read_csv(
        near_tss,
        sep="\t",
        header=None,
        names=[
            "chr",
            "start",
            "end",
            "id",
            "tss_chr",
            "tss_start",
            "tss_end",
            "gene",
            "score",
            "strand",
            "dist",
        ],
    )

    for col in ["start", "end", "tss_start", "tss_end", "dist"]:
        near[col] = pd.to_numeric(near[col], errors="coerce")

    no_tss = near["tss_start"].isna() | near["tss_end"].isna() | (near["dist"] == -1)
    anchor_center = (near["start"] + near["end"]) / 2.0
    tss_center = (near["tss_start"] + near["tss_end"]) / 2.0
    near["center_dist"] = (anchor_center - tss_center).abs()
    near.loc[no_tss, "center_dist"] = float(1e12)

    near = near.sort_values("center_dist").drop_duplicates("id", keep="first")
    rdhs = near[["chr", "start", "end", "id", "center_dist", "gene"]].copy()
    rdhs["tss_core"] = rdhs["center_dist"] <= TSS_CORE_BP
    rdhs["tss_proximal"] = rdhs["center_dist"] <= TSS_PROX_BP
    rdhs = rdhs.set_index("id")

    def load_max_z(mark: str) -> pd.Series:
        path = zsdir / f"{mark}_max_zscore.txt"
        if not path.exists():
            print(f"  [WARN] Missing Z-score file for {mark}; defaulting to 0.")
            return pd.Series(0.0, index=rdhs.index, name=f"z_{mark.lower()}")
        df = pd.read_csv(path, sep="\t").set_index("name")
        return df["max_zscore"].rename(f"z_{mark.lower()}")

    print("Loading max Z-scores...")
    rdhs["z_atac"] = load_max_z("ATAC")
    rdhs["z_h3k4me3"] = load_max_z("H3K4me3")
    rdhs["z_h3k27ac"] = load_max_z("H3K27ac")
    rdhs["z_ctcf"] = load_max_z("CTCF")
    rdhs = rdhs.fillna(0)

    def classify(row: pd.Series) -> str:
        k4me3 = row["z_h3k4me3"] > Z_THRESH
        k27ac = row["z_h3k27ac"] > Z_THRESH
        ctcf = row["z_ctcf"] > Z_THRESH
        core = bool(row["tss_core"])
        prox = bool(row["tss_proximal"])

        if core and k4me3:
            return "PLS"
        if prox and k27ac:
            return "pELS"
        if (not prox) and k27ac:
            return "dELS"
        if k4me3:
            return "CA-H3K4me3"
        if ctcf:
            return "CA-CTCF"
        return "CA"

    print("Classifying cCREs...")
    rdhs["cCRE_class"] = rdhs.apply(classify, axis=1)
    rdhs_out = rdhs.reset_index()

    annotated_out = outdir / "monkey_cCREs_annotated.txt"
    bed_out = outdir / "monkey_cCREs.bed"
    stats_out = outdir / "cCRE_class_stats.txt"

    rdhs_out.to_csv(annotated_out, sep="\t", index=False)
    rdhs_out[["chr", "start", "end", "id", "cCRE_class"]].to_csv(
        bed_out, sep="\t", header=False, index=False
    )
    rdhs_out["cCRE_class"].value_counts().to_csv(
        stats_out, sep="\t", header=["count"]
    )

    print("\n=== cCRE classification summary ===")
    print(rdhs_out["cCRE_class"].value_counts().to_string())
    print(f"\nTotal ATAC anchors retained: {len(rdhs_out)}")
    print(f"Output BED: {bed_out}")


if __name__ == "__main__":
    main()

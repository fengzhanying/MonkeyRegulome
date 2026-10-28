#!/usr/bin/env python3
"""Check released global coordinates, labels and counts; optionally refresh their manifest."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import CLASSES, ROOT, atomic_path, cli, sha256, write_tsv

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, default=ROOT / "results/global")
    p.add_argument("--write-manifest", action="store_true")
    args = p.parse_args()
    folder = args.results
    bed = pd.read_csv(folder / "monkey_cCREs.bed", sep="\t", header=None,
                      names=["chr", "start", "end", "id", "cCRE_class"])
    full = pd.read_csv(folder / "monkey_cCREs_annotated.txt", sep="\t")
    if bed.id.duplicated().any() or full.id.duplicated().any() or set(bed.id) != set(full.id):
        raise ValueError("BED and full annotations have missing or duplicated IDs")
    a = bed.set_index("id").sort_index()
    b = full.set_index("id").sort_index()
    pd.testing.assert_frame_equal(a, b[a.columns], check_dtype=False)
    if (a.start < 0).any() or (a.end <= a.start).any() or not set(a.cCRE_class).issubset(CLASSES):
        raise ValueError("Invalid coordinates or class labels")
    if not (bed.id == bed.chr + "_" + bed.start.astype(str) + "_" + bed.end.astype(str)).all():
        raise ValueError("IDs disagree with coordinates")
    z = b[["z_atac", "z_h3k4me3", "z_h3k27ac", "z_ctcf"]]
    if not np.isfinite(z.to_numpy()).all():
        raise ValueError("Non-finite annotated signals")
    for column in ["tss_core", "tss_proximal"]:
        if not b[column].isin([True, False]).all():
            raise ValueError(f"Invalid boolean column {column}")
    if not (b.tss_core == (b.center_dist <= 200)).all() or not (b.tss_proximal == (b.center_dist <= 2000)).all():
        raise ValueError("Stored distance flags disagree with distances")
    expected = np.select([b.tss_core & (b.z_h3k4me3 > 1.64),
                          b.tss_proximal & (b.z_h3k27ac > 1.64),
                          b.z_h3k27ac > 1.64, b.z_h3k4me3 > 1.64, b.z_ctcf > 1.64],
                         ["PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF"], default="CA")
    if not np.array_equal(expected, b.cCRE_class.to_numpy()):
        raise ValueError("Stored classifications disagree with the documented rule")
    counts = bed.cCRE_class.value_counts().reindex(CLASSES, fill_value=0)
    stats = pd.read_csv(folder / "cCRE_class_stats.txt", sep="\t", index_col=0).iloc[:, 0]
    if not np.array_equal(counts.to_numpy(), stats.reindex(CLASSES).to_numpy()):
        raise ValueError("Class-count file disagrees with coordinate-level data")
    pretty = pd.read_csv(folder / "Monkey_cCRE_Composition_Stats.csv")
    mapping = {"Promoter": "PLS", "Proximal enhancer": "pELS", "Distal enhancer": "dELS"}
    pretty["Class"] = pretty.Class.replace(mapping)
    pretty = pretty.set_index("Class").reindex(CLASSES)
    np.testing.assert_array_equal(pretty.Count, counts)
    np.testing.assert_allclose(pretty.Percentage, counts / len(bed), rtol=1e-12)
    report = {"status": "PASS", "global_anchors": len(bed),
              "class_counts": {k: int(v) for k, v in counts.items()},
              "scope": "Internal consistency and classification from released z-scores/distances; not a raw-data rerun"}
    developmental = folder / "monkey_cCRE_developmental_annotation.tsv"
    if developmental.exists():
        ids = pd.read_csv(developmental, sep="\t", usecols=["id"]).id
        if ids.duplicated().any() or not set(ids).issubset(set(bed.id)):
            raise ValueError("Developmental annotation is not a unique subset of global IDs")
        report["optional_developmental_rows"] = len(ids)
        report["global_ids_without_developmental_rows"] = len(bed) - len(ids)
    if args.write_manifest:
        rows = []
        for path in sorted(folder.iterdir()):
            if path.is_file() and path.suffix.lower() in {".bed", ".txt", ".tsv", ".csv", ".png", ".gz"} and path.name != "MANIFEST.tsv":
                rows.append((path.name, path.stat().st_size, sha256(path)))
        write_tsv(pd.DataFrame(rows, columns=["file", "size_bytes", "sha256"]), folder / "MANIFEST.tsv")
        with atomic_path(folder / "validation_report.json") as tmp:
            tmp.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    cli(main)

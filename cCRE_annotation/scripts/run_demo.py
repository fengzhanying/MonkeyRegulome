#!/usr/bin/env python3
"""Run the portable synthetic core demo and compare against independently specified outputs."""
import argparse
import json
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from common import ROOT, MARKS, atomic_path, cli

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", type=Path, default=ROOT / "demo_output")
    args = p.parse_args()
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Demo output directory must be new or empty")
    fixture = ROOT / "demo"
    started = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir()
    shutil.copyfile(fixture / "input/experiments.tsv", out / "logs/experiments.tsv")
    shutil.copytree(fixture / "input/signals", out / "02.signals")
    commands = [
        ["00_prepare_masterlist.py", "--masterlist", str(fixture / "input/masterlist.bed")],
        ["01_generate_tss.py", "--gtf", str(fixture / "input/transcripts.gtf")],
        ["03_calculate_zscores.py"],
        ["04_classify_ccres.py", "--nearest-tss", str(fixture / "input/nearest_tss.tsv"),
         "--distance-mode", "legacy"],
    ]
    for script, *options in commands:
        subprocess.run([sys.executable, str(ROOT / "scripts" / script),
                        "--workdir", str(out), *options], check=True)
    for name, folder in [("monkey_TSS.bed", "01.masterlist"),
                          ("monkey_cCREs.bed", "04.classification"),
                          ("cCRE_class_stats.txt", "04.classification")]:
        actual = (out / folder / name).read_text().splitlines()
        expected = (fixture / "expected" / name).read_text().splitlines()
        if actual != expected:
            raise ValueError(f"Demo output mismatch: {name}")
    expected = pd.read_csv(fixture / "expected/max_zscores.tsv", sep="\t").set_index("name")
    for mark in MARKS:
        actual = pd.read_csv(out / "03.zscores" / f"{mark}_max_zscore.txt", sep="\t").set_index("name")
        np.testing.assert_allclose(actual.loc[expected.index, "max_zscore"],
                                   expected[mark], rtol=1e-6, atol=1e-7)
    report = {"status": "PASS", "anchors": 20, "classes": 6,
              "elapsed_seconds": round(time.perf_counter() - started, 3),
              "platform": platform.platform(), "python": platform.python_version(),
              "numpy": np.__version__, "pandas": pd.__version__,
              "scope": "Core demo; precomputed synthetic signals and nearest-TSS fixture. No bigWig extraction or SLURM."}
    with atomic_path(out / "demo_report.json") as tmp:
        tmp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    cli(main)

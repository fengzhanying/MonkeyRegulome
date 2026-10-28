"""Input validation, atomic output and provenance shared by the cCRE stages."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

MARKS = ("ATAC", "H3K4me3", "H3K27ac", "CTCF")
CLASSES = ("PLS", "pELS", "dELS", "CA-H3K4me3", "CA-CTCF", "CA")
ROOT = Path(__file__).resolve().parents[1]

def parser(description):
    result = argparse.ArgumentParser(description=description)
    result.add_argument("--workdir", type=Path, default=os.environ.get("CCRE_WORKDIR"),
                        required=not os.environ.get("CCRE_WORKDIR"))
    return result

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

@contextmanager
def atomic_path(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        yield temporary
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)

def write_tsv(frame, path, **kwargs):
    with atomic_path(path) as temporary:
        frame.to_csv(temporary, sep="\t", index=False, lineterminator="\n", **kwargs)

def read_anchors(path):
    frame = pd.read_csv(path, sep="\t", header=None, dtype={0: str, 3: str})
    if frame.shape[1] != 4 or frame.empty:
        raise ValueError(f"{path}: expected nonempty BED4")
    frame.columns = ["chr", "start", "end", "id"]
    if frame.isna().any().any() or frame.id.duplicated().any():
        raise ValueError(f"{path}: missing fields or duplicate anchor IDs")
    for column in ("start", "end"):
        values = pd.to_numeric(frame[column], errors="raise")
        if not np.isfinite(values).all() or (values % 1 != 0).any():
            raise ValueError(f"{path}: invalid integer coordinates")
        frame[column] = values.astype("int64")
    if (frame.start < 0).any() or (frame.end <= frame.start).any():
        raise ValueError(f"{path}: invalid BED intervals")
    return frame

def read_manifest(path):
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(frame.columns) != ["sample", "mark", "bigwig"] or frame.empty:
        raise ValueError(f"{path}: expected sample, mark, bigwig columns")
    if frame['sample'].duplicated().any() or not frame['sample'].str.fullmatch(r"[A-Za-z0-9_.-]+").all():
        raise ValueError("Sample names must be unique, nonempty safe file names")
    if set(frame.mark) != set(MARKS):
        raise ValueError(f"Manifest must contain all and only these marks: {MARKS}")
    return frame

def read_signal(path, ids):
    frame = pd.read_csv(path, sep="\t", header=None)
    if frame.shape[1] != 6:
        raise ValueError(f"{path}: expected six bigWigAverageOverBed columns")
    frame.columns = ["name", "size", "covered", "sum", "mean0", "mean"]
    if frame.name.duplicated().any() or set(frame.name) != set(ids):
        raise ValueError(f"{path}: anchor IDs are duplicated, missing or unexpected")
    numeric = frame.iloc[:, 1:].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric.to_numpy()).all():
        raise ValueError(f"{path}: non-finite signal values")
    if (numeric['size'] <= 0).any() or (numeric.covered < 0).any() or (numeric.covered > numeric['size']).any():
        raise ValueError(f"{path}: invalid size/coverage")
    frame.iloc[:, 1:] = numeric
    return frame.set_index("name").loc[ids]

def record(workdir, stage, inputs, outputs, parameters=None):
    def inventory(paths):
        return [{"name": Path(p).name, "size_bytes": Path(p).stat().st_size,
                 "sha256": sha256(p)} for p in paths]
    payload = {
        "stage": stage, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "code_version": (ROOT / "VERSION").read_text().strip(),
        "platform": platform.platform(), "python": platform.python_version(),
        "numpy": np.__version__, "pandas": pd.__version__,
        "parameters": parameters or {}, "inputs": inventory(inputs),
        "outputs": inventory(outputs),
        "source_sha256": {p.name: sha256(p) for p in sorted((ROOT / "scripts").glob("*.py"))},
    }
    with atomic_path(Path(workdir) / "provenance" / f"{stage}.json") as temporary:
        temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

def cli(main):
    try:
        main()
    except (ValueError, OSError, pd.errors.ParserError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

#!/usr/bin/env python3
"""Build a reviewer ZIP from an explicit allowlist, with internal SHA256 checksums."""
import argparse
import csv
import io
import zipfile
from pathlib import Path
from common import ROOT, atomic_path, cli, sha256

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    version = (ROOT / "VERSION").read_text().strip()
    target = args.output or ROOT / "dist" / f"cCRE_annotation-{version}-review.zip"
    if target.exists():
        raise ValueError(f"Archive exists; choose a new --output: {target}")
    top = ["README.md", "LICENSE", "VERSION", "CHANGELOG.md", "requirements.txt",
           "environment.yml", ".gitignore", ".gitattributes"]
    paths = [ROOT / name for name in top]
    for directory in ["scripts", "demo", "tests", "docs", "results"]:
        paths.extend(p for p in (ROOT / directory).rglob("*")
                     if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
    paths.extend([ROOT / "config/config.template.sh", ROOT / "config/experiments.example.tsv"])
    rows = [(p.relative_to(ROOT).as_posix(), sha256(p)) for p in sorted(paths)]
    manifest = io.StringIO()
    writer = csv.writer(manifest, delimiter="\t", lineterminator="\n")
    writer.writerow(["file", "sha256"])
    writer.writerows(rows)
    with atomic_path(target) as tmp:
        with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for path in sorted(paths):
                archive.write(path, "cCRE_annotation/" + path.relative_to(ROOT).as_posix())
            archive.writestr("cCRE_annotation/SHA256SUMS.tsv", manifest.getvalue())
        with zipfile.ZipFile(tmp) as archive:
            failed = archive.testzip()
            if failed:
                raise ValueError(f"ZIP CRC failure: {failed}")
    print(f"Archive: {target}")
    print(f"Files: {len(paths) + 1}; bytes: {target.stat().st_size}")
    print(f"SHA256: {sha256(target)}")

if __name__ == "__main__":
    cli(main)

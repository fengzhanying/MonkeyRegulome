"""Regression tests for scientific thresholds, malformed input and demo execution."""
import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from common import atomic_path, read_anchors, read_manifest, read_signal

def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / (name + ".py"))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

zscores = module("03_calculate_zscores")
classification = module("04_classify_ccres")

class CoreTests(unittest.TestCase):
    def test_threshold_is_strict_and_priority_is_preserved(self):
        base = dict(z_h3k4me3=1.64, z_h3k27ac=1.64, z_ctcf=1.64,
                    tss_core=True, tss_proximal=True)
        self.assertEqual(classification.classify(base), "CA")
        base.update(z_h3k4me3=2, z_h3k27ac=3, z_ctcf=4)
        self.assertEqual(classification.classify(base), "PLS")
        base.update(tss_core=False)
        self.assertEqual(classification.classify(base), "pELS")
        base.update(tss_proximal=False)
        self.assertEqual(classification.classify(base), "dELS")
        base.update(z_h3k27ac=0)
        self.assertEqual(classification.classify(base), "CA-H3K4me3")
        base.update(z_h3k4me3=0)
        self.assertEqual(classification.classify(base), "CA-CTCF")

    def test_standardization_ddof_and_negative_clipping(self):
        z = zscores.standardize([-5, 0, 0, 100])
        np.testing.assert_allclose(z, [-0.5, -0.5, -0.5, 1.5])
        self.assertEqual(z.dtype, np.float32)

    def test_constant_zero_signal(self):
        np.testing.assert_array_equal(zscores.standardize([0, 0, 0]), [0, 0, 0])

    def test_bad_normalization_inputs(self):
        for values in [[1], [0, float("nan")], [0, float("inf")]]:
            with self.assertRaises(ValueError):
                zscores.standardize(values)

    def test_signed_minus_one_modes_and_missing_tss(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "near.tsv"
            anchors = pd.DataFrame([["chr1", 100, 101, "one"]], columns=["chr", "start", "end", "id"])
            path.write_text("chr1\t100\t101\tone\tchr1\t99\t100\tGENE\t0\t+\t-1\n")
            legacy = classification.parse_nearest(path, anchors, "legacy")
            corrected = classification.parse_nearest(path, anchors, "corrected")
            self.assertEqual(legacy.center_dist.iloc[0], 1e12)
            self.assertEqual(corrected.center_dist.iloc[0], 1)
            path.write_text("chr1\t100\t101\tone\tnone\t-1\t-1\t.\t-1\t.\t-1\n")
            self.assertEqual(classification.parse_nearest(path, anchors, "corrected").center_dist.iloc[0], 1e12)

    def test_distance_boundaries(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "near.tsv"
            anchors = pd.DataFrame([["chr1", 1200, 1201, "a"], ["chr1", 3000, 3001, "b"],
                                    ["chr1", 3001, 3002, "c"]],
                                   columns=["chr", "start", "end", "id"])
            path.write_text("".join(f"chr1\t{s}\t{e}\t{i}\tchr1\t1000\t1001\tG\t0\t+\t{1000-s}\n"
                                    for _, s, e, i in anchors.itertuples(index=False, name=None)))
            distances = classification.parse_nearest(path, anchors, "corrected").center_dist
            self.assertEqual(list(distances <= 200), [True, False, False])
            self.assertEqual(list(distances <= 2000), [True, True, False])

    def test_duplicate_and_invalid_anchors(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "anchors.bed"
            for text in ["chr1\t0\t1\ta\nchr1\t0\t1\ta\n", "chr1\t2\t1\ta\n"]:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_anchors(path)

    def test_signal_ids_and_nan_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "signal.tab"
            for text in ["a\t1\t1\t1\t1\t1\na\t1\t1\t1\t1\t1\n",
                         "b\t1\t1\t1\t1\t1\n",
                         "a\t1\t0\t0\tNaN\t0\n"]:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_signal(path, ["a"])

    def test_manifest_requires_all_marks(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "manifest.tsv"
            path.write_text("sample\tmark\tbigwig\none\tATAC\tx\n")
            with self.assertRaises(ValueError):
                read_manifest(path)

    def test_failed_atomic_write_keeps_original(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "result.txt"
            path.write_text("complete")
            with self.assertRaises(RuntimeError):
                with atomic_path(path) as tmp:
                    tmp.write_text("partial")
                    raise RuntimeError("simulated extraction failure")
            self.assertEqual(path.read_text(), "complete")
            self.assertEqual(len(list(Path(folder).iterdir())), 1)

    def test_demo_and_failure_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder) / "demo"
            command = [sys.executable, str(ROOT / "scripts/run_demo.py"), "--outdir", str(out)]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            # Required mark cannot silently become all zeros at classification.
            (out / "03.zscores/CTCF_max_zscore.txt").unlink()
            command = [sys.executable, str(ROOT / "scripts/04_classify_ccres.py"),
                       "--workdir", str(out), "--nearest-tss", str(ROOT / "demo/input/nearest_tss.tsv")]
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            # Nor may a partial extraction silently drop an experiment.
            (out / "02.signals/demo_CTCF.tab").unlink()
            command = [sys.executable, str(ROOT / "scripts/03_calculate_zscores.py"), "--workdir", str(out)]
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)

if __name__ == "__main__":
    unittest.main()

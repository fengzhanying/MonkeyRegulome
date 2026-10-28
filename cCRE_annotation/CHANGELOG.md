# Changes

## 1.1.0 - 2026-10-28 submission preparation

- Add a portable synthetic core demo, independent expected outputs and regression tests.
- Add versioned dependencies, installation/runtime guidance, MIT license copy,
  exact method pseudocode and mapping to the supplied Nature checklist.
- Validate experiments, all four marks, finite signals and exact anchor IDs.
  Missing data no longer silently produces zero scores or skipped experiments.
- Replace path hardcoding with CLI/environment inputs and explicit experiment manifests.
- Stream per-mark maxima using the original float32 z-score representation;
  full matrices require --save-matrix.
- Require a fresh production workdir; use atomic signal/TSS outputs, no implicit
  signal/TSS cache reuse, and record hashes, versions and parameters.
- Make local/SLURM execution explicit; use afterok dependencies and report queued
  jobs as pending rather than completed.
- Expose legacy (default) and corrected missing-TSS modes. Legacy preserves the
  historical signed -1 defect; corrected mode can change biological labels.
  Interval-based nearest-TSS selection is unchanged.
- Sort newly generated annotations by genomic coordinates. Archived results
  retain their original contents/order.
- Verify all 702,876 archived classes against their stored scores/distances,
  coordinates and summaries; refresh stale result inventory.
- Document the optional 687,576-row developmental subset and custom BED-like
  five-column schema accurately.

These changes do not establish a raw-data rerun or validation of the production
Linux/SLURM/R environment. See docs/validation.md.

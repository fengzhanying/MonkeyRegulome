# Validation record

Date: 2026-10-28. Scope: submission preparation of cCRE_annotation 1.1.0.
This record distinguishes executed checks from intended production environments.

## Executed environment

- Windows 11 x86-64, build 26200.
- Python 3.12.14, NumPy 2.3.5, pandas 3.0.1.
- Fresh venv dependencies: python-dateutil 3.0.0, six 1.17.0, tzdata 2026.5.
- Processor reported as Intel64 Family 6 Model 183 Stepping 1; 28 logical CPUs.
- No GPU used. Installed RAM was not measured; the 4 GB recommendation is a
  conservative capacity guideline, not a measured peak-memory result.
- Native Bash, bedtools, UCSC tools, R and SLURM were unavailable for execution.

## Results

| Check | Result |
| --- | --- |
| Clean isolated dependency install from requirements | Passed; 31.179 seconds for package installation, about 3 seconds for venv creation |
| pip check in fresh venv | No broken requirements |
| Initial core demo in provided runtime | Passed in 5.503 seconds |
| Core demo in clean venv | Passed in 2.132 seconds |
| Regression suite in clean venv | 11 tests passed in 3.514 seconds |
| Existing global result validation | Passed for 702,876 unique anchors |
| Optional developmental ID subset | 687,576 unique IDs; 15,300 global IDs absent |

Times describe this host/network and include no native track extraction.
A dependency download cache and machine load affect them. They do not establish
timing on every normal desktop. Typical installation/runtime expectations in
the README are explicitly estimates.

The demo checks the actual anchor-preparation, TSS, z-score and classification
stages against known labels, counts, TSS records and closed-form standardized
values. It uses precomputed synthetic signals and a nearest-TSS fixture.
The regression suite also tests strict 1.64 comparisons and class priority,
inclusive 200/2000 bp distances, signed -1 versus absent TSS records, constant
signals, negative-signal clipping, non-finite values, duplicate/missing IDs,
missing marks, stale output refusal and preservation of complete files after
an interrupted atomic write.

The global validator compares both coordinate tables, all stored class labels
against the stored z-scores/distance flags, six class counts, readable summary
proportions, and the developmental ID subset. This is internal consistency
validation; it cannot independently validate the source signals or TSS choice.

## Remaining external validation

A native Linux run (including the provided conda recipe), actual bigWig extraction,
bedtools closest behavior on the original inputs, shell/SLURM execution, and
the optional R composition plot have not been executed in this validation.
The original cluster software versions and complete input files were unavailable.

No fresh full-data corrected-mode results have been generated. The existing
global biological data have not been replaced. An independent colleague
installation and reviewer access checks remain author actions, as recommended
by the supplied checklist.

The first sandboxed test attempt was blocked by temporary-directory permissions;
the suite was subsequently run successfully outside that sandbox. This was a
local execution-environment issue, not a passing result for the failed attempt.

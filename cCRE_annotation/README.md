# Macaque global cCRE annotation

Global candidate cis-regulatory element (cCRE) annotation for
[MonkeyRegulome](https://github.com/fengzhanying/MonkeyRegulome).
This is an ENCODE-inspired, ATAC-anchor-based workflow, not the ENCODE registry
pipeline. Version: **1.1.0**. License: [MIT](LICENSE).

Start with the self-contained core demo below. It needs no cluster, genome
download, bigWig files or network after installation.

## 1. System requirements

| Component | Version / status |
| --- | --- |
| Core demo and Python stages | Tested on Windows 11, build 26200, x86-64 |
| Python | Tested: 3.12.14 |
| NumPy / pandas | Tested: 2.3.5 / 3.0.1; see requirements.txt |
| Production shell pipeline | Target: Linux x86-64 (Ubuntu 22.04 or compatible HPC); not executed on this Windows host |
| Bash | Production target: 5.1 or newer |
| bedtools | Production recipe: 2.31.1; native invocation not tested here |
| UCSC bigWigAverageOverBed | Production recipe: 482; native invocation not tested here |
| SLURM | Optional site-managed scheduler; no tested version available |
| Optional plotting | R 4.5.2, ggplot2 4.0.0, scales 1.4.0 are recipe targets, not verified in this validation |

No GPU or other non-standard hardware is needed. A normal desktop with 4 GB RAM
is sufficient for the core demo. For the full atlas, start with 8-32 GB RAM and
disk space for the input tracks and signal tables; resource use depends on the
number of experiments. The historical 30 GB SLURM request is a starting setting,
not a measured minimum. Maximum z-scores are now accumulated one experiment at
a time; full per-experiment matrices are optional.

Actual tests and their limits are recorded in [validation.md](docs/validation.md).
The versions used in the original cluster analysis were not recorded in the
available files; the versions above must not be described as those original versions.

## 2. Installation

Run commands from this directory. Install Python 3.12, then create an isolated
environment. No compilation is required for the Python demo.

Windows PowerShell:

~~~powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts/run_demo.py
~~~

Linux:

~~~bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/run_demo.py
~~~

Allow approximately 1-5 minutes for a typical desktop installation with binary
wheels and a working network (an estimate, not a universal benchmark). See the
validation record for the measured installation on this host.

For production on Linux, the provided conda recipe adds the external tools:

~~~bash
conda env create -f environment.yml
conda activate ccre
~~~

The Linux environment is a versioned installation recipe, not a verified
platform lockfile. A typical solve/download may take 5-15 minutes, depending on
the network and solver; this estimate has not been measured here.

Optional plot dependencies can be installed in a separate conda environment:

~~~bash
conda create -n ccre-plot -c conda-forge r-base=4.5.2 r-ggplot2=4.0.0 r-scales=1.4.0
conda activate ccre-plot
Rscript scripts/05_plot_ccre_composition.R results/global/cCRE_class_stats.txt plot_output
~~~

This optional recipe has not been executed here. It is not needed for the demo
or global classification. The plot script records counts and proportions, and
writes PDF and PNG outputs in the chosen output directory.

## 3. Demo

~~~bash
python scripts/run_demo.py --outdir demo_output
~~~

Use the Python executable from your environment. On Windows without activation,
replace "python" with ".\.venv\Scripts\python.exe".

The demo uses 20 synthetic anchors, a tiny GTF and four synthetic signal tables.
It runs the actual anchor-preparation, TSS, z-score and classification scripts.
It supplies a precomputed nearest-TSS fixture, so it does **not** test
bigWigAverageOverBed, bedtools, shell orchestration or SLURM.

Expected outputs:

| Class | Count |
| --- | ---: |
| PLS | 1 |
| pELS | 1 |
| dELS | 1 |
| CA-H3K4me3 | 1 |
| CA-CTCF | 1 |
| CA | 15 |

Outputs appear under demo_output/01.masterlist, 03.zscores and 04.classification,
with stage provenance JSON and demo_report.json. The runner checks all labels,
TSS records, counts and maximum z-scores against demo/expected; failure returns
a nonzero exit status. Expect less than one minute on a normal desktop.
The first measured run on this host took 5.503 seconds.

Use a new --outdir for each run. Existing nonempty directories are rejected.
See [demo/README.md](demo/README.md) for fixture definitions and expected numbers.

Run regression tests:

~~~bash
python -m unittest discover -s tests -v
~~~

## 4. Run on your data

Inputs must use the same genome assembly and chromosome naming. The historical
GTF filename specifies Macaca mulatta Mmul_10, Ensembl release 114. Confirm that
the ATAC master list and every track use that assembly; renaming chromosomes
cannot convert between assemblies.

1. Supply the already merged global ATAC peak universe as BED (zero-based,
   half-open; at least three columns). Peak calling and peak merging are upstream.
2. Supply a plain or gzipped GTF with transcript records.
3. Create a TSV using config/experiments.example.tsv: sample, mark, bigwig.
   Use unique sample IDs and all four marks: ATAC, H3K4me3, H3K27ac, CTCF.
   Relative bigWig paths are resolved relative to this TSV.
4. Copy config/config.template.sh to config/config.sh and set absolute paths.
   Set CCRE_EXPERIMENTS to the manifest path and CCRE_WORKDIR to a new/empty directory.
5. Activate the Linux environment and run:

~~~bash
bash scripts/run_global_ccre.sh config/config.sh
~~~

The default executor is local (serial extraction). Set CCRE_EXECUTOR=slurm
explicitly to submit an array plus a downstream job dependent on successful
extraction. Submission is asynchronous; inspect the jobs and logs before using
results. All paths and the activated environment must be accessible on compute nodes.

As an alternative to an explicit manifest, the driver discovers the historical
patterns ATAC-*.fc.signal.bw, CUT-Tag-*-H3K4me3.fc.signal.bw,
CUT-Tag-*-H3K27ac.fc.signal.bw and CUT-Tag-*-CTCF.fc.signal.bw in the two configured
track directories. Selected inputs are frozen in logs/experiments.tsv.

The driver refuses a nonempty workdir to prevent stale intermediates. Individual
stages can be rerun explicitly after diagnosing a failure; they validate their
inputs. Outputs and SHA256 provenance are written within the workdir.
The default z-score stage omits the large experiment matrix; add --save-matrix
when invoking that stage directly if it is required.

## 5. Method, outputs and published results

See [method_summary.md](docs/method_summary.md) for precise rules and pseudocode,
[code_summary.md](docs/code_summary.md) for script responsibilities, and
[results/global/README.md](results/global/README.md) for formats and provenance.

Files in `results/global/` include:

- `monkey_cCREs.bed`: global macaque cCRE coordinates and class labels.
- `monkey_cCREs_annotated.txt`: global macaque cCRE coordinates, nearest
  gene/TSS information, mark z-scores, and class labels.
- `cCRE_class_stats.txt`: class counts.
- `Monkey_cCRE_Composition_Stats.csv`: class counts and percentages using
  readable class names.
- `monkey_cCRE_developmental_annotation.tsv`: sample, tissue, and time-point
  support summary for global cCRE anchors.

The default distance mode is legacy for compatibility with the archived
classification. It preserves an upstream signed-distance -1 handling defect.
The optional corrected mode fixes that missing-TSS test; it can change labels.
Both modes retain the historical interval-nearest TSS selection. See the method
notes before using corrected results in a manuscript.

The existing global results contain 702,876 anchors and have been checked for
internal consistency. To repeat that check:

~~~bash
python scripts/verify_results.py
~~~

This check verifies coordinate tables, stored scores, classes and counts.
It is not a rerun from the original raw tracks. A complete reanalysis also needs
the original merged BED, exact GTF, selected experiment manifest and bigWigs;
those inputs are not bundled. Downloading the archived classification alone
does not reproduce signal extraction.

The developmental annotation table is optional and is not produced or consumed
by this global workflow. It contains 687,576 unique IDs, a subset of the global
universe. Missing rows must not be interpreted as inactive elements.

## 6. Reviewer package and availability

~~~bash
python scripts/make_submission.py
~~~

This creates a ZIP containing source, demo, expected outputs, tests, documentation,
license and the available global results, with a file checksum manifest.
It excludes virtual environments, personal config, caches and run directories.
See [upload_checklist.md](docs/upload_checklist.md).

Repository: https://github.com/fengzhanying/MonkeyRegulome/tree/main/cCRE_annotation
The submitted ZIP or an accessible immutable release should accompany review.
Local changes are not automatically published. Repository accessibility and an
independent colleague installation still require author verification.

[Nature checklist mapping](docs/nature_checklist.md) documents which attachment
requirements are addressed and which validation steps remain.

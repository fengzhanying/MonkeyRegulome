# Macaque cCRE annotation

This directory contains the global macaque candidate cis-regulatory element
(cCRE) annotation workflow used in MonkeyRegulome.

The workflow is ENCODE-inspired. It starts from the accumulated macaque ATAC
peak universe, extracts ATAC/H3K4me3/H3K27ac/CTCF signals, computes
per-experiment z-scores, and assigns each accessible anchor to one global cCRE
class.

## Directory layout

- `scripts/`: code for the global macaque cCRE workflow.
- `config/`: template environment file for local or cluster paths.
- `results/global/`: processed global cCRE result files that can be released
  with the repository.
- `docs/`: method notes and upload checklist.

## Main workflow

1. Prepare the macaque rDHS/cCRE anchor master list from the merged ATAC peak
   set.
2. Generate a TSS BED file from the macaque GTF annotation.
3. Extract mean signal over each anchor from ATAC-seq and CUT&Tag bigWig files.
4. Compute mark-wise maximum z-scores across experiments.
5. Classify global macaque cCREs as `PLS`, `pELS`, `dELS`, `CA-H3K4me3`,
   `CA-CTCF`, or `CA`.

## Required inputs

- Merged macaque ATAC peak BED file.
- Macaca mulatta GTF annotation, matching the genome assembly used for peak
  calling and signal tracks.
- Fold-change signal bigWig files for ATAC, H3K4me3, H3K27ac, and CTCF.

Expected bigWig filename patterns:

- `ATAC-*.fc.signal.bw`
- `CUT-Tag-*-H3K4me3.fc.signal.bw`
- `CUT-Tag-*-H3K27ac.fc.signal.bw`
- `CUT-Tag-*-CTCF.fc.signal.bw`

## Software

- `bash`
- `python3`
- `pandas`
- `numpy`
- `bedtools`
- UCSC `bigWigAverageOverBed`
- Optional: `R`, `ggplot2`, `scales` for the composition plot
- Optional: SLURM `sbatch` for array-job execution

## Run

Copy and edit the configuration template:

```bash
cp config/config.template.sh config/config.sh
vim config/config.sh
```

Then run:

```bash
bash scripts/run_global_ccre.sh config/config.sh
```

The workflow writes intermediate files under `$CCRE_WORKDIR`:

- `01.masterlist/monkey_rDHS.bed`
- `01.masterlist/monkey_TSS.bed`
- `02.signals/*.tab`
- `03.zscores/*_max_zscore.txt`
- `04.classification/monkey_cCREs_annotated.txt`
- `04.classification/monkey_cCREs.bed`
- `04.classification/cCRE_class_stats.txt`

## Published global results

The preferred files to publish under `results/global/` are:

- `monkey_cCREs.bed.gz`: global macaque cCRE coordinates and class labels.
- `monkey_cCREs_annotated.txt.gz`: global macaque cCRE coordinates, nearest
  gene/TSS information, mark z-scores, and class labels.
- `cCRE_class_stats.txt`: class counts.
- `Monkey_cCRE_Composition_Stats.csv`: class counts and percentages using
  readable class names.
- `monkey_cCRE_developmental_annotation.tsv`: sample, tissue, and time-point
  support summary for global cCRE anchors.

Large files should be gzip-compressed before upload. If any result file exceeds
GitHub's file-size limits, keep a small manifest in this directory and release
the large file through GitHub Releases, Zenodo, or another permanent data
repository.

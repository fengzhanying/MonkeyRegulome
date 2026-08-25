# cCRE annotation upload checklist

Upload only the global cCRE annotation package to the MonkeyRegulome repository:

```text
cCRE_annotation/
  README.md
  config/
    config.template.sh
  scripts/
    run_global_ccre.sh
    00_prepare_masterlist.sh
    01_generate_tss.py
    02_extract_signals.sh
    03_calculate_zscores.py
    04_classify_ccres.py
    05_plot_ccre_composition.R
  results/
    global/
      README.md
      Monkey_cCRE_Composition_Stats.csv
      monkey_cCRE_developmental_annotation.tsv
      monkey_cCREs.bed.gz
      monkey_cCREs_annotated.txt.gz
      cCRE_class_stats.txt
  docs/
    method_summary.md
    upload_checklist.md
```

Do not upload:

- Raw ATAC-seq or CUT&Tag bigWig files.
- Intermediate signal tables in `02.signals/`.
- Intermediate z-score matrices in `03.zscores/`.
- SLURM logs.
- ENCODE reference PDFs.
- Sample-specific, tissue-specific, or cross-species comparison scripts for
  this first global-only release.

Before pushing to GitHub:

- Confirm that the coordinate assembly is stated consistently.
- Confirm that `monkey_cCREs.bed.gz` and
  `monkey_cCREs_annotated.txt.gz` are generated from the same run.
- Add checksums for large result files.
- Keep files larger than GitHub's limit out of the normal git history and
  publish them through a release or data repository.

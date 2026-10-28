# Code inventory

| File in scripts/ | Responsibility |
| --- | --- |
| run_global_ccre.sh | Configured Linux entry point; explicit local/SLURM execution |
| prepare_run.py | Require fresh output directory, validate inputs, freeze manifest |
| 00_prepare_masterlist.py / .sh | Normalize and validate merged ATAC anchors |
| 01_generate_tss.py | GTF/gzipped GTF to deduplicated BED6 TSS |
| extract_signal.py / 02_extract_signals.sh | Extract and validate one bigWig; atomic output |
| 03_calculate_zscores.py | Validated log1p(mean0) normalization and streaming mark maxima |
| 04_classify_ccres.py | TSS selection, explicit legacy/corrected modes, ordered six-class rules |
| run_downstream.sh | Run z-score and classification stages after successful extraction |
| 05_plot_ccre_composition.R | Optional composition plot and human-readable summary |
| common.py | Input validation, atomic file publication, SHA256 and provenance |
| run_demo.py | Execute real core stages on fixtures; compare expected outputs |
| verify_results.py | Check existing global results against each other and their stored scores |
| make_submission.py | Package reviewer ZIP with a per-file checksum manifest |

All Python command-line entry points support --help. Stage scripts accept
--workdir or the CCRE_WORKDIR environment variable. GTF and master-list stages
also accept explicit input arguments. The production shell entry point takes
a trusted config file. No network access is used by the analysis scripts.

Tests are in tests/test_core.py. Input and independently specified expected
demo files are under demo/. The optional large developmental table is neither
an input nor an output of these global scripts.

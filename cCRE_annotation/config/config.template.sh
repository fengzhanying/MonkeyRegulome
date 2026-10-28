# Source only trusted config files. Use absolute paths accessible on compute nodes.
export CCRE_WORKDIR="/path/to/new_empty_ccre_workdir"
export CCRE_MASTERLIST="/path/to/Monkey_dTAC.bed"
export CCRE_GTF="/path/to/Macaca_mulatta.Mmul_10.114.chr.gtf"

# Preferred: explicit TSV (sample, mark, bigwig). Relative bigwig paths are
# resolved against the manifest directory. Uncomment to use it.
# export CCRE_EXPERIMENTS="/path/to/experiments.tsv"

# Otherwise use the historical filename patterns described in the README.
export CCRE_ATAC_BW_DIR="/path/to/atac/bw"
export CCRE_CUTTAG_BW_DIR="/path/to/cuttag/bw"

# Activate your environment before running, or supply an optional conda setup.
export CCRE_CONDA_SH=""
export CCRE_CONDA_ENV="ccre"
export CCRE_PYTHON="python3"
export CCRE_EXECUTOR="local"  # local or slurm (no automatic submission)
export CCRE_DISTANCE_MODE="legacy"  # historical reproduction; see method_summary.md

export CCRE_ARRAY_CONCURRENCY="50"
export CCRE_SIGNAL_MEM="5G"
export CCRE_SIGNAL_TIME="00:30:00"
export CCRE_DOWNSTREAM_MEM="30G"
export CCRE_DOWNSTREAM_TIME="04:00:00"

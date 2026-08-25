# Template configuration for macaque global cCRE annotation.
# Copy this file to config/config.sh and edit all paths before running.

export CCRE_WORKDIR="/path/to/macaque_cCRE_workdir"
export CCRE_MASTERLIST="/path/to/Monkey_dTAC.bed"
export CCRE_GTF="/path/to/Macaca_mulatta.Mmul_10.114.chr.gtf"

export CCRE_ATAC_BW_DIR="/path/to/atac/bw"
export CCRE_CUTTAG_BW_DIR="/path/to/cuttag/bw"

# Optional conda activation. Leave CCRE_CONDA_SH empty if the environment is
# already active.
export CCRE_CONDA_SH="/path/to/anaconda3/etc/profile.d/conda.sh"
export CCRE_CONDA_ENV="ucsctools"

# SLURM settings. Used only when sbatch is available.
export CCRE_ARRAY_CONCURRENCY="50"
export CCRE_SIGNAL_MEM="5G"
export CCRE_SIGNAL_TIME="00:30:00"
export CCRE_DOWNSTREAM_MEM="30G"
export CCRE_DOWNSTREAM_TIME="04:00:00"

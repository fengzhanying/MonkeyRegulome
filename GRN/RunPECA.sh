#!/bin/bash
#SBATCH --job-name=PECA
#SBATCH --partition=CU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --array=2-139
#SBATCH --error=./Log/PECA_%a.err
#SBATCH --output=./Log/PECA_%a.out

Name=`cat ../MakePrior/SampleNameFile.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
source PECA.sh ${Name}

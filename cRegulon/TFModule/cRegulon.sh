#!/bin/bash

#SBATCH -p CU
#SBATCH --job-name=cRegulon
#SBATCH --output=./Log/cRegulon_%a.out
#SBATCH --error=./Log/cRegulon_%a.err
#SBATCH --array=5-30
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --time=01-00:00:00

/lustre/home/zhangfy/Pipeline/Anaconda3/bin/python cRegulon.py --k $SLURM_ARRAY_TASK_ID

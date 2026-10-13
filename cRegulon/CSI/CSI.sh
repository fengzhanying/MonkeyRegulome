#!/bin/bash
#SBATCH --job-name=CSI
#SBATCH --output=CSI.out
#SBATCH --error=CSI.err
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -p CU
#SBATCH --cpus-per-task=16

/lustre/home/zhangfy/Pipeline/Anaconda3/bin/python CSI.py

#!/bin/sh

#SBATCH --job-name="ump"
#SBATCH --gres=lscratch:128
#SBATCH --cpus-per-task=128
#SBATCH --mem=40g
#SBATCH --mail-type=ALL

/data/kanferg/conda_v1/envs/spatiocore_spatial_env/bin/python3 run_ump.py
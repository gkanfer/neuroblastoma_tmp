#!/bin/sh

#SBATCH --job-name="liedne"
#SBATCH --gres=lscratch:128
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=60g
#SBATCH --mail-type=ALL

/data/kanferg/conda_v1/envs/spatiocore_spatial_env/bin/python3 run_lieden.py
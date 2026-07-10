#!/bin/sh
#SBATCH --job-name="sdata gener"
#SBATCH --gres=lscratch:128
#SBATCH --cpus-per-task=64
# #SBATCH --time=10:00:00
#SBATCH --mem=360g
#SBATCH --mail-type=ALL

/data/kanferg/conda_v1/envs/spatialdata_env/bin/python3 run_sdata_highres_v1.py
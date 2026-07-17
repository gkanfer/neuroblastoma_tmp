#!/bin/bash
#SBATCH --job-name="sdata gener"
#SBATCH --gres=lscratch:128
#SBATCH --cpus-per-task=64
# #SBATCH --time=10:00:00
#SBATCH --mem=200g 
# SBATCH --mem=360g 

#SBATCH --mail-type=ALL

source myconda_v1
mamba activate spatialdata_env
python run_sdata_highres_v1.py
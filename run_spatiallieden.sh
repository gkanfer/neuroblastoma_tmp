#!/bin/sh

#SBATCH --job-name="nb ratio sl"
#SBATCH --gres=lscratch:128
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=60g
#SBATCH --mail-type=ALL

source myconda 
mamba activate spatialleiden-cupy 

/data/kanferg/conda/envs/spatialleiden-cupy/bin/python3 run_spatiallieden_ratio.py
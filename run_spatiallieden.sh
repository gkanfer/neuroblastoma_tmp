#!/bin/sh

#SBATCH --job-name="nb spatiallieden"
#SBATCH --gres=lscratch:128
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=60g
#SBATCH --mail-type=ALL

source myconda 
mamba activate spatialleiden-cupy 

/data/kanferg/conda/envs/spatialleiden-cupy/bin/python3 run_spatiallieden.py
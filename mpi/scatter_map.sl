#!/bin/bash -e

#SBATCH --job-name=scatter
#SBATCH --error=scatter-%j.error
#SBATCH --output=scatter-%j.output
#SBATCH --time=00:10:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1

ml Python Boost
ml forge
map --profile srun python scatter.py

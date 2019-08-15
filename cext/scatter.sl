#!/bin/bash -e

#SBATCH --job-name=scatter
#SBATCH --error=scatter-%j.error
#SBATCH --output=scatter-%j.output
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread

module load Python Boost
python setup.py build --force
srun time python scatter.py -checksum

#!/bin/bash -e

#SBATCH --job-name=scatter
#SBATCH --error=scatter-%j.error
#SBATCH --output=scatter-%j.output
#SBATCH --time=00:10:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1

module load Python # mahuika
module load Boost

srun time python scatter.py -checksum

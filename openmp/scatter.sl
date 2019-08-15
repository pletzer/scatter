#!/bin/bash -e

#SBATCH --job-name=scatter
#SBATCH --error=scatter-%j.error
#SBATCH --output=scatter-%j.output
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread

module load Python >& /dev/null    # mahuika, ignore elsewhere
module load Anaconda3 >& /dev/null # maui, ignore elsewhere
module load Boost

python setup.py build --force

srun python scatter.py -checksum -nx 256 -ny 256

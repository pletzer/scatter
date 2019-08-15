#!/bin/bash -e

#SBATCH --job-name=scatter
#SBATCH --error=scatter-%j.error
#SBATCH --output=scatter-%j.output
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load Python >& /dev/null    # mahuika, ignore elsewhere
module load Anaconda3 >& /dev/null # maui, ignore elsewhere
module load Boost

srun time python -m cProfile -o output.pstats scatter.py
gprof2dot --colour-nodes-by-selftime -f pstats output.pstats | \
    dot -Tpng -o output.png

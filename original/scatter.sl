#!/bin/bash -e

#SBATCH --job-name=scatter
#SBATCH --error=scatter-%j.error
#SBATCH --output=scatter-%j.output
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# load modules according to host
host=$(echo $(hostname) | awk -F. '{print $1}')
if [[ $host =~ *mauivlab.* ]]; then
    # maui_ancil
    module load Anaconda3 Boost
else
    # default is mahuika
    module load Python Boost
fi

srun time python scatter.py -checksum

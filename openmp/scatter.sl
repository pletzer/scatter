#!/bin/bash -e

#SBATCH --job-name=scatter
#SBATCH --error=scatter-%j.error
#SBATCH --output=scatter-%j.output
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread

# load modules according to host
host=$(echo $(hostname) | awk -F. '{print $1}')
if [[ $host =~ *mauivlab.* ]]; then
    # maui_ancil
    module load Anaconda3 Boost
else
    # default is mahuika
    module load Python Boost
fi

python setup.py build --force

srun python scatter.py -checksum -nx 256 -ny 256

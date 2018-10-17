#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=Hills_demultiplex
#SBATCH --output=/projects/bgmp/mhills/demultiplex/slurm_demultiplex.out
#SBATCH --error=/projects/bgmp/mhills/demultiplex/slurm_demultiplex.err
#SBATCH --time=0-09:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28

# Move to my demultiplex folder.
cd /projects/bgmp/mhills/demultiplex
# Purge loaded modules.
module purge
# Load Python 3.6.1
module load easybuild intel/2017a Python/3.6.1 icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 ifort/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 matplotlib/2.0.1-Python-3.6.1 
# Run the python script 'demulti.py'.
python demulti.py


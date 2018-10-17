#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=hills_deMULTI
#SBATCH --output=/projects/bgmp/mhills/Demulti/deMULTI.out
#SBATCH --error=/projects/bgmp/mhills/Demulti/deMULTI.err
#SBATCH --time=0-05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28

# Move to the folder containing the files.
cd /projects/bgmp/mhills/Demulti/

# Load required modules.
module purge
ml python3/3.6.1
# Call my python script.
/usr/bin/time -v python3 deMULTI.py -r1 1294_S1_L008_R1_001.fq -r2 1294_S1_L008_R4_001.fq -i1 \
1294_S1_L008_R2_001.fq -i2 1294_S1_L008_R3_001.fq -k indexes.txt -c 20


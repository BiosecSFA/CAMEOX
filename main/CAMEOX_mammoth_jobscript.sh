#!/bin/bash
#SBATCH -M mammoth
#SBATCH -p pbatch
#SBATCH -t 1-00:00:00
#SBATCH -J CAMEOX
#SBATCH -N 1

# USAGE:
# sbatch CAMEOX_mammoth_jobscript.sh <CAMEOX_params_file>

export OMP_NUM_THREADS=256
cd /usr/workspace/kpath/cameos/CAMEOX/main/
echo 'CAMEOX starts...'
/usr/gapps/julia/bin/julia-1.8.3 -t 256 main.jl $1
echo 'CAMEOX run DONE!'

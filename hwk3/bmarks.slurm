#!/bin/bash
#SBATCH --partition=mic
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=4

setpkgs -a intel_cluster_studio_compiler

# rebuild executable 
make cleanall
make

# run using 8 MPI processes on system size of 1000
bash simulate-slurm.sh 16 1000

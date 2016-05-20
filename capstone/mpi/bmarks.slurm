#!/bin/bash
#SBATCH --partition=mic
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1

time bash simulate-slurm.sh 8

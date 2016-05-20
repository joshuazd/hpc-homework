#!/bin/bash

[ $# -lt 1 ] && { echo "Usage: $0 numprocs"; exit 999; }

numprocs=$1
time mpirun -n $numprocs ./image_denoise ~/bigmoon.noisy.bmp bigmoon.out.bmp 10 2.1 2.0
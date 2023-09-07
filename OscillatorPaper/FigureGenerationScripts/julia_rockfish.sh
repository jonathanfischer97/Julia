#!/bin/bash -l
#SBATCH --job-name=julia
#SBATCH --time=72:0:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --account=mjohn218_condo
#SBATCH --export=ALL

ml intel-mkl

julia -t 48 $1 > julia_output_$1.txt

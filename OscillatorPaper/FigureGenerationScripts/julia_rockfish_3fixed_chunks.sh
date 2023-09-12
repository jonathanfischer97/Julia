#!/bin/bash -l
#SBATCH --job-name=3fixed_chunks
#SBATCH --time=72:0:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --account=mjohn218_condo
#SBATCH --export=ALL

ml intel-mkl

# Concatenating all arguments and exported variables to form a unique name
unique_filename=$(echo "$@" "$START_IDX" "$END_IDX" | tr ' ' '_')

julia --threads=48 -- "$@" "$START_IDX $END_IDX" > "JULIA_OUT_${unique_filename}.txt"


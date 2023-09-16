#!/bin/bash

# Loop through values from 1 to 36 with steps of 4
for i in $(seq 1 4 36); do
  # Submit a new SLURM job with $i cores
  sbatch --ntasks=$i --export=JULIA_NUM_THREADS=$i run_core_test.sh
done


#!/bin/bash

# Loop through values from 1 to 36 with steps of 4
for i in $(seq 1 4 48); do
  # Generate a unique filename for each job's output
  unique_output_file="TimerResults/timeresults_${i}.txt"

  # Submit the job to SLURM with the specified number of threads
  sbatch --output=$unique_output_file --job-name="C_$i" --ntasks-per-node=$i run_core_test.sh

  # Add a separator line for readability in the output file
  echo "------------------" >> $unique_output_file
done

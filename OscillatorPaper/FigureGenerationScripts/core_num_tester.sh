#!/bin/bash

# Loop through values from 1 to 36 with steps of 4
for i in $(seq 1 4 36); do
  # Generate a unique filename for each job's output
  unique_output_file="TimerResults/timeresults_${i}.txt"

  # Submit the job to SLURM with the specified number of threads
  sbatch --job-name="CoreTest_$i" --ntasks-per-node=$i --wrap="julia --threads=$i CoreNumTester.jl >> $unique_output_file"

  # Add a separator line for readability in the output file
  echo "------------------" >> $unique_output_file
done

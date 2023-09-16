#!/bin/bash

# Output file where the Julia script outputs will be stored
output_file="TimerResults/timeresults.txt"

# Loop through values from 1 to 36 with steps of 4
for i in $(seq 1 4 36); do
  # Set the number of threads for Julia
  export JULIA_NUM_THREADS=$i

  # Run the Julia script and append output to file
  julia CoreNumTester.jl >> $output_file

  # Add a separator line for readability
  echo "------------------" >> $output_file
done

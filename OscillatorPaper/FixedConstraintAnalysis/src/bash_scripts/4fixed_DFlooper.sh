#!/bin/bash -l

# Define additional arguments
julia_file="SCRIPT_4FixedInitialConditionsCSVMaker.jl"
range_number=$1
population_number=$2
df_values=$3

# Remove brackets and quote marks, if any
df_values=${df_values//[\[\]\'\"]/}

# Convert comma-separated string to array
IFS=',' read -ra df_array <<< "$df_values"

# Loop through each number
for df in "${df_array[@]}"; do
    sbatch julia_rockfish.sh $julia_file $range_number $population_number $df
done



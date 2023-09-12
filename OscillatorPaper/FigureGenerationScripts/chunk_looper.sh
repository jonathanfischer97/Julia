#!/bin/bash -l

# Define additional arguments
julia_file="SCRIPT_3FixedParam&ICsFixedDF_CSVMaker.jl"
range_number=$1
population_number=$2
num_chunks=$3

# Generate chunks using Julia script
chunks=$(julia ./GenerateChunks_ROCKFISH.jl $num_chunks)

for chunk in $chunks; do
    start_idx=${chunk%:*}
    end_idx=${chunk#*:}
    sbatch --export=START_IDX=$start_idx,END_IDX=$end_idx julia_rockfish_3fixed_chunks.sh $julia_file $range_number $population_number
done


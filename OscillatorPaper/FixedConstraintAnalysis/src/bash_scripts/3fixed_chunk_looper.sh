#!/bin/bash -l

# Define additional arguments
cd "$(dirname "$0")"
julia_file="../data_generation_scripts/SCRIPT_3FixedParam&ICsFixedDF_CSVMaker.jl"
range_number=$1
population_number=$2
num_chunks=$3
output_dir="../../job_outputs"  # Relative to the directory where this script is called, so need to remember to call from inside the bash_scripts directory

# Create output directory if it does not exist
mkdir -p $output_dir

# Generate chunks using Julia script
chunks=$(julia "../data_generation_scripts/GenerateChunks_ROCKFISH.jl" $num_chunks)

# Initialize job counter
job_counter=1

# Loop through each chunk and submit a SLURM job
for chunk in $chunks; do
    start_idx=${chunk%:*}
    end_idx=${chunk#*:}
    unique_job_name="3F$job_counter"
    unique_output_file="$output_dir/jout_${unique_job_name}.txt"
    sbatch --job-name=$unique_job_name --output=$unique_output_file --export=START_IDX=$start_idx,END_IDX=$end_idx julia_rockfish_3fixed_chunks.sh $julia_file $range_number $population_number
    job_counter=$((job_counter + 1))
done




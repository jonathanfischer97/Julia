#!/bin/bash -l

#SBATCH --job-name=3fixed_chunks
#SBATCH --time=72:0:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --account=mjohn218_condo
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=jfisch27@jhu.edu


ml intel-mkl

# Create a unique filename based on all arguments and exported variables
unique_filename=$(echo "$@" "$START_IDX" "$END_IDX" | tr ' ' '_')

# Divide the start and end indices for 6 Julia processes
chunk_size=$(( (END_IDX - START_IDX + 1) / 6 ))

# Start 6 Julia processes in the background
for i in $(seq 0 5); do
    sub_start=$(( START_IDX + i * chunk_size ))
    sub_end=$(( sub_start + chunk_size - 1 ))
    if [ $i -eq 5 ]; then
        sub_end=$END_IDX
    fi
    julia --threads=8 -- "$@" "$sub_start $sub_end" > "JULIA_OUT_${unique_filename}_${i}.txt" &
done

# Wait for all background jobs to finish
wait



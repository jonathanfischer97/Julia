#!/bin/bash -l

#SBATCH --time=72:0:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --account=mjohn218_condo
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=jfisch27@jhu.edu


ml intel-mkl

# Divide the start and end indices for 6 subchunks for each 8-threaded Julia process
chunk_size=$(( (END_IDX - START_IDX + 1) / 6 ))

# Start 6 Julia processes in the background
for i in $(seq 0 5); do
    sub_start=$(( START_IDX + i * chunk_size ))
    sub_end=$(( sub_start + chunk_size - 1 ))
    if [ $i -eq 5 ]; then
        sub_end=$END_IDX
    fi
    mkdir -p $OUTFILE
    julia --threads=8 -- "$@" "$sub_start $sub_end" > "$OUTFILE/julia_out_${i}.txt" 2> "$OUTFILE/julia_err_${i}.txt" &
done

# Wait for all background jobs to finish
wait



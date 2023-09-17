#!/bin/bash
#SBATCH --job-name=CoreTest
#SBATCH --time=72:0:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --account=mjohn218_condo
#SBATCH --export=ALL
#SBATCH --output=TimerResults/timeresults_%j.txt


# Run the Julia script and write output to file
julia --threads=$SLURM_NTASKS_PER_NODE -- CoreNumTester.jl




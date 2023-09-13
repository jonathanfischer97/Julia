#!/bin/bash -l

# Define additional arguments
julia_file="SCRIPT_4FixedInitialConditionsCSVMaker.jl"


for df in $@; do
    sbatch julia_rockfish.sh $julia_file $range_number $population_number $df
done


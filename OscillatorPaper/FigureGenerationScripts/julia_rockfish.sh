#!/bin/bash -l
#SBATCH --job-name=NERDSS_restart
#SBATCH --time=72:0:0
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --account=mjohn218_condo
#SBATCH --export=ALL

ml intel-mkl

# Loop over command-line arguments
for dir in "$@"; do
  # Change to the directory
  cd $dir

  # Call the executable with the correct input file
  nerdss.sh -r restart.dat > output.txt &

  # Change back to the original directory
  cd -
done
wait

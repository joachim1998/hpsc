#!/bin/bash

#  module load gcc/4.8.1
#  gcc  -O2  -lm  -fopenmp  -o conducting  main.c

#SBATCH --job-name=conducting
#SBATCH --mail-user=Joachim.Paquay@student.uliege.be
#SBATCH --mail-type=ALL
#SBATCH  --ntasks=16
#SBATCH  --cpus-per-task=16
#        ##################

#SBATCH --mem-per-cpu=2000
#SBATCH --time=30:00

# on NIC4

module load gcc/4.8.1

gcc  -O0  -lm  -std=c99 -fopenmp  -o conducting  main.c

for ii in 1 2 3 4 5 6 7 8 9 10 11; do echo -e "\n*** " $ii; time OMP_NUM_THREADS=$ii ./conducting 1 500 0.23 10000 ; done
#!/bin/bash
# Submission script for NIC4
#SBATCH --job-name=TestRun
#SBATCH --time=10:20:00 # hh:mm:ss
#
#SBATCH --mem-per-cpu=500 # megabytes
#SBATCH --partition=defq
#
#SBATCH --mail-user=Joachim.Paquay@student.uliege.be
#SBATCH --mail-type=ALL
#
#SBATCH --comment=HPCProject-2
#


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_NTASKS

parameter_file=parameters.txt
map_file=sriLanka.dat

mpirun -np $SLURM_NTASKS ./waves $parameter_file $map_file 0

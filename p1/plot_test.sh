#!/bin/bash

#SBATCH --job-name=conducting
#SBATCH --mail-user=Joachim.Paquay@student.uliege.be
#SBATCH --mail-type=ALL
#SBATCH  --ntasks=16
#SBATCH  --cpus-per-task=16
#        ##################

#SBATCH --mem-per-cpu=2000
#SBATCH --time=30:00

# on NIC4

export OMP_NUM_THREADS=10

module load gcc/9.2.0
module load gnuplot

gcc main.c -fopenmp --std=c99 -O3 -lm -o conducting

touch tmpfile.dat
touch tmpfile2.dat

for d in $(seq 0 0.001 0.5)
do
	proba=`./conducting 1 25 $d 3000 | egrep -o "[0-9]{1,}\.[0-9]{2}"`
	echo "$d $proba" >> tmpfile.dat
	proba=`./conducting 1 50 $d 3000 | egrep -o "[0-9]{1,}\.[0-9]{2}"`
	echo "$d $proba" >> tmpfile2.dat
done

gnuplot -e "set term svg enhanced background rgb 'white';
		set output 'plot.svg';
		set yrange[0:100];
		set ylabel 'probability of conduction';
		set xrange[0:$d+0.1];
		set xlabel 'fiber densities';
		set title 'Evolution of the probability of conduction for fibers densities with M = 3000';
		set style line 1 lc rgb 'orange' lt 1;
		set style line 2 lc rgb 'purple' lt 1;
		plot 'tmpfile.dat' with lines ls 1 title 'N=25', 'tmpfile2.dat' with lines ls 2 title 'N=50';";



rm tmpfile.dat
rm tmpfile2.dat
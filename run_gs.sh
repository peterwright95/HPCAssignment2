#!/bin/bash
# 02614 - High-Performance Computing, January 2022
# 
# batch script to run matmult on a decidated server in the hpcintro
# queue
#
# Author: Group 14
#
#BSUB -J jacobi
#BSUB -o mm_batch_%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 15
# uncomment the following line, if you want to assure that your job has
# a whole CPU for itself (shared L3 cache)
# #BSUB -R "span[hosts=1] affinity[socket(1)]"
# choose the server to collect
#BSUB -R "select[model == XeonE5_2650v4]"

module load gcc

# output cpu info
lscpu

EXECUTABLE=./poisson_gs
OUTFILE=gs_serial.dat

# define the grid size here
SIZES="10 20 50 100 200 500 1000"

# define the iteration time
ITER="5000"

# define the tolerance bound
TOL="0.1"

# define the start T
START="1"


# start the collect command with the above settings

for S in $SIZES
do
    $EXECUTABLE $S $ITER $TOL $START >> DATAFILE/$OUTFILE
done




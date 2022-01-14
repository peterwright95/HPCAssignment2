#!/bin/bash
# 02614 - High-Performance Computing, January 2022
# 
# batch script to run matmult on a decidated server in the hpcintro
# queue
#
# Author: Group 14
#
#BSUB -J gsParallel
#BSUB -o mm_batch_%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 60
# uncomment the following line, if you want to assure that your job has
# a whole CPU for itself (shared L3 cache)
# #BSUB -R "span[hosts=1] affinity[socket(1)]"
# choose the server to collect
#BSUB -R "select[model == XeonE5_2650v4]"

module load gcc
module load mpi

# output cpu info
lscpu

EXECUTABLE=./poisson_gs
OUTFILE=gs_parallel.dat

# define the grid size here
SIZES="10 20 50 100 200"

# define the iteration time
ITER="5000"

# define the tolerance bound
TOL="0.1"

# define the start T
START="1"

# define number of threads for parallel computing
THREADS="1 2 4 8 16"


# start the collect command with the above settings
#Clear the file before appending to it

for N in $THREADS
do
    >DATAFILE/gs_parallel_$N.dat
    for S in $SIZES
    do
        OMP_NUM_THREADS=$N $EXECUTABLE $S $ITER $TOL $START >> DATAFILE/gs_parallel_$N.dat
    done
done


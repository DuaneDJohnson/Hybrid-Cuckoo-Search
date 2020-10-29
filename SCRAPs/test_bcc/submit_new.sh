#!/bin/bash
#PBS -u 
#PBS -j 
#PBS -l 
#PBS -l 
#PBS -N 
#############################################################################

export PATH=/opt/gcc-6.1.0/bin:$PATH LD_LIBRARY_PATH=/opt/gcc-6.1.0/lib64:$LD_LIBRARY_PATH

MPIRUN="~/mpirun"
SCRAPS="~/HybridCS.out"

HERE=`pwd`
DIR0=`pwd | sed "s/\//_/g" | sed "s/_home_//g"`
DIR="/scratch/$DIR0"

${MPIRUN}  ${SCRAPS} -p -np xx >& HybridCS.log

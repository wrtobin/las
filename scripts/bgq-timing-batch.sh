#!/bin/bash
#SBATCH --job-name=las-timing
#SBATCH --partition=debug
#SBATCH -t 05
#SBATCH -D /gpfs/u/barn/PASC/PASCtbnw/las/build/test/
#SBATCH --mail-type=ALL
#SBATCH

if [ "$#" -ne  3 ] ; then
  echo "Usage: " $0 " [model] [mesh] "
fi

EXES="petsc_raw petsc_lasops petsc_call petsc_cvirt petsc_virtual"
for EXE in ${EXES} ; do
    echo "running " $EXE
    srun -n 1 -o timing_${EXE}.log $EXE $1 $2
done
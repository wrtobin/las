#!/bin/sh

TIMING_EXES="petsc_raw petsc_lasops petsc_call petsc_cvirt petsc_virtual"
EXE_DIR=$DEVROOT/las/build/test/
for EXE in $TIMING_EXES; do
    mpirun -np 1 $EXE_DIR/$EXE $DEVROOT/problems/simple_cube/models/cube_test.smd $EXE_DIR/cube_24.smb
done

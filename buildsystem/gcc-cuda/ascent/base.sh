#!/bin/bash

export MY_CLUSTER=ascent

module reset

# Load system modules
module load gcc/11.2.0
module load spectrum-mpi/10.4.0.3-20210112

export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)

export MY_PETSC_DIR=$PETSC_DIR
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='jsrun -g 1 -n 1'"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_ENABLE_PYTHON=OFF"

#!/bin/bash

export MY_CLUSTER=ascent

module reset

# Load system modules
module load gcc/10.2.0
module load spectrum-mpi/10.4.0.3-20210112

export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)

export MY_PETSC_DIR=$PETSC_DIR
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='jsrun -g 1 -n 1'"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DLAPACK_LIBRARIES:STRING=/sw/ascent/spack-envs/base/opt/linux-rhel8-ppc64le/gcc-10.2.0/openblas-0.3.17-6te4qwdzetkoyitdryljjstei6jw77gg/lib/libopenblas.so"

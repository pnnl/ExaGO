#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module reset

# System modules
module load PrgEnv-cray
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load cce/15.0.0
module load rocm/5.4.0
module load cray-mpich
module load libfabric

# Consider changing to $(which clang) as for deception
export CC=$(which cc)
export CXX=$(which CC)
export FC=$(which ftn)

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"

#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module reset

# System modules
module load PrgEnv-cray/8.3.3
module load craype-x86-trento
module load craype-accel-amd-gfx90a
module load cce/14.0.2
module load rocm/5.2.0
module load cray-mpich/8.1.23
module load libfabric/1.15.0.0

# Consider changing to $(which clang) as for deception
export CC=$(which craycc)
export CXX=$(which crayCC)
export FC=$(which crayftn)

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"

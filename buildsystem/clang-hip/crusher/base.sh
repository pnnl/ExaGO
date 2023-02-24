#!/bin/bash

export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module reset

# System modules
module load PrgEnv-amd/8.3.3
module load amd/5.2.0
module load cray-mpich/8.1.23

# Consider changing to $(which clang) as for deception
export CC=$(which cc)
export CXX=$(which CC)
export FC=$(which ftn)

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"

#!/bin/bash

export MY_CLUSTER=incline

module purge

# System modules
module load rocm/5.3.0

# Consider changing to $(which clang) as for deception
export CC=$(which clang)
export CXX=$(which clang++)
export FC=$(which flang)

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='srun'"

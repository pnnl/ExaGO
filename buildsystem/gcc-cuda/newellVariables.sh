#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Shared system configuration
source $SRCDIR/buildsystem/gcc-cuda/newell/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/$MY_CLUSTER/modules/dependencies.sh

# Platform specific CMake
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=70 -DEXAGO_ENABLE_IPOPT=ON"


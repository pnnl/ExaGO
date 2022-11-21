#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Shared platform config
source $SRCDIR/buildsystem/gcc-cuda/deception/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/$MY_CLUSTER/modules/dependencies.sh

# Platform specific CMake configuration for building manually

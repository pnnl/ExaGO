#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Shared system configuration
source $SRCDIR/buildsystem/gcc-cuda/newell/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/$MY_CLUSTER/modules/dependencies.sh

# Platform specific CMake

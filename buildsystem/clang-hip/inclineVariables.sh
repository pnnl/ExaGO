#!/bin/bash

SRCDIR=${SRCDIR:-$PWD}

# Platform specific configuration
source $SRCDIR/buildsystem/clang-hip/incline/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/incline/modules/dependencies.sh 

#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Platform specific configuration
source $SRCDIR/buildsystem/clang-hip/crusher/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/crusher/modules/exago-optimized.sh

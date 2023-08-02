#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Platform specific configuration
source $SRCDIR/buildsystem/clang-hip/frontier/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/frontier/modules/exago-optimized.sh


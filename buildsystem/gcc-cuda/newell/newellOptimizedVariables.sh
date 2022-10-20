#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Shared system configuration
source $SRCDIR/buildsystem/gcc-cuda/newell/env.sh

# Spack modules
source $SRCDIR/buildsystem/spack/$MY_CLUSTER/modules/optimized-dependencies.sh


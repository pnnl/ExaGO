#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}
export MY_CLUSTER=ascent

# Spack modules
source $SRCDIR/buildsystem/spack/$MY_CLUSTER/modules/dependencies.sh

# Shared system configuration
source $SRCDIR/buildsystem/gcc-cuda/ascent/base.sh


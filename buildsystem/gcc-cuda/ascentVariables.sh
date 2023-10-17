#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Shared system configuration
source $SRCDIR/buildsystem/gcc-cuda/ascent/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/$MY_CLUSTER/modules/dependencies.sh

module load spectrum-mpi/10.4.0.3-20210112

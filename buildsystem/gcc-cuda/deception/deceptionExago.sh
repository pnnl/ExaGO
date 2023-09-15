#!/bin/bash

export SRCDIR=${SRCDIR:-$PWD}

# Shared platform config
source $SRCDIR/buildsystem/gcc-cuda/deception/base.sh

# Spack modules
source $SRCDIR/buildsystem/spack/$MY_CLUSTER/modules/exago.sh

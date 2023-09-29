#!/bin/bash

# Configure python
module load python/3.9-anaconda3

BASE=/gpfs/wolf/proj-shared/csc359/exago/spack-ci
export SPACK_INSTALL=$BASE/install
export SPACK_CACHE=$BASE/../$(whoami)/spack-cache
export SPACK_MIRROR=/gpfs/wolf/csc359/world-shared/exago/spack-ci/mirror
export SPACK_MODULES=modules
export SPACK_PYTHON=$(which python)
export SPACK_USER_CACHE_PATH=$BASE/../$(whoami)
export SPACK_DISABLE_LOCAL_CONFIG=true

export tempdir=$SPACK_CACHE
export TMP=$SPACK_CACHE
export TMPDIR=$SPACK_CACHE

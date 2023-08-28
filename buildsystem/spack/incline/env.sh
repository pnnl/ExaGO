#!/bin/bash

# Only need to source this file in CI...
. /etc/profile.d/modules.sh
# Clear out any unecessary other modules
module purge

# MPI module is finnicky on incline
modules=$(module list 2>&1)
if echo $modules | grep -q 'openmpi'; then
  module load gcc/8.4.0
  module rm openmpi
fi

# Configure python and other system modules
module load gcc/8.4.0
# ROCm module overrides compiler version used for logic in mpi module...
module load openmpi/4.1.4
# This ROCm module inherently uses 8.4.0
module load rocm/5.3.0
module load python/3.7.0

# Try forcing behaviour of spack compiler
export SPACK_CC=$(which gcc)
export SPACK_CXX=$(which g++)

BASE=/qfs/projects/exasgd/src/ci-incline

export SPACK_INSTALL=$BASE/install
export SPACK_MODULES=modules
export SPACK_CACHE=$BASE/../$(whoami)/cache
export SPACK_MIRROR=$BASE/mirror
export SPACK_PYTHON=$(which python)
export SPACK_DISABLE_LOCAL_CONFIG=1

export tempdir=$SPACK_CACHE
export TMP=$SPACK_CACHE
export TMPDIR=$SPACK_CACHE

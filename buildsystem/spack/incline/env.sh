#!/bin/bash

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
module load cmake/3.21.4

# Try forcing behaviour of spack compiler
export SPACK_CC=$(which gcc)
export SPACK_CXX=$(which g++)


BASE=/vast/projects/exasgd/spack

export SPACK_INSTALL=$BASE/install
export SPACK_MODULES=modules
export SPACK_MIRROR=$BASE/mirror
export SPACK_CACHE=$BASE/cache
export SPACK_PYTHON=/share/apps/python/3.7.0
export SPACK_DISABLE_LOCAL_CONFIG=1


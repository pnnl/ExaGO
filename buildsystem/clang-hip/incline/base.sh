#!/bin/bash

. /etc/profile.d/modules.sh

module purge

# MPI module is finnicky on incline
modules=$(module list 2>&1)
if echo $modules | grep -q 'openmpi'; then
  module load gcc/8.4.0
  module rm openmpi
fi

# System modules
module load gcc/8.4.0
module load openmpi/4.1.4
module load rocm/5.3.0
module load cmake/3.21.4
module load python/3.11.4

# Consider changing to $(which clang) as for deception
export CC=$(which clang)
export CXX=$(which clang++)
export FC=$(which gfortran)

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DAMDGPU_TARGETS='gfx908'"


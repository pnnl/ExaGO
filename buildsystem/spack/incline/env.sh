#!/bin/bash

# Configure python
module load python/3.7.0
module load rocm/5.3.0
module load cmake/3.21.4

BASE=/vast/projects/exasgd/spack/

export SPACK_INSTALL=$BASE/install
export SPACK_MODULES=modules
export SPACK_CACHE=$BASE/cache
export SPACK_PYTHON=/share/apps/python/3.7.0
export SPACK_DISABLE_LOCAL_CONFIG=1

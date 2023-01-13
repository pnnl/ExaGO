#!/bin/bash

module reset

# Configure python
module load cray-python/3.9.12.1 

export SPACK_INSTALL=/gpfs/alpine/proj-shared/csc359/cameron/spack-install
export SPACK_MODULES=test-modules
export SPACK_CACHE=/gpfs/alpine/proj-shared/csc359/$(whoami)/spack-cache
export SPACK_MIRROR=/gpfs/alpine/proj-shared/csc359/$(whoami)/spack-mirror
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1


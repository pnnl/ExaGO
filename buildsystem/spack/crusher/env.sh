#!/bin/bash

# Configure python
module load cray-python/3.9.12.1 

export SPACK_INSTALL=/gpfs/alpine/proj-shared/csc359/cameron/spack-install/amd-03212023
export SPACK_MODULES=test-modules-amd
export SPACK_CACHE=/gpfs/alpine/proj-shared/csc359/$(whoami)/spack-cache
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1


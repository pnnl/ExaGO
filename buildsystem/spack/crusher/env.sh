#!/bin/bash

# Configure python
module load cray-python/3.9.12.1 

BASE=/gpfs/alpine/proj-shared/csc359/cameron
export SPACK_INSTALL=$BASE/spack-install
export SPACK_MODULES=test-modules
export SPACK_CACHE=$BASE/$(whoami)/spack-cache
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1
export SPACK_MIRROR=$BASE/mirror


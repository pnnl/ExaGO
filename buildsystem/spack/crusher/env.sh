#!/bin/bash

# Configure python
module load cray-python/3.9.12.1 

BASE=/lustre/orion/csc359/proj-shared/$(whoami)

export SPACK_INSTALL=/lustre/orion/csc359/proj-shared/nkouk/spack-install
export SPACK_MODULES=modules
export SPACK_CACHE=$BASE/spack-cache
export SPACK_MIRROR=$BASE/spack-mirror
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1


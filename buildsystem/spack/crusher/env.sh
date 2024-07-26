#!/bin/bash

# Configure python
module load cray-python/3.11.5

BASE=/lustre/orion/eng145/world-shared/$(whoami)

export SPACK_INSTALL=/lustre/orion/eng145/world-shared/spack-install
export SPACK_MODULES=modules
export SPACK_CACHE=$BASE/spack-cache
export SPACK_MIRROR=$BASE/spack-mirror
export SPACK_PYTHON=$CRAY_PYTHON_PREFIX
export SPACK_DISABLE_LOCAL_CONFIG=1

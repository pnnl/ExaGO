#!/bin/bash

# Configure python
module load python/3.8-anaconda3

BASE=/gpfs/alpine2/stf006/world-shared/nkouk
export SPACK_INSTALL=$BASE/exago-spack-install
export SPACK_MODULES=summit-modules
export SPACK_CACHE=$BASE/exago-spack-cache
export SPACK_PYTHON=$OLCF_PYTHON_ROOT
export SPACK_DISABLE_LOCAL_CONFIG=1
export SPACK_MIRROR=$BASE/exago-spack-mirror

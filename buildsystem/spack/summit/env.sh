#!/bin/bash

# Configure python
module load python/3.8-anaconda3

BASE=/gpfs/alpine/proj-shared/csc359/cameron
export SPACK_INSTALL=$BASE/spack-install
export SPACK_MODULES=summit-modules
export SPACK_CACHE=$BASE/$(whoami)/spack-cache
export SPACK_PYTHON=$OLCF_PYTHON_ANACONDA3_ROOT
export SPACK_DISABLE_LOCAL_CONFIG=1
export SPACK_MIRROR=$BASE/mirror

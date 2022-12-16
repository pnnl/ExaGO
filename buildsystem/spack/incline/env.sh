#!/bin/bash

# Configure python
module load python/3.7.0

export SPACK_INSTALL=/qfs/projects/exasgd/src/ci-incline
export SPACK_MODULES=test-modules
export SPACK_CACHE=/qfs/projects/exasgd/src/$(whoami)/cache
export SPACK_PYTHON=/share/apps/python/3.7.0
export SPACK_DISABLE_LOCAL_CONFIG=1


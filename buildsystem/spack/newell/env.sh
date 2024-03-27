#!/bin/bash

. /etc/profile.d/modules.sh

# Load system python
module load python/miniconda3.8
. /share/apps/python/miniconda3.8/etc/profile.d/conda.sh

# Load compiler/system modules
module load gcc/8.5.0
module load cuda/11.4
module load openmpi/4.1.4

# Define environment variables for where spack stores key files
# For now, SPACK_INSTALL is the path where everything spack related is installed
# If you want to modify the module install path, edit the spack.yaml manually
BASE=/qfs/projects/earthshot/src/newell-ci
export SPACK_INSTALL=$BASE/spack-install
export SPACK_MODULES=ci-modules
export SPACK_CACHE=$BASE/../$(whoami)/spack-cache
export SPACK_PYTHON=$(which python)
export SPACK_DISABLE_LOCAL_CONFIG=1
export SPACK_MIRROR=$BASE/mirror

export tempdir=$SPACK_CACHE
export TMP=$SPACK_CACHE
export TMPDIR=$SPACK_CACHE

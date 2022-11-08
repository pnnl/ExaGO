#!/bin/bash

source /etc/profile.d/modules.sh
module purge

# Load system python
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh

# Load system modules
module load gcc/10.2.0
module load cuda/11.4
module load openmpi/4.1.0mlx5.0
module load cmake/3.21.4

# Define environment variables for where spack stores key files
# For now, SPACK_INSTALL is the path where everything spack related is installed
# If you want to modify the module install path, edit the spack.yaml manually
export SPACK_INSTALL=/qfs/projects/exasgd/src/ci-deception
export SPACK_MODULES=ci-modules
export SPACK_CACHE=/qfs/projects/exasgd/src/$(whoami)/spack-cache
export SPACK_DISABLE_LOCAL_CONFIG=1
export SPACK_PYTHON=/share/apps/python/miniconda3.9


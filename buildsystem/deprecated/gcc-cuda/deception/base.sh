#!/bin/bash
. /etc/profile.d/modules.sh

module purge

# Load system modules
module load gcc/9.1.0
module load openmpi/4.1.0mlx5.0
module load cuda/11.4

# System Python
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh

export CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)
export MY_CLUSTER=deception

#!/bin/bash
source /etc/profile.d/modules.sh
export MY_CLUSTER=deception

module purge

# Load system modules
module load gcc/10.2.0
module load openmpi/4.1.0mlx5.0
module load cuda/11.4
module load python/miniconda3.8

source /share/apps/python/miniconda3.8/etc/profile.d/conda.sh

export CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)

[ -f $PWD/nvblas.conf ] && rm $PWD/nvblas.conf
cat > $PWD/nvblas.conf <<-EOD
NVBLAS_LOGFILE  nvblas.log
NVBLAS_CPU_BLAS_LIB $OPENBLAS_LIBRARY_DIR/libopenblas.so
NVBLAS_GPU_LIST ALL
NVBLAS_TILE_DIM 2048
NVBLAS_AUTOPIN_MEM_ENABLED
EOD
export NVBLAS_CONFIG_FILE=$PWD/nvblas.conf
echo "Generated $PWD/nvblas.conf"


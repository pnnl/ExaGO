#!/bin/bash

. /etc/profile.d/modules.sh

export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
export MY_CLUSTER=newell

module purge

# Load system modules
module load gcc/8.5.0
module load openmpi/4.1.4
module load cuda/11.4
module load python/miniconda3.8

source /share/apps/python/miniconda3.8/etc/profile.d/conda.sh

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


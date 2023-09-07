#!/bin/bash

export MY_CLUSTER=summit

module reset

# Load system modules
module load gcc/10.2.0
module load spectrum-mpi/10.4.0.3-20210112

export CC=/sw/summit/gcc/10.2.0-2/bin/gcc
export CXX=/sw/summit/gcc/10.2.0-2/bin/g++
export FC=/sw/summit/gcc/10.2.0-2/bin/gfotran

export MY_PETSC_DIR=$PETSC_DIR
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='jsrun -g 1 -n 1'"

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

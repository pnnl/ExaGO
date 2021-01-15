export MY_CLUSTER=ascent
export PROJ_DIR=/gpfs/wolf/proj-shared/csc359
source $PROJ_DIR/src/spack/share/spack/setup-env.sh
module purge
module load cuda/11.0.2
module use $PROJ_DIR/$MY_CLUSTER/Modulefiles/Core
module load exasgd-base
module load gcc-ext/7.4.0
module load spectrum-mpi-ext
module load magma/2.5.3-cuda11
module load raja
module load umpire
module load openblas
module load metis/5.1.0
module load cmake/3.18.2
export MY_RAJA_DIR=$RAJA_ROOT
export MY_UMPIRE_DIR=`spack location -i umpire +cuda`
export MY_UMPIRECPU_DIR=`spack location -i umpire -cuda`
export MY_MAGMA_DIR=$MAGMA_ROOT
export MY_PETSC_DIR=$PROJ_DIR/$MY_CLUSTER/petsc-3.13.5
export MY_IPOPT_DIR=$PROJ_DIR/$MY_CLUSTER/ipopt
export MY_HIOP_DIR=$PROJ_DIR/installs/hiop-gpu-new-interface
export MY_NVCC_ARCH="sm_70"
extraCmakeArgs="\
  $extra_cmake_args \
  -DHIOP_NVCC_ARCH=$MY_NVCC_ARCH \
  -DEXAGO_TEST_WITH_BSUB=ON \
  "
if [[ ! -f $builddir/nvblas.conf ]]; then
  cat > $builddir/nvblas.conf <<-EOD
NVBLAS_LOGFILE  nvblas.log
NVBLAS_CPU_BLAS_LIB  /gpfs/wolf/proj-shared/csc359/ascent/Compiler/gcc-7.4.0/openblas/0.3.10/lib/libopenblas.so
NVBLAS_GPU_LIST ALL
NVBLAS_TILE_DIM 2048
NVBLAS_AUTOPIN_MEM_ENABLED
EOD
fi
export NVBLAS_CONFIG_FILE=$builddir/nvblas.conf
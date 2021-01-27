export MY_CLUSTER=ascent
export PROJ_DIR=/gpfs/wolf/proj-shared/csc359
source $PROJ_DIR/src/spack/share/spack/setup-env.sh
module purge
module load cuda/11.0.2
module use $PROJ_DIR/$MY_CLUSTER/Modulefiles/Core
module load exasgd-base
module load gcc-ext/7.4.0
module load spectrum-mpi-ext
module load openblas
module load cmake/3.18.2

spack env activate exago-v0-99-1 --without-view
spack load petsc umpire raja hiop ipopt magma metis parmetis camp

export MY_RAJA_DIR=`spack location -i raja`
export MY_UMPIRE_DIR=`spack location -i umpire`
export MY_MAGMA_DIR=`spack location -i magma`
export MY_PETSC_DIR=`spack location -i petsc`
export MY_IPOPT_DIR=`spack location -i ipopt`
export MY_HIOP_DIR=`spack location -i hiop`
export MY_NVCC_ARCH="sm_70"
extraCmakeArgs="\
  $extra_cmake_args \
  -DHIOP_NVCC_ARCH=$MY_NVCC_ARCH \
  -DEXAGO_TEST_WITH_BSUB=ON \
  "
if [[ ! -f $builddir/nvblas.conf ]]; then
  cat > $builddir/nvblas.conf <<-EOD
NVBLAS_LOGFILE  nvblas.log
NVBLAS_CPU_BLAS_LIB $(spack location -i openblas)/lib/libopenblas.so
NVBLAS_GPU_LIST ALL
NVBLAS_TILE_DIM 2048
NVBLAS_AUTOPIN_MEM_ENABLED
EOD
  export NVBLAS_CONFIG_FILE=$builddir/nvblas.conf
fi

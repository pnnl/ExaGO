source /etc/profile.d/modules.sh
module purge
export MY_CLUSTER=marianas
export PROJ_DIR=/qfs/projects/exasgd
export APPS_DIR=/share/apps

# When using system packages, we load them directly
module load gcc/7.3.0
module load openmpi/3.1.3
module load cuda/10.2.89
module load cmake/3.15.3

# Configure spack
source $PROJ_DIR/src/spack/share/spack/setup-env.sh
spack env activate exago-ci --without-view

# When loading spack packages, use spack and not modulefiles directly
spack load petsc+mpi umpire+cuda raja+cuda hiop+cuda+mpi ipopt magma openmpi metis parmetis camp+cuda

export MY_PETSC_DIR=`spack location -i petsc+mpi`
export MY_RAJA_DIR=`spack location -i raja+cuda`
export MY_HIOP_DIR=`spack location -i hiop+cuda+mpi`
export MY_UMPIRE_DIR=`spack location -i umpire+cuda`
export MY_IPOPT_DIR=`spack location -i ipopt`
export MY_MAGMA_DIR=`spack location -i magma`
export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf
export extraCmakeArgs="$extraCmakeArgs -DHIOP_NVCC_ARCH=sm_60"

# Variables specified only for build matrix:
export MY_UMPIRECPU_DIR=`spack location -i umpire~cuda`
export MY_PETSC_NOMPI_DIR=`spack location -i petsc~mpi`
export MY_HIOP_NOMPI_DIR=`spack location -i hiop~mpi`

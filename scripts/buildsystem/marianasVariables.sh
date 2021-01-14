source /etc/profile.d/modules.sh
module purge
export MY_CLUSTER=marianas
export PROJ_DIR=/qfs/projects/exasgd
export APPS_DIR=/share/apps
module load gcc/7.3.0
module load cmake/3.15.3
module load openmpi/3.1.3
module load magma/2.5.2_cuda10.2
module load cuda/10.2.89
module load metis/5.1.0
export MY_PETSC_DIR=$PROJ_DIR/$MY_CLUSTER/petsc
export MY_UMPIRE_DIR=$PROJ_DIR/$MY_CLUSTER/umpire
export MY_RAJA_DIR=$PROJ_DIR/$MY_CLUSTER/raja
export MY_HIOP_DIR=$PROJ_DIR/$MY_CLUSTER/hiop-gpu-new-interface
export MY_UMPIRE_DIR=$PROJ_DIR/$MY_CLUSTER/umpire
export MY_IPOPT_DIR=$PROJ_DIR/$MY_CLUSTER/ipopt
export MY_MAGMA_DIR=$APPS_DIR/magma/2.5.2/cuda10.2
export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf

# Variables specified only for build matrix:
export SPACK_INSTALL_DIR=$PROJ_DIR/spack-install/linux-centos7-broadwell/gcc-7.3.0
export MY_UMPIRECPU_DIR=$PROJ_DIR/$MY_CLUSTER/umpire-cpu
export MY_PETSC_NOMPI_DIR=$SPACK_INSTALL_DIR/petsc/3.14.1-mzxn4
export MY_HIOP_NOMPI_DIR=$PROJ_DIR/$MY_CLUSTER/hiop-gpu-nompi-new-interface

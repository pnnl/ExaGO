source /etc/profile.d/modules.sh
module purge
export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
export MY_CLUSTER=newell
export PROJ_DIR=/qfs/projects/exasgd
export APPS_DIR=/share/apps
module load gcc/7.4.0
module load cmake/3.16.4
module load openmpi/3.1.5
module load cuda/10.2

# Configure spack
source $PROJ_DIR/src/spack/share/spack/setup-env.sh
spack env activate exago-v0-99-1-newell --without-view

# When loading spack packages, use spack and not modulefiles directly
spack load petsc umpire raja hiop ipopt magma openmpi metis parmetis camp

export MY_PETSC_DIR=`spack location -i petsc+mpi`
export MY_UMPIRE_DIR=`spack location -i umpire+cuda`
export MY_RAJA_DIR=`spack location -i raja+cuda`
export MY_HIOP_DIR=`spack location -i hiop+cuda+shared`
export MY_IPOPT_DIR=`spack location -i ipopt`
export MY_MAGMA_DIR=`spack location -i magma`
export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf

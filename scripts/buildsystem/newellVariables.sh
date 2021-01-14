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
module load magma/2.5.2_cuda10.2
module load metis/5.1.0
module load cuda/10.2
export MY_PETSC_DIR=$PROJ_DIR/$MY_CLUSTER/petsc
export MY_UMPIRE_DIR=$PROJ_DIR/$MY_CLUSTER/umpire
export MY_RAJA_DIR=$PROJ_DIR/$MY_CLUSTER/raja
export MY_HIOP_DIR=$PROJ_DIR/$MY_CLUSTER/hiop-gpu-new-interface
export MY_UMPIRE_DIR=$PROJ_DIR/$MY_CLUSTER/umpire
export MY_UMPIRECPU_DIR=$PROJ_DIR/$MY_CLUSTER/umpire-cpu
export MY_IPOPT_DIR=$PROJ_DIR/$MY_CLUSTER/ipopt
export MY_MAGMA_DIR=$APPS_DIR/magma/2.5.2/cuda10.2
export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf

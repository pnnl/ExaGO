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
source ~/software/ExaSGD_Spack/spack/share/spack/setup-env.sh

# Dirty workaround to fix permissions errors
# see https://github.com/spack/spack/issues/17407
ls ~/software/ExaSGD_Spack/spack/var/spack/environments/*

#spack env activate exago-v0-99-2-hiop-v0-3-99-2-newell
spack env activate exago-hiop-pridecomp

# Petsc is the only dependency that needs an explicit path
export MY_PETSC_DIR=`spack location -i petsc`
#export MY_HIOP_DIR=/people/abhy245/software/hiop-github/hiop/_dist-Debug/
export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf

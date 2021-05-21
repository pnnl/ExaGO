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
ls $PROJ_DIR/src/spack/var/spack/*
#spack env activate exago-v0-99-2-hiop-v0-3-99-2-marianas
spack env activate exago-v1-0-0-deps-marianas

export MY_PETSC_DIR=`spack location -i petsc`
export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DHIOP_NVCC_ARCH=sm_60"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DLAPACK_LIBRARIES=$(spack location -i openblas)"

source /etc/profile.d/modules.sh
module purge
export MY_CLUSTER=newell
export PROJ_DIR=/qfs/projects/exasgd
export APPS_DIR=/share/apps

# Configure spack
source $PROJ_DIR/src/spack/share/spack/setup-env.sh

module use -a $PROJ_DIR/src/spack/share/spack/modules/linux-rhel7-power9le
module load llvm/13.0.0
module load cmake-3.18.4-clang-13.0.0-s7aylhf

# Dirty workaround to fix permissions errors
# see https://github.com/spack/spack/issues/17407
ls $PROJ_DIR/src/spack/var/spack/environments/*

spack env activate exago-v1-0-0-hiop-v0-4-0-clang

# Petsc is the only dependency that needs an explicit path
export MY_PETSC_DIR=`spack location -i petsc`

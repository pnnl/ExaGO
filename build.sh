#!/bin/bash

cleanup() {
  echo
  echo Exit code $1 caught in build script.
  echo
  if [[ "$1" == "0" ]]; then
    echo BUILD_STATUS:0
  else
    echo
    echo Failure found on line $2 in build script.
    echo
    echo BUILD_STATUS:1
  fi
}

trap 'cleanup $? $LINENO' EXIT

set -xv
make_args="-j 8"
export CXXFLAGS='-I/share/apps/magma/2.5.2/cuda10.2/include -I/share/apps/cuda/10.2/include/'
export OMPI_MCA_btl="^vader,tcp,openib,uct"

if [[ ! -v MY_CLUSTER ]]
then
  export MY_CLUSTER=`uname -n | sed -e 's/[0-9]//g' -e 's/\..*//'`
fi

ulimit -s unlimited || echo 'Could not set stack size to unlimited.'
ulimit -l unlimited || echo 'Could not set max locked memory to unlimited.'

. /etc/profile.d/modules.sh
module purge

case $MY_CLUSTER in
newell*)
  export OMP_CANCELLATION=true
  export OMP_PROC_BIND=true
  export OMPI_MCA_pml="ucx"
  export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
  export MY_CLUSTER=newell
  module load gcc/7.4.0
  module load cmake/3.16.4
  module load openmpi/3.1.5
  module load magma/2.5.2_cuda10.2
  module load metis/5.1.0
  module load cuda/10.2
  ;;
dl|shared|marianas)
  export MY_CLUSTER=marianas
  module load gcc/7.3.0
  module load cmake/3.15.3
  module load openmpi/3.1.3
  module load magma/2.5.2_cuda10.2
  module load metis/5.1.0
  module load cuda/10.2.89
  ;;
*)
  echo
  echo Script not configured for this cluster.
  echo
  exit 1
  ;;
esac

builddir=$(pwd)/build
installdir=$(pwd)/install
petsc_dir=/qfs/projects/exasgd/$MY_CLUSTER/petsc
hiop_dir=/qfs/projects/exasgd/$MY_CLUSTER/hiop
raja_dir=/qfs/projects/exasgd/$MY_CLUSTER/raja
umpire_dir=/qfs/projects/exasgd/$MY_CLUSTER/umpire
magma_dir=/share/apps/magma/2.5.2/cuda10.2/
metis_dir=/share/apps/metis/5.1.0/
cuda_dir=/share/apps/cuda/10.2/
export NVBLAS_CONFIG_FILE=/qfs/projects/exasgd/$MY_CLUSTER/nvblas.conf
export LD_LIBRARY_PATH="$magma_dir/lib:$hiop_dir/lib:$LD_LIBRARY_PATH"
export CTEST_OUTPUT_ON_FAILURE=1

module list

declare -a builds=(
  " \
  -DCMAKE_INSTALL_PREFIX=$installdir/ \
  -DCMAKE_BUILD_TYPE=Debug \
  -DEXAGO_ENABLE_GPU=ON \
  -DEXAGO_ENABLE_HIOP=ON \
  -DEXAGO_ENABLE_IPOPT=ON \
  -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_PETSC=ON \
  -DEXAGO_RUN_TESTS=ON \
  -DEXAGO_ENABLE_RAJA=ON \
  -DRAJA_DIR=$raja_dir \
  -Dumpire_DIR=$umpire_dir \
  -DHIOP_DIR=$hiop_dir \
  -DIPOPT_DIR=/qfs/projects/exasgd/$MY_CLUSTER/ipopt \
  -DMAGMA_DIR=$magma_dir \
  -DPETSC_DIR=$petsc_dir"
)


#  NOTE: The following is required when running from Gitlab CI via slurm job
base_path=`dirname $0`
if [ -z "$SLURM_SUBMIT_DIR" ]; then
    cd $base_path          || exit 1
fi

for ((i=0; i<${#builds[@]}; i++))
do
  echo
  echo Build $((i+1)) / ${#builds[@]}
  echo
  args=${builds[i]}
  echo Building with args $args

  [ -d $builddir ] && rm -rf $builddir
  mkdir -p $builddir

  [ -d $installdir ] && rm -rf $installdir
  mkdir -p $installdir

  mkdir $builddir/datafiles
  for f in case118.m case9mod.m case_ACTIVSg200.m
  do
    [ -f ./datafiles/$f ] || {
      echo Could not find needed data files.
      exit 1
    }
    cp ./datafiles/$f $builddir/datafiles/$f
  done

  pushd $builddir

  echo
  echo Configuring
  echo
  cmake $args .. || exit 1

  echo
  echo Building
  echo
  make $make_args || exit 1

  echo
  echo Installing
  echo
  make install || exit 1

  echo
  echo Testing
  echo
  ctest -VV || exit 1
  popd

  echo
  echo Build $((i+1)) / ${#builds[@]} successful
  echo
done

exit 0

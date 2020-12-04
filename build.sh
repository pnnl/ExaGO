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
ctest_args=" -VV "
export CXXFLAGS='-I/share/apps/magma/2.5.2/cuda10.2/include -I/share/apps/cuda/10.2/include/'
export OMPI_MCA_btl="^vader,tcp,openib,uct"
export BUILD=1
export TEST=1

while [[ $# -gt 0 ]]; do
  case $1 in
  --short)
    echo "Only running representative subset of unit tests."
    ctest_args="$ctest_args -R UNIT_TEST"
    shift
    ;;
  --build-only)
    export TEST=0 BUILD=1
    shift
    ;;
  --test-only)
    export TEST=1 BUILD=0
    shift
    ;;
  *)
    echo "Argument $1 not recognized."
    exit 1
    ;;
  esac
done

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
  export MY_HIOP_DIR=$PROJ_DIR/$MY_CLUSTER/hiop
  export MY_UMPIRE_DIR=$PROJ_DIR/$MY_CLUSTER/umpire
  export MY_IPOPT_DIR=$PROJ_DIR/$MY_CLUSTER/ipopt
  export MY_MAGMA_DIR=$APPS_DIR/magma/2.5.2/cuda10.2
  export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf
  ;;
dl|shared|marianas)
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
  export MY_HIOP_DIR=$PROJ_DIR/$MY_CLUSTER/hiop
  export MY_UMPIRE_DIR=$PROJ_DIR/$MY_CLUSTER/umpire
  export MY_IPOPT_DIR=$PROJ_DIR/$MY_CLUSTER/ipopt
  export MY_MAGMA_DIR=$APPS_DIR/magma/2.5.2/cuda10.2
  export NVBLAS_CONFIG_FILE=$PROJ_DIR/$MY_CLUSTER/nvblas.conf
  ;;
ascent)
  export MY_CLUSTER=ascent
  export PROJ_DIR=/gpfs/wolf/proj-shared/csc359
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
  export MY_UMPIRE_DIR=$UMPIRE_ROOT
  export MY_MAGMA_DIR=$MAGMA_ROOT
  export MY_PETSC_DIR=$PROJ_DIR/$MY_CLUSTER/petsc-3.13.5
  export MY_IPOPT_DIR=$PROJ_DIR/$MY_CLUSTER/ipopt
  export MY_HIOP_DIR=$PROJ_DIR/installs/hiop
  export MY_NVCC_ARCH="sm_70"
  extra_cmake_args="\
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
export LD_LIBRARY_PATH="$magma_dir/lib:$hiop_dir/lib:$LD_LIBRARY_PATH"
export CTEST_OUTPUT_ON_FAILURE=1

module list

cmake_args=" \
  -DCMAKE_INSTALL_PREFIX=$installdir/ \
  -DCMAKE_BUILD_TYPE=Debug \
  -DEXAGO_ENABLE_GPU=ON \
  -DEXAGO_ENABLE_HIOP=ON \
  -DEXAGO_ENABLE_IPOPT=ON \
  -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_PETSC=ON \
  -DEXAGO_RUN_TESTS=ON \
  -DEXAGO_ENABLE_RAJA=ON \
  -DEXAGO_ENABLE_IPOPT=ON \
  -DIPOPT_DIR=$MY_IPOPT_DIR \
  -DRAJA_DIR=$MY_RAJA_DIR \
  -Dumpire_DIR=$MY_UMPIRE_DIR \
  -DHIOP_DIR=$MY_HIOP_DIR \
  -DMAGMA_DIR=$MY_MAGMA_DIR \
  -DPETSC_DIR=$MY_PETSC_DIR \
  $extra_cmake_args"

#  NOTE: The following is required when running from Gitlab CI via slurm job
base_path=`dirname $0`
if [ -z "$SLURM_SUBMIT_DIR" ]; then
    cd $base_path          || exit 1
fi

if [[ $BUILD -eq 1 ]]; then
  echo Building with args $cmake_args

  [ -d $builddir ] && rm -rf $builddir
  mkdir -p $builddir

  [ -d $installdir ] && rm -rf $installdir
  mkdir -p $installdir

  mkdir $builddir/datafiles
  for f in case118.m case9/case9mod.m case_ACTIVSg200.m
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
  cmake $cmake_args .. || exit 1

  echo
  echo Building
  echo
  make $make_args || exit 1

  echo
  echo Installing
  echo
  make install || exit 1
  popd
fi

if [[ $TEST -eq 1 ]]; then
  pushd $builddir
  echo
  echo Testing
  echo
  ctest $ctest_args || exit 1
  popd
fi

exit 0

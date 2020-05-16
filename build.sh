#!/bin/bash

make_args="-j 8"
builddir=$(realpath $(pwd))/build
petsc_dir=/qfs/projects/exasgd/newell/petsc
hiop_dir=/qfs/projects/exasgd/newell/hiop_shared_gpu
magma_dir=/share/apps/magma/2.5.2/cuda10.2/
metis_dir=/share/apps/metis/5.1.0/
cuda_dir=/share/apps/cuda/10.2/
export NVBLAS_CONFIG_FILE=/qfs/projects/exasgd/newell/nvblas.conf
export LD_LIBRARY_PATH="$magma_dir/lib:$hiop_dir/lib:$LD_LIBRARY_PATH"

[[ $(uname -n) =~ newell* ]] || echo 'Script configured for newell. May not build correctly.'

module purge
module load gcc/7.4.0
module load cmake/3.16.4
module load openmpi/3.1.5
module load cuda/10.2
module load magma/2.5.2_cuda10.2
module load metis/5.1.0

declare -a builds=(
  "-DSCOPFLOW_ENABLE_HIOP=ON \
    -DHIOP_DIR=$hiop_dir \
    -DCMAKE_INSTALL_PREFIX=$builddir -DCMAKE_BUILD_TYPE=Debug \
    -DPETSC_DIR=$petsc_dir \
    -DCMAKE_INSTALL_PREFIX=$builddir \
    -DMAGMA_DIR=$magma_dir \
    -DSCOPFLOW_RUN_TESTS=ON"
)

for i in $(seq 0 0)
do
  args=${builds[i]}
  echo Building with args $args
  [ -d $builddir ] && rm -rf $builddir
  mkdir -p $builddir
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
  ctest || exit 1

  popd $builddir

  echo
  echo Build successful
  echo
done

echo
echo All builds successul
echo
exit 0

#!/bin/bash

make_args="-j 8"
builddir=$(realpath $(pwd))/build
installdir=$(realpath $(pwd))/install
petsc_dir=/qfs/projects/exasgd/newell/petsc
hiop_dir=/qfs/projects/exasgd/newell/hiop
magma_dir=/share/apps/magma/2.5.2/cuda10.2/
metis_dir=/share/apps/metis/5.1.0/
cuda_dir=/share/apps/cuda/10.2/
export NVBLAS_CONFIG_FILE=/qfs/projects/exasgd/newell/nvblas.conf
export LD_LIBRARY_PATH="$magma_dir/lib:$hiop_dir/lib:$LD_LIBRARY_PATH"

[[ $(uname -n) =~ newell* ]] || echo 'Script configured for newell. May not build correctly.'

# Needed to supress/resolve openmpi warnings
export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl="^vader,tcp,openib,uct"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1

source /etc/profile.d/modules.sh
module purge
module load gcc/7.4.0
module load cmake/3.16.4
module load openmpi/3.1.5
module load magma/2.5.2_cuda10.2
module load metis/5.1.0
module load cuda/10.2

module list

export CXXFLAGS='-I/share/apps/magma/2.5.2/cuda10.2/include -I/share/apps/cuda/10.2/include/'

declare -a builds=(
  "-DSCOPFLOW_ENABLE_HIOP=ON \
  -DHIOP_DIR=/qfs/projects/exasgd/newell/hiop \
  -DCMAKE_INSTALL_PREFIX=$installdir/ \
  -DCMAKE_BUILD_TYPE=Debug \
  -DSCOPFLOW_ENABLE_IPOPT=OFF \
  -DIPOPT_DIR=/qfs/projects/exasgd/newell/ipopt \
  -DMAGMA_DIR=/share/apps/magma/2.5.2/cuda10.2/ \
  -DSCOPFLOW_RUN_TESTS=ON \
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

echo
echo All builds successul
echo
exit 0

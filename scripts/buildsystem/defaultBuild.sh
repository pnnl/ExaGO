#
# Builds default configuration of ExaGO
#

for requiredVar in MY_IPOPT_DIR MY_RAJA_DIR MY_UMPIRE_DIR MY_HIOP_DIR MY_MAGMA_DIR MY_PETSC_DIR
do
  if [[ ! -v $requiredVar ]]; then
    echo "Required CMake variable $requiredVar was not set!"
    echo "script defaultBuild.sh may not function as you intend."
  fi
done

cmakeArgs=" \
  -DCMAKE_INSTALL_PREFIX=$installdir/ \
  -DCMAKE_BUILD_TYPE=Debug \
  -DEXAGO_HAVE_GPU=ON \
  -DEXAGO_HAVE_HIOP=ON \
  -DEXAGO_HAVE_IPOPT=ON \
  -DEXAGO_HAVE_MPI=ON \
  -DEXAGO_HAVE_PETSC=ON \
  -DEXAGO_RUN_TESTS=ON \
  -DEXAGO_HAVE_RAJA=ON \
  -DEXAGO_HAVE_IPOPT=ON \
  -DIPOPT_DIR=$MY_IPOPT_DIR \
  -DRAJA_DIR=$MY_RAJA_DIR \
  -Dumpire_DIR=$MY_UMPIRE_DIR \
  -DHIOP_DIR=$MY_HIOP_DIR \
  -DMAGMA_DIR=$MY_MAGMA_DIR \
  -DPETSC_DIR=$MY_PETSC_DIR \
  $EXTRA_CMAKE_ARGS"

for requiredVar in builddir installdir makeArgs ctestArgs
do
  if [[ ! -v $requiredVar ]]
  then
    echo "Required variable $requiredVar was not set!"
    echo "script defaultBuild.sh may not function as you intend."
  fi
done
  
function defaultBuild {
  
  if [[ $BUILD -eq 1 ]]; then
    echo Building with args $cmakeArgs
  
    [ -d $builddir ] && rm -rf $builddir
    mkdir -p $builddir
  
    [ -d $installdir ] && rm -rf $installdir
    mkdir -p $installdir
  
    mkdir $builddir/datafiles
    for f in case118.m case9/case9mod.m case_ACTIVSg200.m
    do
      [ -f ./datafiles/$f ] || {
        echo Could not find needed data files.
        return 1
      }
      cp ./datafiles/$f $builddir/datafiles/$f
    done
  
    pushd $builddir
  
    echo
    echo Configuring
    echo
    cmake $cmakeArgs .. || return 1
  
    echo
    echo Building
    echo
    make $makeArgs || return 1
  
    echo
    echo Installing
    echo
    make install || return 1
    popd
  fi
  
  if [[ $TEST -eq 1 ]]; then
    pushd $builddir
    echo
    echo Testing
    echo
    ctest $ctestArgs || return 1
    popd
  fi

  return 0
}

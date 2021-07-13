#
# Builds default configuration of ExaGO
#

function doBuild {

  for requiredVar in SRCDIR BUILDDIR INSTALLDIR MAKEARGS CTESTARGS
  do
    if [[ ! -v $requiredVar ]]
    then
      echo "Required variable $requiredVar was not set!"
      echo "script defaultBuild.sh may not function as you intend."
    fi
  done
  
  if [[ $BUILD -eq 1 ]]; then
    echo Building with args $CMAKEARGS
  
    [ -d $BUILDDIR ] && rm -rf $BUILDDIR
    mkdir -p $BUILDDIR
  
    [ -d $INSTALLDIR ] && rm -rf $INSTALLDIR
    mkdir -p $INSTALLDIR
  
    pushd $BUILDDIR
  
    echo
    echo Configuring
    echo
    cmake \
      -C $SRCDIR/buildsystem/gcc-cuda/cache.cmake \
      $EXTRA_CMAKE_ARGS \
      .. || return 1
  
    echo
    echo Building
    echo
    make $MAKEARGS || return 1
  
    echo
    echo Installing
    echo
    make install || return 1
    popd
  fi
  
  if [[ $TEST -eq 1 ]]; then
    pushd $BUILDDIR
    echo
    echo Testing
    echo
    ctest $CTESTARGS || return 1
    popd
  fi

  return 0
}

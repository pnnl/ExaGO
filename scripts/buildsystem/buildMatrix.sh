#
# Full matrix script
#
# To redirect the logs for the build matrix, set variable LOGFILE to a full path.

for requiredVar in MY_IPOPT_DIR MY_RAJA_DIR MY_UMPIRE_DIR MY_HIOP_DIR MY_MAGMA_DIR MY_PETSC_DIR
do
  if [[ ! -v $requiredVar ]]; then
    echo "Required CMake variable $required_var was not set!"
    echo "script defaultBuild.sh may not function as you intend."
  fi
done

# Configure logging
export logFile=${LOGFILE:-"$srcdir/buildmatrix-log.csv"}
if [[ -f $logFile ]]; then
  rm $logFile
fi
touch $logFile
export logCols="index;build_status;cmake_args"
echo "$logCols" >> $logFile
export index=0

# Error handling to make reporting easier
errorMessages=("SUCCESS" "CMAKE ERROR" "BUILD ERROR" "INSTALL ERROR" "TEST ERROR" "OTHER ERROR")
export success=0 cmakeError=1 buildError=2 installError=3 testError=4 otherError=5

function logRun {
  local retcode=${1:-"-1"}
  local cmakeArgs=${2:-"-1"}
  if [[ $retcode -eq -1 ]] || [[ $cmakeArgs -eq -1 ]]; then
    echo 'Must pass two arguments to function `logRun`'
    exit 1
  fi
  echo "Run number $index"
  local buildStatus="${errorMessages[retcode]}"
  local cmakeArgs="$(echo "$cmakeArgs" | sed 's/\s\+/ /g')"
  echo "$index;$buildStatus;$cmakeArgs" | tee -a $logFile
  echo 'Full logs so far:'
  cat $logFile
  let index++
}
# End logging code

# If BUILD_MATRIX_PARALLEL == 1, a single combination from the build matrix
# will be run. Otherwise, the entire matrix will be run.
function buildMatrix {
  baseCmakeArgs=" \
    -DCMAKE_INSTALL_PREFIX=$installdir/ \
    -DCMAKE_BUILD_TYPE=Debug \
    -DEXAGO_RUN_TESTS=ON \
    -DEXAGO_ENABLE_PETSC=ON \
    $extra_cmake_args"

  # Setting Umpire dir will not be a problem if umpire is not enabled.
  # We just need to enable the correct umpire installation.
  gpuOptions=("-DEXAGO_ENABLE_GPU=ON -DMAGMA_DIR=$MY_MAGMA_DIR -Dumpire_DIR=$MY_UMPIRE_DIR"
              "-DEXAGO_ENABLE_GPU=OFF -Dumpire_DIR=$MY_UMPIRECPU_DIR")
  
  ipoptOptions=("-DEXAGO_ENABLE_IPOPT=ON -DIPOPT_DIR=$MY_IPOPT_DIR"
                "-DEXAGO_ENABLE_IPOPT=OFF")

  rajaOptions=("-DEXAGO_ENABLE_RAJA=ON -DRAJA_DIR=$MY_RAJA_DIR"
               "-DEXAGO_ENABLE_RAJA=OFF")

  hiopOptions=("-DEXAGO_ENABLE_HIOP=ON -DHIOP_DIR=$MY_HIOP_DIR"
               "-DEXAGO_ENABLE_HIOP=OFF")

  # We have MPI and PETSc variables in this array, but their ENABLE_ cmake
  # variables are toggled elsewhere, so this should not be a problem.
  mpiOptions=("-DEXAGO_ENABLE_MPI=ON -DHIOP_DIR=$MY_HIOP_DIR -DPETSC_DIR=$MY_PETSC_DIR" 
              "-DEXAGO_ENABLE_MPI=OFF -DHIOP_DIR=$MY_HIOP_NOMPI_DIR -DPETSC_DIR=$MY_PETSC_NOMPI_DIR")

  if [[ $BUILD_MATRIX_PARALLEL -eq 1 ]]
  then
    echo 'Running parallel build matrix'
    echo "CI_GPUOPT  : $CI_GPUOPT"
    echo "CI_IPOPTOPT: $CI_IPOPTOPT"
    echo "CI_RAJAOPT : $CI_RAJAOPT"
    echo "CI_HIOPOPT : $CI_HIOPOPT"
    echo "CI_MPIOPT  : $CI_MPIOPT"
    gpuOpt="${gpuOptions[CI_GPUOPT]}"
    ipoptOpt="${ipoptOptions[CI_IPOPTOPT]}"
    rajaOpt="${rajaOptions[CI_RAJAOPT]}"
    hiopOpt="${hiopOptions[CI_HIOPOPT]}"
    mpiOpt="${mpiOptions[CI_MPIOPT]}"
    cmakeArgs="$baseCmakeArgs $gpuOpt $ipoptOpt $rajaOpt $hiopOpt $mpiOpt"
    export cmakeArgs
    runBuild
    return $?
  else
    local fail=0
    for gpuOpt in "${gpuOptions[@]}"; do
      for ipoptOpt in "${ipoptOptions[@]}"; do
        for rajaOpt in "${rajaOptions[@]}"; do
          for hiopOpt in "${hiopOptions[@]}"; do
            for mpiOpt in "${mpiOptions[@]}"; do
              cmakeArgs="$baseCmakeArgs $gpuOpt $ipoptOpt $rajaOpt $hiopOpt $mpiOpt"
              export cmakeArgs
              runBuild
              ret=$?
              logRun $ret "$cmakeArgs"
              let fail+=$ret
            done
          done
        done
      done
    done
    return $fail
  fi
}

function runBuild {
  if [[ ! -v cmakeArgs ]]; then
    echo "Variable 'cmakeArgs' not set!"
    return $otherError
  fi

  echo Building with args $cmakeArgs

  [ -d $builddir ] && rm -rf $builddir
  mkdir -p $builddir
  [ -d $installdir ] && rm -rf $installdir
  mkdir -p $installdir
  mkdir $builddir/datafiles
  for f in case118.m case9/case9mod.m case_ACTIVSg200.m
  do
    [ -f $srcdir/datafiles/$f ] || {
      echo Could not find needed data files.
      return $otherError
    }
    cp $srcdir/datafiles/$f $builddir/datafiles/$f
  done

  pushd $builddir

  echo
  echo Configuring
  echo
  cmake $cmakeArgs $extraCmakeArgs .. || return $cmakeError

  echo
  echo Building
  echo
  make $makeArgs || return $buildError

  echo
  echo Installing
  echo
  make install || return $installError

  popd

  pushd $builddir
  echo
  echo Testing
  echo
  ctest $ctestArgs || return $testError
  popd

  return $success
}

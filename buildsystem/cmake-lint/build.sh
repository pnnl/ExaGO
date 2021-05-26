#!/usr/bin/env bash

cmakeFormat=${CMAKEFORMAT_EXECUTABLE:-'cmake-format'}
which $cmakeFormat
if [[ $? -eq 1 ]]; then
  echo
  echo No cmake-format script was found!
  echo cmake-format may be installed via pip with 'pip install cmake-format'.
  echo Please install cmake-format and add its directory to your path and call
  echo this again.
  echo
  exit 1
fi
cmakeFormat="$cmakeFormat --log-level info"

export srcdir=${CMAKE_SOURCE_DIR:-$PWD}

set -x

function doBuild {
  files="$srcdir/CMakeLists.txt"
  for d in buildsystem/cmake src tests applications; do
    files="$files
$(find $srcdir/$d -name '*.cmake' -or -name CMakeLists.txt)"
  done

  local ret=0
  while read f; do
    $cmakeFormat --check $f
    ((ret=$ret + $?))
  done <<< "$files"

  return $ret
}

if [ $1 = "--run" ]; then
  doBuild
  ret=$?

  if [[ $ret -ne 0 ]]; then
    echo 'Linting failed!'
  fi
  exit $?
fi

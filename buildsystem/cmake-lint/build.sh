#!/usr/bin/env bash

cmakeFormat="cmake-format --log-level info"

cd ..
export srcdir=${srcdir:-$PWD}
export findCmakeFiles=`find $srcdir  \
  -type f -name '*.cmake' -o \
  -name CMakeLists.txt \
  | egrep -v "/toml11/|/pybind11/|/build/|/install/"`

set -x

function doBuild {
  local ret=0
  for f in $findCmakeFiles; do
    $cmakeFormat --check $f
    ret=$(($ret + $?))
  done
  return $ret
}

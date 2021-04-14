#!/usr/bin/env bash

cmakeFormat="cmake-format --log-level info"

export srcdir=${srcdir:-$PWD}
export findCmakeFiles=`find $srcdir -name '*.cmake' -or -name CMakeLists.txt`

set -x

function doBuild {
  local ret=0
  for f in $findCmakeFiles; do
    $cmakeFormat --check $f
    ret=$(($ret + $?))
  done
  return $ret
}

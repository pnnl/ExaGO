#!/usr/bin/env bash

cmakeFormat="cmake-format --log-level info"

cd ..
export srcdir=${srcdir:-$PWD}
export findCmakeFiles=`find $srcdir -type f -name '*.cmake' -o \
  -name CMakeLists.txt \
  -a -not -path '*/tests/toml11/*'\
  -a -not -path '*/build/*'\
  -a -not -path '*/install/*'`

set -x

function doBuild {
  local ret=0
  for f in $findCmakeFiles; do
    $cmakeFormat --check $f
    ret=$(($ret + $?))
  done
  return $ret
}

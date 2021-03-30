#!/usr/bin/env bash

module load python
set -x
pip install cmake-format --user
cmakeFormat="~/.local/bin/cmake-format --log-level info"
srcdir=${srcdir:-$PWD}

export findCmakeFiles="find $srcdir -name *.cmake \
  -or -name CMakeLists.txt"

if [[ $1 = "--help" ]]; then
  echo '------------------------------------------------------------------------'
  echo
  echo 'To format all cmake files in ExaGO, run the following command:'
  echo
  echo "$findCmakeFiles -exec $cmakeFormat --in-place {} ';'"
  echo
  echo '------------------------------------------------------------------------'
  exit 0
fi

function checkCmakeFormat {
  local ret=0
  for f in `$findCmakeFiles`; do
    ~/.local/bin/cmake-format --log-level info --check $f
    ret=$(($ret + $?))
  done
  return $ret
}

function applyCmakeFormat {
  for f in `$findCmakeFiles`; do
    ~/.local/bin/cmake-format --log-level info --in-place $f
  done
}

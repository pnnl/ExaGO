#!/usr/bin/env bash
(return 0 2>/dev/null) && isSourced=1 || isSourced=0

module load python
set -x
pip install cmake-format --user
export PATH="$PATH:~/.local/bin"
export srcdir=${$srcdir:-$PWD}
cmakeFormat="cmake-format --log-level info"
srcdir=${srcdir:-$PWD}

export findCmakeFiles=`find $srcdir -name '*.cmake' -or -name CMakeLists.txt`

if [[ $1 = "--help" ]]; then
  set +x
  echo '------------------------------------------------------------------------'
  echo
  echo 'To format all cmake files in ExaGO, run the following command:'
  echo
  echo "for f in '$findCmakeFiles'; do $cmakeFormat --in-place \$f; done"
  echo
  echo '------------------------------------------------------------------------'
  exit 0
fi

function checkCmakeFormat {
  local ret=0
  for f in $findCmakeFiles; do
    $cmakeFormat --check $f
    ret=$(($ret + $?))
  done
  return $ret
}

function applyCmakeFormat {
  for f in $findCmakeFiles; do
    $cmakeFormat --in-place $f
  done
}

if [[ isSourced -eq 0 ]]; then
  checkCmakeFormat
  exit $?
fi

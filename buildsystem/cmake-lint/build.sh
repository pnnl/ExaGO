#!/usr/bin/env bash
(return 0 2>/dev/null) && isSourced=1 || isSourced=0

cmakeFormat="cmake-format --log-level info"

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

function doBuild {
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

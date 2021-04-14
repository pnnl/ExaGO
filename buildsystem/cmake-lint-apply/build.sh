
set -x
which cmake-format
if [[ $? -eq 1 ]]; then
  echo
  echo No cmake-format script was found!
  echo cmake-format may be installed via pip with 'pip install cmake-format'.
  echo Please install cmake-format and add its directory to your path and call
  echo this again.
  echo
  exit 1
fi

export cmakeFormat="cmake-format --log-level info"
export srcdir=${srcdir:-$PWD}
export findCmakeFiles=`find $srcdir -name '*.cmake' -or -name CMakeLists.txt`

function doBuild {
  for f in $findCmakeFiles; do
    $cmakeFormat --in-place $f
  done
}

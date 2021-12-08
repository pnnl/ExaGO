
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

cd ..
export srcdir=${srcdir:-$PWD}
function doBuild {
  $srcdir/scripts/cmake-format.pl -vi
}

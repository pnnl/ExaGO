#!/bin/bash

source buildsystem/spack/load_spack.sh &&
SPACK_MIRROR="${SPACK_MIRROR:?SPACK_MIRROR is unset. Please use the load_spack script first.}" &&
spack -e $SPACKENV develop --path=$(pwd) exago@develop &&
spack -e $SPACKENV bootstrap now &&
spack -e $SPACKENV concretize -f &&
spack -e $SPACKENV mirror create -a --exclude-specs exago@develop --directory $SPACK_MIRROR &&
spack -e $SPACKENV mirror add local file://$SPACK_MIRROR &&
spack -e $SPACKENV mirror list

res=$?

chmod -R ugo+wrx $SPACK_MIRROR/* &

exit $res

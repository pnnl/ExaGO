#!/bin/bash

source buildsystem/spack/load_spack.sh &&
SPACK_MIRROR="${SPACK_MIRROR:?SPACK_MIRROR is unset. Please use the load_spack script first.}" &&
spack develop --path=$(pwd) exago@develop &&
spack bootstrap now &&
spack concretize -f &&
# This fails for resolve at the moment since it is a private repo
(spack mirror create -a --directory $SPACK_MIRROR || true) &&
spack mirror add local file://$SPACK_MIRROR &&
spack mirror list

res=$?

chmod -R ugo+wrx $SPACK_MIRROR &

exit $res

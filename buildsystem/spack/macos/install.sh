#!/bin/bash -e

export MY_CLUSTER=macos

. buildsystem/spack/load_spack.sh
spack compiler find
spack develop --no-clone --path=$(pwd) exago@git.$(git branch --show-current)=develop
spack install -j 24


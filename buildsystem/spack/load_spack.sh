#!/bin/bash

#assume running from exago root dir

source ./tpl/spack/share/spack/setup-env.sh

SPACKENV=$(pwd)/spack-env

mkdir -p $SPACKENV
cp ./buildsystem/clang-hip/crusher/spack.yaml $SPACKENV

spack env create -d $SPACKENV
spack env activate -p $SPACKENV


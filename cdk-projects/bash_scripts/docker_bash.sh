#!/bin/bash

set -e -o pipefail
#condition for executing a file
if [[ ! -f $PWD/buildsystem/build.sh ]]; then
  echo 'Please run this script from the top-level ExaGO source directory.'
  exit 1
fi
#check spack clone existence 
if [[ ! -f $PWD/tpl/spack/share/spack/setup-env.sh ]]; then
  echo 'Spack not found, please run the command : git submodule update --init --recursive' 
  exit 1
fi
#check for spack.yaml
if [[ ! -f $PWD/cdk-projects/bash_scripts/spack.yaml ]]; then
  echo 'spack.yaml file not found'  
  exit 1
fi
# source spack env 
source tpl/spack/share/spack/setup-env.sh
echo "source to spack env is done"
# Create new directory for spack environment and debugging
SPACK_ENV=spack_cloud_env
mkdir -p $SPACK_ENV
echo "created a dir named $SPACK_ENV"
# Clean out the dir just in case
rm -rf $SPACK_ENV/*
echo "Cleaned out the dir just in case"
# Create spack directory based environment
spack env create -d $SPACK_ENV
echo "Created spack directory based environment"
spack env activate -p $SPACK_ENV
echo "Activated spack directory based environment"

cp ./cdk-projects/bash_scripts/spack.yaml $SPACK_ENV/


cd $SPACK_ENV
#create a docker file
spack containerize > ../Dockerfile
cd -

set +e

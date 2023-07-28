#!/bin/bash

set -e

# source spack env 
source ../../../tpl/spack/share/spack/setup-env.sh
echo "source to spack env is done"
# Create new directory for spack environment and debugging
SPACK_ENV=spack_env
mkdir -p $SPACK_ENV
echo "created a dir named $SPACK_ENV"
# Clean out the dir just in case
rm -rf $SPACK_ENV/*
echo "Cleaned out the dir just in case"
# Create spack directory based environment
spack env create -d $SPACK_ENV
echo "Created spack directory based environment"
spack env activate -p spack_env
echo "Activated spack directory based environment"

# Content to be added to spack.yaml
cat << EOF > ./spack_env/spack.yaml
spack:
  specs:
    - cmake@3.26.3 arch=linux-amzn2-x86_64_v3
    # - exago~mpi~hiop~ipopt+python
  mirrors:
    spack: https://binaries.spack.io/v0.20.1
  container:
    format: docker
    images:
      os: "amazonlinux:2"
      spack: v0.20.1
EOF

cd spack_env
#create a docker file
spack containerize > Dockerfile

# THIS ONLY WORKS ON MAC!!!
sed -i "" "s/spack install/spack buildcache keys --install --trust \&\& spack install/g" Dockerfile

echo "Dockerfile created"
# Build Docker image with Amazon Linux
IMG_TAG=amazonlinux_image
docker build --progress plain -t $IMG_TAG .
echo "Docker image '$IMG_TAG' built."
# Run Docker container
docker run --rm -i -t $IMG_TAG
echo "Docker container 'exago_cont' ran successfully."
cd -

set +e

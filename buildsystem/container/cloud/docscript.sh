#!/bin/bash

set -e

# source spack env 
source ../../../tpl/spack/share/spack/setup-env.sh
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

cp spack.yaml $SPACK_ENV/

cd $SPACK_ENV
#create a docker file
spack containerize > Dockerfile

#copy the coinhsl, make sure to put the coinhsl  tar in the exago root directory
cp ../../../../coinhsl-archive-2019.05.21.tar.gz $SPACK_ENV

COIN_COPY="COPY coinhsl-archive-2019.05.21.tar.gz /opt/spack-environment/coinhsl-archive-2019.05.21.tar.gz"
INSTALL_SOFTWARE_NOTE="# Install the software, remove unnecessary deps"
sed -i "" "s|$INSTALL_SOFTWARE_NOTE|$COIN_COPY\n\n$INSTALL_SOFTWARE_NOTE|" Dockerfile

# THIS ONLY WORKS ON MAC!!!
sed -i "" "s/spack install/spack buildcache keys --install --trust \&\& spack install/g" Dockerfile

sed -i "" "s|CMD \[ \"/bin/bash\" \]|\nRUN yum -y install libgomp \&\& yum -y install libgfortran|" Dockerfile
echo "CMD [ \"/bin/bash\" ]" >> Dockerfile

echo "Dockerfile created"
# Build Docker image with Amazon Linux
IMG_TAG=amazonlinux_image
docker build --progress plain -t $IMG_TAG .
echo "Docker image '$IMG_TAG' built."
# Run Docker container
docker run --rm -i -t $IMG_TAG
echo "Docker container ran successfully."

cd -

set +e

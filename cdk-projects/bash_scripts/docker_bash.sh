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
#copy the coinhsl, make sure to put the coinhsl  tar in the exago root directory
#cp ../../../../coinhsl-archive-2019.05.21.tar.gz $SPACK_ENV

#COIN_COPY="COPY coinhsl-archive-2019.05.21.tar.gz /opt/spack-environment/coinhsl-archive-2019.05.21.tar.gz"
#INSTALL_SOFTWARE_NOTE="# Install the software, remove unnecessary deps"
#sed -i "" "s|$INSTALL_SOFTWARE_NOTE|$COIN_COPY\n\n$INSTALL_SOFTWARE_NOTE|" Dockerfile

# THIS ONLY WORKS ON MAC!!!
sed -i "" "s/spack install/spack buildcache keys --install --trust \&\& spack install/g" Dockerfile

# sed -i "" "s|CMD \[ \"/bin/bash\" \]|\nRUN yum -y install libgomp \&\& yum -y install libgfortran|" Dockerfile
# echo "CMD [ \"/bin/bash\" ]" >> Dockerfile
# jupyter lab --ip 0.0.0.0 --no-browser --allow-root --NotebookApp.token='' --NotebookApp.password='' --ServerApp.root_dir='/home/performance_analysis'

# copy the contents of the datafile
tar -cvf datafiles.tar datafiles 
# copy the contents of the perf file
tar -cvf perf.tar performance_analysis 

echo "Dockerfile created"
# Build Docker image with Amazon Linux
IMG_TAG=amazonlinux_image
docker build --progress plain -t $IMG_TAG .
echo "Docker image '$IMG_TAG' built."
# Run Docker container
docker run --rm -i -t -p 8888:8888 $IMG_TAG
echo "Docker container ran successfully."

rm datafiles.tar
rm perf.tar

cd -

set +e

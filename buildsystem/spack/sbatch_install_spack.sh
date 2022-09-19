#!/bin/bash

#SBATCH -A csc359_crusher
#SBATCH --reservation=hack4
#SBATCH -p batch
#SBATCH -t 240
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -J exasgd_spack_install
#SBATCH -o spack_install.%J
#SBATCH -e spack_install.%J

# Configure https proxy because spack is going to do some things with git
export all_proxy="socks://proxy.ccs.ornl.gov:3128"
export ftp_proxy="ftp://proxy.ccs.ornl.gov:3128"
export http_proxy="http://proxy.ccs.ornl.gov:3128"
export https_proxy="http://proxy.ccs.ornl.gov:3128"
export HTTP_PROXY="http://proxy.ccs.ornl.gov:3128"
export HTTPS_PROXY="http://proxy.ccs.ornl.gov:3128"
export proxy="proxy.ccs.ornl.gov:3128"
export no_proxy='localhost,127.0.0.0/8,*.ccs.ornl.gov,*.olcf.ornl.gov,*.ncrc.gov'

. buildsystem/spack/load_spack.sh
spack config blame config
spack develop --no-clone --path=$(pwd) exago@develop
buildsystem/spack/configure_modules.sh 32


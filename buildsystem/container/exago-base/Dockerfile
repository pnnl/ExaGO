FROM centos:centos8

RUN yum install -y epel-release
RUN dnf install -y 'dnf-command(config-manager)'
RUN yum update -y
RUN dnf config-manager -y --set-enabled powertools plus

RUN yum install -y \
  gcc gcc-c++ \
  openmpi openmpi-devel \
  git \
  cmake m4 make pkgconf diffutils patch libarchive\
  python3 openssl openblas \
  wget curl \
  perl \
  gmp gmp-c++ gmp-devel \
  mpfr mpfr-devel \
  hdf5-openmpi hdf5-openmpi-devel \
  zlib zlib-devel \
  suitesparse suitesparse-devel

# Centos RPM package installs mpi headers into a strange location
RUN cp -R /usr/include/openmpi-x86_64/ /usr/lib64/openmpi/include

# OpenMPI is installed to /usr/lib64/openmpi/ by default, but we always want the
# mpi binaries in our path
RUN ln -s /usr/lib64/openmpi/bin/* /usr/bin/

# Openblas is not found since the RPM only adds libopenblas.so.0
RUN ln -s /usr/lib64/libopenblas.so.0 /usr/lib64/libopenblas.so

# Clone ExaGO, HiOp, and Spack
RUN git clone --depth 1 --branch develop https://github.com/spack/spack.git /spack
 
# Set up bashrc so user has access to everything we've set up
RUN echo 'source /spack/share/spack/setup-env.sh' >> ~/.bashrc
RUN echo 'export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1' >> ~/.bashrc

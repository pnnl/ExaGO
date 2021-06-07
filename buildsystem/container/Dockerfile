FROM centos:centos8

# Uncomment if you're on PNNL's network
# ENV HTTPS_PROXY=http://proxy01.pnl.gov:3128

RUN yum install -y \
  gcc gcc-c++ \
  openmpi openmpi-devel \
  git \
  cmake \
  python3 suitesparse openssl openblas \
  diffutils patch libarchive make pkgconf \
  zlib \
  wget curl

# Centos RPM package installs mpi headers into a strange location
RUN cp -R /usr/include/openmpi-x86_64/ /usr/lib64/openmpi/include

# OpenMPI is installed to /usr/lib64/openmpi/ by default, but we always want the
# mpi binaries in our path
RUN ln -s /usr/lib64/openmpi/bin/* /usr/bin/

# Openblas is not found since the RPM only adds libopenblas.so.0
RUN ln -s /usr/lib64/libopenblas.so.0 /usr/lib64/libopenblas.so

# Initialize spack environment and workspace
RUN mkdir /workspace
COPY centos8.yaml /workspace/environment.yaml

RUN git clone --depth 1 --branch develop https://github.com/spack/spack.git /spack

# Add spack to path by default
RUN echo 'source /spack/share/spack/setup-env.sh' >> ~/.bashrc

# Set up and install spack environment
# RUN /bin/bash -c 'source /spack/share/spack/setup-env.sh; spack env create exago-devel /workspace/environment.yaml'
# RUN /bin/bash -c 'source /spack/share/spack/setup-env.sh; spack env activate exago-devel; spack concretize -f'
# RUN /bin/bash -c 'source /spack/share/spack/setup-env.sh; spack env activate exago-devel; spack install'

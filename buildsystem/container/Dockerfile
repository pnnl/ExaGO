FROM registry.access.redhat.com/ubi8:latest
RUN yum install -y gcc gcc-c++ gcc-gfortran \
  python3 git cmake make tar bzip2 xz patch \
  automake libevent openssl openssl-devel openssl-libs

ENV SPACK_INSTALL_LOCATION=/opt/spack
RUN git clone --branch develop https://github.com/LLNL/spack.git ${SPACK_INSTALL_LOCATION}
RUN . "${SPACK_INSTALL_LOCATION}/share/spack/setup-env.sh"

RUN mkdir /opt/spack-env \
&& (echo "spack:" \
&& echo "  specs:" \
&& echo "    - exago@develop%gcc@8.4.1~ipopt+hiop+mpi~cuda" \
&& echo "     ^hiop@0.4.1 ^openmpi ^petsc@3.13.5+mpi~hypre~hdf5~metis~superlu-dist" \
&& echo "  config:" \
&& echo "    clingo: true" \
&& echo "    install_tree: /opt/sw" \
&& echo "  concretization: together" \
&& echo "  view: false") > /opt/spack-env/spack.yaml

# Find a few system packages to improve build times
RUN for pkg in cmake gmake python git; do \
    spack external find $pkg \
  done

RUN spack compiler find

RUN cd /opt/spack-env \
&& spack env activate . \
&& spack concretize -f \
&& spack install --source --keep-stage -j `nproc`

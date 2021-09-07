FROM spack/ubuntu-bionic:latest

RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - exago@develop~ipopt+hiop+mpi~cuda ^hiop@0.4.1 ^openmpi ^petsc@3.13.5+mpi~hypre~hdf5~metis~superlu-dist" \
&&   echo "  concretization: together" \
&&   echo "  view: /opt/view" \
&&   echo "  config:" \
&&   echo "    clingo: true" \
&&   echo "    build_jobs: 1" \
&&   echo "    install_tree: /opt/software" \
&&   echo "  mirrors:" \
&&   echo "    e4s: https://cache.e4s.io" \
&&   echo "  packages:" \
&&   echo "    all:" \
&&   echo "      target:" \
&&   echo "      - x86_64" \
&&   echo "    petsc:" \
&&   echo "      version:" \
&&   echo "      - 3.13.5" \
&&   echo "      variants: +mpi~hypre~hdf5~metis~superlu-dist") > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && \
  spack env activate . && \
  spack gpg init && \
  spack buildcache keys -it &&
  spack install --fail-fast


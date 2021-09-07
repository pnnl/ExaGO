FROM spack/ubuntu-bionic:latest

RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - exago@develop~ipopt~hiop+mpi~cuda target=linux-ubuntu18.04-x86_64" \
&&   echo "  config:" \
&&   echo "    clingo: true" \
&&   echo "    build_jobs: 1" \
&&   echo "    install_tree: /opt/software" \
&&   echo "  concretization: together" \
&&   echo "  view: /opt/view" \
&&   echo "  mirrors:" \
&&   echo "    e4s: https://cache.e4s.io" \
&&   echo "  packages:" \
&&   echo "    all:" \
&&   echo "      providers:" \
&&   echo "        mpi:" \
&&   echo "        - openmpi" \
&&   echo "      target:" \
&&   echo "      - x86_64" \
&&   echo "    petsc:" \
&&   echo "      version:" \
&&   echo "      - 3.13.5" \
&&   echo "      variants: +mpi~hypre~hdf5~metis~superlu-dist" \
&&   echo "    openmpi:" \
&&   echo "      version:" \
&&   echo "      - 4.1.1" \
&&   echo "      variants: ~atomics~cuda~cxx~cxx_exceptions+gpfs~internal-hwloc~java~legacylaunchers~lustre~memchecker~pmi~singularity~sqlite3+static~thread_multiple+vt+wrapper-rpath fabrics=none schedulers=none") > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && \
  spack env activate . && \
  spack gpg init && \
  spack buildcache keys -it && \
  spack install --fail-fast


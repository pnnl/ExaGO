FROM spack/ubuntu-bionic:latest
ARG AWS_ACCESS_KEY_ID
ARG AWS_SECRET_ACCESS_KEY

RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - pkgconf" \
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
  spack gpg create 'Asher Mancinelli' 'ashermancinelli@gmail.com' && \
  spack buildcache keys -it && \
  export S3_ENDPOINT_URL=http://cache.exasgd.pnl.gov && \
  spack mirror add minio s3://spack && \
  spack mirror add local file:///cache && \
  spack install --fail-fast && \
  mkdir /cache && \
  for ii in $(spack find --format "yyy {version} /{hash}" | \
        grep -v -E "^(develop^master)" | \
        grep "yyy" | \
        cut -f3 -d" "); \
  do \
    spack buildcache create -af -d /cache --only=package $ii ; \
  done && \
  spack buildcache sync --src-directory /cache --dest-mirror-url s3://spack && \
  spack gpg export key.pub && \
  cat key.pub


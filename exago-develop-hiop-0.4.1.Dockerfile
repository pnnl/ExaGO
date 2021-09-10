FROM spack/ubuntu-bionic:latest
ARG AWS_ACCESS_KEY_ID
ARG AWS_SECRET_ACCESS_KEY

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

ENV S3_ENDPOINT_URL=http://cache.exasgd.pnl.gov

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && \
  . /kaniko/s3env.sh && \
  spack env activate . && \
  spack gpg init && \
  spack gpg create 'Asher Mancinelli' 'ashermancinelli@gmail.com' && \
  spack buildcache keys -it && \
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


# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-jammy as builder

# Install OS packages needed to build the software
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install autoconf \
 && rm -rf /var/lib/apt/lists/*

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo spack: \
&&   echo '  specs:' \
&&   echo '  - exago@1.6.0+ipopt+python+mpi' \
&&   echo '  - petsc~fortran~hdf5~hypre+metis' \
&&   echo '  - ipopt@3.12.10+coinhsl~mumps+metis' \
&&   echo '  - coinhsl@2019.05.21+blas' \
&&   echo '  - py-jupyterlab' \
&&   echo '  concretizer:' \
&&   echo '    unify: true' \
&&   echo '    reuse: true' \
&&   echo '  mirrors:' \
&&   echo '    spack: https://binaries.spack.io/develop' \
&&   echo '  packages:' \
&&   echo '    all:' \
&&   echo '      providers:' \
&&   echo '        mpi:' \
&&   echo '        - openmpi' \
&&   echo '        zlib-api:' \
&&   echo '        - zlib' \
&&   echo '  config:' \
&&   echo '    install_tree: /opt/software' \
&&   echo '  view: /opt/views/view') > /opt/spack-environment/spack.yaml

# This is manual, and required to build Ipopt
COPY coinhsl-archive-2019.05.21.tar.gz /opt/spack-environment/coinhsl-archive-2019.05.21.tar.gz

# Do this separate of install to cache keys...
RUN spack mirror add develop https://binaries.spack.io/develop && spack buildcache keys --install --trust

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/views/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . > activate.sh


# Bare OS image to run the installed executables
FROM mcr.microsoft.com/devcontainers/python:3.11-bookworm

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

RUN { \
      echo '#!/bin/sh' \
      && echo '.' /opt/spack-environment/activate.sh \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh \
&& ln -s /opt/views/view /opt/view


RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install gfortran \
 && rm -rf /var/lib/apt/lists/*
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

# Make sure devcontainer user gets spack packages
RUN echo "source /entrypoint.sh" >> /home/vscode/.bashrc

# Install the Jupyter kernel
RUN . /entrypoint.sh && \
       python3 -m ipykernel install \
        --name py311-exago\
        --display-name ExaGO \
        --prefix=

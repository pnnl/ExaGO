FROM amazonlinux:2 as bootstrap

ENV SPACK_ROOT=/opt/spack \
    CURRENTLY_BUILDING_DOCKER_IMAGE=1 \
    container=docker

RUN yum update -y \
 && yum groupinstall -y "Development Tools" \
 && yum install -y \
        curl \
        findutils \
        gcc-c++ \
        gcc \
        gcc-gfortran \
        git \
        gnupg2 \
        hg \
        hostname \
        iproute \
        make \
        patch \
        python3 \
        python3-pip \
        python3-setuptools \
        unzip \
        zstd \
 && pip3 install boto3 \
 && rm -rf /var/cache/yum \
 && yum clean all

RUN mkdir $SPACK_ROOT && cd $SPACK_ROOT && \
    git init --quiet && git remote add origin https://github.com/spack/spack.git && git fetch --depth=1 origin develop && git checkout --detach FETCH_HEAD && \
    mkdir -p $SPACK_ROOT/opt/spack

RUN ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/docker-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/interactive-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/spack-env

RUN mkdir -p /root/.spack \
 && cp $SPACK_ROOT/share/spack/docker/modules.yaml \
        /root/.spack/modules.yaml \
 && rm -rf /root/*.* /run/nologin

# [WORKAROUND]
# https://superuser.com/questions/1241548/
#     xubuntu-16-04-ttyname-failed-inappropriate-ioctl-for-device#1253889
RUN [ -f ~/.profile ]                                               \
 && sed -i 's/mesg n/( tty -s \\&\\& mesg n || true )/g' ~/.profile \
 || true


WORKDIR /root
SHELL ["docker-shell"]

# Creates the package cache
RUN spack bootstrap now \
    && spack bootstrap status --optional \
    && spack spec hdf5+mpi

ENTRYPOINT ["/bin/bash", "/opt/spack/share/spack/docker/entrypoint.bash"]
CMD ["interactive-shell"]

# Build stage with Spack pre-installed and ready to be used
FROM bootstrap as builder


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo spack: \
&&   echo '  specs:' \
&&   echo '    # - cmake@3.26.3 arch=linux-amzn2-x86_64_v3' \
&&   echo '  - exago@develop+python~mpi' \
&&   echo '  concretizer:' \
&&   echo '    unify: when_possible' \
&&   echo '    reuse: true' \
&&   echo '  mirrors:' \
&&   echo '    spack: https://binaries.spack.io/develop' \
&&   echo '  packages:' \
&&   echo '    exago:' \
&&   echo '      require: +raja+hiop+ipopt' \
&&   echo '    hiop:' \
&&   echo '      require: '"'"'@develop+sparse+mpi+ginkgo+kron'"'"'' \
&&   echo '    magma:' \
&&   echo '      require: '"'"'@2.6.2'"'"'' \
&&   echo '    coinhsl:' \
&&   echo '      require: '"'"'@2019.05.21'"'"'' \
&&   echo '    ipopt:' \
&&   echo '      require: '"'"'@3.12.10~metis+coinhsl~mumps'"'"'' \
&&   echo '    raja:' \
&&   echo '      require: ~examples~exercises' \
&&   echo '    umpire:' \
&&   echo '      require: ~openmp~examples' \
&&   echo '    petsc:' \
&&   echo '      require: ~hypre~superlu-dist~hdf5~metis' \
&&   echo '    python:' \
&&   echo '      require: '"'"'@3.9.12'"'"'' \
&&   echo '    all:' \
&&   echo '      compiler:' \
&&   echo '      - gcc@11.2.0' \
&&   echo '      providers:' \
&&   echo '        mpi:' \
&&   echo '        - openmpi' \
&&   echo '  config:' \
&&   echo '    install_tree: /opt/software' \
&&   echo '  view: /opt/views/view') > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack env activate . && spack buildcache keys --install --trust && spack install --fail-fast && spack gc -y

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
FROM amazonlinux:2

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


ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

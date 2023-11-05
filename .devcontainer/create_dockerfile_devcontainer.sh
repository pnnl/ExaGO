#!/bin/bash
# Generates Dockerfile for use as Devcontainer...

# Only necessary for HTTPS operations on PNNL VPN
# export HTTPS_PROXY=http://proxy01.pnl.gov:3128 &&
# export https_proxy=http://proxy01.pnl.gov:3128 &&

# Get CoinHSL tarball if not present
if [ ! -f "./coinhsl-archive-2019.05.21.tar.gz" ]; then
    echo "No coinhsl-archive-2019.05.21.tar.gz found - make sure to get your own license and copy it to this directory"
    echo "Note that the tarball must be named coinhsl-archive-2019.05.21.tar.gz"
    exit 1
fi

# Add read/execute permissions to CoinHSL tarball
chmod +rx ./coinhsl-archive-2019.05.21.tar.gz &&

# If "-f" flag passed, force rebuild
if [ "$1" == "-f" ]; then
    echo "Forcing rebuild" &&
    rm -rf ~/.spack &&
    rm -rf ./spack &&
    rm -rf ./spack-env
fi

#  or no ./spack dir, force a rebuild
if [ ! -d "./tpl/spack/share" ]; then
    echo "No spack submodule directory found, downloading spack"
    git submodule update --init --recursive
fi

export SPACK_PYTHON=$(which python3) &&
# Setup spack
source ./tpl/spack/share/spack/setup-env.sh &&

# Create heredoc for spack.yaml with environment
cat > $(pwd)/spack.yaml <<EOF
spack:
  specs:
    - exago@1.6.0+ipopt+python+mpi
    - petsc~fortran~hdf5~hypre+metis
    - ipopt@3.12.10+coinhsl~mumps+metis
    - coinhsl@2019.05.21+blas
    - py-jupyterlab
    - py-ipyparallel
  concretizer:
    unify: true
    reuse: true
  mirrors:
    spack: https://binaries.spack.io/develop
  container:
    format: docker
    images:
      build: spack/ubuntu-jammy
      final: mcr.microsoft.com/devcontainers/python:3.11-bookworm
    os_packages:
      command: apt
      build:
      - autoconf
      final:
      - gfortran
      - nodejs
      - npm
  packages:
    all:
      providers:
        mpi:
        - openmpi
        zlib-api:
        - zlib
    zlib-ng:
      buildable: false
EOF

# Containerize into a Dockerfile
# Eventually put this straight into ExaGO repo...
spack containerize > ./.devcontainer/Dockerfile

# THIS ONLY WORKS ON MAC!!! - maybe get gsed wit brew working...
# Add Docker command before "# Install" that copies CoinHSL into image
sed -i "" "s|# Install the software|# This is manual, and required to build Ipopt\nCOPY coinhsl-archive-2019.05.21.tar.gz /opt/spack-environment/coinhsl-archive-2019.05.21.tar.gz\n\n# Install the software|" ./.devcontainer/Dockerfile

# Find external packages
sed -i "" "s|# Install the software|# Find external packages\nRUN cd /opt/spack-environment \&\& spack env activate . \&\& spack external find --all\n\n# Install the software|" ./.devcontainer/Dockerfile

# Also trusts the build cache and create source mirror
sed -i "" "s|# Install the software|# Do this separate of install to cache keys...\nRUN cd /opt/spack-environment \&\& spack env activate . \&\& spack mirror add develop https://binaries.spack.io/develop \&\& spack buildcache keys --install --trust \&\& spack concretize -f \&\& spack mirror create -a\n\n# Install the software|" ./.devcontainer/Dockerfile

# Install Ipopt and packages in separate stages to cache builds
sed -i "" "s|# Install the software|# Install Ipopt w/ CoinHSL (and other deps) in stages to cache builds\nRUN cd /opt/spack-environment \&\& spack env activate . \&\& spack install --fail-fast ipopt\n\n# Install the software|" ./.devcontainer/Dockerfile
sed -i "" "s|# Install the software|# Install PETSc \nRUN cd /opt/spack-environment \&\& spack env activate . \&\& spack install --fail-fast petsc\n\n# Install the software|" ./.devcontainer/Dockerfile
sed -i "" "s|# Install the software|# Install py-jupyterlab\nRUN cd /opt/spack-environment \&\& spack env activate . \&\& spack install --fail-fast py-jupyterlab\n\n# Install the software|" ./.devcontainer/Dockerfile
sed -i "" "s|# Install the software|# Install py-ipyparallel\nRUN cd /opt/spack-environment \&\& spack env activate . \&\& spack install --fail-fast py-ipyparallel\n\n# Install the software|" ./.devcontainer/Dockerfile

# Configure environment for VSCode User in devcontainer
echo "# Make sure devcontainer user gets spack packages" >> ./.devcontainer/Dockerfile
echo "RUN echo \"source /entrypoint.sh\" >> /home/vscode/.bashrc" >> ./.devcontainer/Dockerfile

# Finally install the Jupyter kernel
echo "" >> ./.devcontainer/Dockerfile
echo "# Install the Jupyter kernel" >> ./.devcontainer/Dockerfile
echo "RUN . /entrypoint.sh && \\" >> ./.devcontainer/Dockerfile
echo "       python3 -m ipykernel install \\" >> ./.devcontainer/Dockerfile
echo "        --name py311-exago \\" >> ./.devcontainer/Dockerfile
echo "        --display-name \"ExaGO\" \\" >> ./.devcontainer/Dockerfile
echo "        --prefix=$(which jupyter)" >> ./.devcontainer/Dockerfile

# Install the mpi4py Jupyter kernel
echo "" >> ./.devcontainer/Dockerfile
echo "# Install the mpi4py Jupyter kernel" >> ./.devcontainer/Dockerfile
echo "RUN . /entrypoint.sh && \\" >> ./.devcontainer/Dockerfile
echo "       python3 -m ipykernel install \\" >> ./.devcontainer/Dockerfile
echo "        --name py311-mpi4py-exago \\" >> ./.devcontainer/Dockerfile
echo "        --display-name \"ExaGO w/ MPI\" \\" >> ./.devcontainer/Dockerfile
echo "        --prefix=$(which jupyter)" >> ./.devcontainer/Dockerfile

# Modify mpi4py kernel to actually launch process with mpiexec
echo "" >> ./.devcontainer/Dockerfile
echo "# Make modifications that are necessary to run mpi4py in kernelspec" >> ./.devcontainer/Dockerfile
echo "RUN sed -i 's/\\\"\/opt\/views\/view\/bin\/python3\\\"/\\\"mpiexec\\\", \\\"-n\\\", \\\"1\\\", \\\"\/opt\/views\/view\/bin\/python3\\\"/g' /usr/local/share/jupyter/kernels/py311-mpi4py-exago/kernel.json" >> ./.devcontainer/Dockerfile

# Make sure vscode user is correct
echo "" >> ./.devcontainer/Dockerfile
echo "# Configure user for container" >> ./.devcontainer/Dockerfile
echo "USER vscode" >> ./.devcontainer/Dockerfile

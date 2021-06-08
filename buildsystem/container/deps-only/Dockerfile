FROM gitlab.pnnl.gov:4567/exasgd/frameworks/exago/exago-base:latest

# Initialize spack environment and workspace
COPY ./environment.yaml /environment.yaml
COPY ./coinhsl-archive-2015.06.23.tar.gz /root/coinhsl-archive-2015.06.23.tar.gz

# Set up and install spack environment
RUN /bin/bash -c '\
  source /spack/share/spack/setup-env.sh; \
  spack env create exago-devel /environment.yaml'

RUN /bin/bash -c '\
  source /spack/share/spack/setup-env.sh; \
  spack env activate exago-devel; \
  spack concretize -f'

RUN /bin/bash -c '\
  source /spack/share/spack/setup-env.sh; \
  cd /root; \
  spack env activate exago-devel; \
  spack install --verbose -j 4'
 
# Set up bashrc so user has access to everything we've set up
RUN echo 'echo Loading spack environment "exago-devel"...' >> ~/.bashrc
RUN echo 'spack env activate exago-devel' >> ~/.bashrc

FROM gitlab.pnnl.gov:4567/exasgd/frameworks/exago/deps-only:latest

# Clone, build, and install HiOp
RUN git clone --depth 1 --branch v0.4.3 \
  https://github.com/LLNL/hiop.git ~/hiop
ENV HIOP_SRC_DIR=/root/hiop
ENV HIOP_BUILD_DIR=/root/build-hiop
ENV HIOP_INSTALL_DIR=/root/install-hiop
RUN mkdir $HIOP_BUILD_DIR $HIOP_INSTALL_DIR
RUN cd $HIOP_BUILD_DIR && source ~/.bashrc && cmake $HIOP_SRC_DIR \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$HIOP_INSTALL_DIR \
  -DHIOP_BUILD_SHARED=ON -DHIOP_BUILD_STATIC=ON \
  -DHIOP_USE_COINHSL=ON -DHIOP_USE_METIS=ON -DHIOP_USE_MPI=ON -DHIOP_USE_RAJA=ON
RUN cd $HIOP_BUILD_DIR && source ~/.bashrc && make -j 4 install

# Clone, build, and install ExaGO
RUN git clone --depth 1 --branch develop \
  https://gitlab.pnnl.gov/exasgd/frameworks/exago.git ~/exago
ENV EXAGO_SRC_DIR=/root/exago
ENV EXAGO_BUILD_DIR=/root/build-exago
ENV EXAGO_INSTALL_DIR=/root/install-exago
RUN mkdir $EXAGO_BUILD_DIR $EXAGO_INSTALL_DIR
RUN cd $EXAGO_BUILD_DIR && source ~/.bashrc && cmake $EXAGO_SRC_DIR \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$EXAGO_INSTALL_DIR \
  -DEXAGO_BUILD_STATIC=ON -DEXAGO_BUILD_SHARED=ON \
  -DEXAGO_ENABLE_HIOP=ON -DEXAGO_ENABLE_IPOPT=ON -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_PETSC=ON -DEXAGO_ENABLE_RAJA=ON -DHIOP_DIR=$HIOP_INSTALL_DIR
RUN cd $EXAGO_BUILD_DIR && source ~/.bashrc && make -j 4 install

RUN echo 'echo Setting environment variables EXAGO_SRC_DIR, EXAGO_BUILD_DIR, and EXAGO_INSTALL_DIR...' >> ~/.bashrc
RUN echo "export EXAGO_SRC_DIR=$EXAGO_SRC_DIR" >> ~/.bashrc
RUN echo "export EXAGO_BUILD_DIR=$EXAGO_BUILD_DIR" >> ~/.bashrc
RUN echo "export EXAGO_INSTALL_DIR=$EXAGO_INSTALL_DIR" >> ~/.bashrc
RUN echo 'echo Setting environment variables HIOP_SRC_DIR, HIOP_BUILD_DIR, and HIOP_INSTALL_DIR...' >> ~/.bashrc
RUN echo "export HIOP_SRC_DIR=$HIOP_SRC_DIR" >> ~/.bashrc
RUN echo "export HIOP_BUILD_DIR=$HIOP_BUILD_DIR" >> ~/.bashrc
RUN echo "export HIOP_INSTALL_DIR=$HIOP_INSTALL_DIR" >> ~/.bashrc
RUN echo 'cd ~' >> ~/.bashrc

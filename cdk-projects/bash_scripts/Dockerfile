# Base image with ExaGO and Python wrapper installed / configured
FROM ghcr.io/pnnl/exago:exago-develop-mcpznzjpfrag4giqpzisvz5rsd6yl7rf.spack
RUN apt-get -yqq update \
    && apt-get -yqq install curl \
    && curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.11 get-pip.py \
    && rm get-pip.py \
    && apt-get -yqq purge curl \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /opt/spack-environment \
&&  (echo spack: \
&&   echo '  specs:' \
&&   echo '  - exago@1.6.0+ipopt+python+mpi' \
&&   echo '    # This is currently a bug in Spack(?)' \
&&   echo '  - libffi@3.4.2' \
&&   echo '  - petsc~fortran~hdf5~hypre+metis' \
&&   echo '  - ipopt@3.12.10+coinhsl~mumps+metis' \
&&   echo '  - coinhsl@2019.05.21+blas' \
&&   echo '  - py-jupyterlab' \
&&   echo '  - py-ipyparallel' \
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
&&   echo '    zlib-ng:' \
&&   echo '      buildable: false' \
&&   echo '  config:' \
&&   echo '    install_tree: /opt/software' \
&&   echo '  view: /opt/views/view') > /opt/spack-environment/spack.yaml

RUN python3.11 -m pip install --upgrade pip && \
    python3.11 -m pip install awslambdaric boto3

# Define lambda task root
RUN mkdir -p /var/task 
ENV LAMBDA_TASK_ROOT /var/task

# THIS DIRECTORY NEEDS TO BE CREATED, ALBEIT IT IS EMPTY
RUN mkdir -p /opt/extensions 

# Set PYTHONPATH to include /var/task
ENV PYTHONPATH "${PYTHONPATH}:${LAMBDA_TASK_ROOT}"
# Copy function code
RUN  mkdir -p ${LAMBDA_TASK_ROOT}
COPY lambda_handler.py ${LAMBDA_TASK_ROOT}
ADD python_wrapper ${LAMBDA_TASK_ROOT}/python_wrapper

# Replicate Spack package path structure as per PR suggestion
# Ensure full path structure for Spack package is created#RUN mkdir -p /path/to/prefix \


# Set up the entry point and command
WORKDIR /var/task
ENTRYPOINT ["python3.11", "-m", "awslambdaric"]
CMD [ "lambda_handler.lambda_handler" ]

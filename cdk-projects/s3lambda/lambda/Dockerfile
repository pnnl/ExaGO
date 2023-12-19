# Base image with ExaGO and Python wrapper installed / configured
FROM ghcr.io/pnnl/exago:exago-develop-f6cnvd7ffkspgdys3myedu4ha5wnanba.spack

RUN apt-get -yqq update \
    && apt-get -yqq install curl \
    && curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.11 get-pip.py \
    && rm get-pip.py \
    && apt-get -yqq purge curl \
    && rm -rf /var/lib/apt/lists/*


RUN python3.11 -m pip install --upgrade pip && \
    python3.11 -m pip install awslambdaric boto3

EXPOSE 9000

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
# Ensure full path structure for Spack package is created
RUN mkdir -p /path/to/prefix \
    && mkdir -p /path/to/prefix/bin \
    && mkdir -p /path/to/prefix/.spack/binary_distribution \
    && mkdir -p /opt/extensions

# Replace the following path with the correct one as per your Spack package
RUN mkdir -p /opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_p/linux-ubuntu22.04-x86_64_v2/gcc-13.1.0/openblas-0.3.25-szvnjls57k2pezgoewgal43airoxuvt7/lib
RUN mkdir -p /var/task/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_p/linux-ubuntu22.04-x86_64_v2/gcc-13.1.0/openblas-0.3.25-szvnjls57k2pezgoewgal43airoxuvt7/lib

# Hotfixes
RUN mkdir -p /opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_p/linux-ubuntu22.04-x86_64_v2/gcc-13.1.0/openblas-0.3.25-szvnjls57k2pezgoewgal43airoxuvt7/lib
RUN mkdir -p /var/task/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_p/linux-ubuntu22.04-x86_64_v2/gcc-13.1.0/openblas-0.3.25-szvnjls57k2pezgoewgal43airoxuvt7/lib
RUN echo "Running LDD"
RUN ldd /opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/linux-ubuntu22.04-x86_64_v2/gcc-13.1.0/exago-develop-f6cnvd7ffkspgdys3myedu4ha5wnanba/lib/python3.11/site-packages/exago.cpython-311-x86_64-linux-gnu.so

# Set up the entry point and command
WORKDIR /var/task
ENTRYPOINT ["python3.11", "-m", "awslambdaric"]
CMD [ "lambda_handler.lambda_handler" ]

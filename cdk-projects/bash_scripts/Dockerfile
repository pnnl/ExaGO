# Base image with ExaGO and Python wrapper installed / configured
FROM ghcr.io/pnnl/exago:exago-develop-mcpznzjpfrag4giqpzisvz5rsd6yl7rf.spack
RUN apt-get -yqq update \
    && apt-get -yqq install curl \
    && curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.11 get-pip.py \
    && rm get-pip.py \
    && apt-get -yqq purge curl \
    && rm -rf /var/lib/apt/lists/*


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

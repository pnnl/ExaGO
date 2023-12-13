FROM ghcr.io/pnnl/exago:exago-develop-x436raytwpfhkihapsbpto3jd3dpfyqg.spack
# Copy requirements.txt
#COPY requirements.txt ${LAMBDA_TASK_ROOT}
 #Define lambda task root
ENV LAMBDA_TASK_ROOT /var/task

# Set PYTHONPATH to include /var/task
ENV PYTHONPATH "${PYTHONPATH}:${LAMBDA_TASK_ROOT}"
# Copy function code
RUN  mkdir -p ${LAMBDA_TASK_ROOT}
COPY lambda_handler.py ${LAMBDA_TASK_ROOT}
ADD python_wrapper ${LAMBDA_TASK_ROOT}/python_wrapper
 
#Setup Proxy
#RUN export HTTPS_PROXY=http://proxy01.pnl.gov:3128
#RUN export https_proxy=http://proxy01.pnl.gov:3128
# Install the specified packages
#RUN export HTTP_PROXY=http://proxy01.pnl.gov:3128
#RUN export http_proxy=http://proxy01.pnl.gov:3128
#RUN apt-get -yqq update && apt-get -yqq upgrade \
#&& apt-get -yqq install python3-pip \
#&& rm -rf /var/lib/apt/lists/*
RUN apt-get -yqq update \
    && apt-get -yqq install curl \
    && curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.11 get-pip.py \
    && rm get-pip.py \
    && apt-get -yqq purge curl \
    && rm -rf /var/lib/apt/lists/*

# Assuming Python 3.11.6 has been installed
#RUN curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    #&& python3.11 get-pip.py \
    #&& rm get-pip.py

RUN python3.11 -m pip install --upgrade pip
#RUN pip install -r ${LAMBDA_TASK_ROOT}/requirements.txt
# Make port 9000 available to the world outside this container
RUN python3.11 -m pip install awslambdaric
EXPOSE 9000
RUN ls /
RUN ls /opt/spack/
RUN python3.11 -v
#ENV LD_LIBRARY_PATH, PYTHONPATH
#RUN python -v
RUN ls /opt/

# Setting environment variables
#ENV PYTHONPATH "${PYTHONPATH}:/path/to/python3.11.6/site-packages"
ENV PYTHONPATH "${PYTHONPATH}:/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_p/linux-ubuntu22.04-x86_64_v2/gcc-13.1.0/python-3.11.6-qgqyp4wogl4b7h4tzepizevq4bopvjry/lib/python3.11/site-packages"
#ENV LD_LIBRARY_PATH "/path/to/ld/library"
#ADD aws-lambda-rie /usr/local/bin/aws-lambda-rie
#ENV PYTHONPATH "${PYTHONPATH}:/usr/local/lib/python3.11.6/dist-packages"
ENTRYPOINT ["python3.11", "-m", "awslambdaric"]

 
# Set the CMD to your handler- filename. lambda_handler-function name (could also be done as a parameter override outside of the Dockerfile)
CMD [ "lambda_handler.lambda_handler" ]

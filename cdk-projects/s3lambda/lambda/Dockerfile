FROM ghcr.io/pnnl/exago:exago-develop-x436raytwpfhkihapsbpto3jd3dpfyqg.spack
# Copy requirements.txt
#COPY requirements.txt ${LAMBDA_TASK_ROOT}
 
# Copy function code
COPY lambda_handler.py ${LAMBDA_TASK_ROOT}
ADD python_wrapper ${LAMBDA_TASK_ROOT}/python_wrapper
 
#Setup Proxy
#RUN export HTTPS_PROXY=http://proxy01.pnl.gov:3128
#RUN export https_proxy=http://proxy01.pnl.gov:3128
# Install the specified packages
#RUN export HTTP_PROXY=http://proxy01.pnl.gov:3128
#RUN export http_proxy=http://proxy01.pnl.gov:3128
RUN apt-get -yqq update && apt-get -yqq upgrade \
&& apt-get -yqq install python3-pip \
&& rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip
#RUN pip install -r ${LAMBDA_TASK_ROOT}/requirements.txt
RUN pip install awslambdaric
ENTRYPOINT ["python", "-m", "awslambdaric"]
ENV PYTHONPATH "${PYTHONPATH}:/usr/local/lib/python3.10/dist-packages"
 
# Set the CMD to your handler- filename. lambda_handler-function name (could also be done as a parameter override outside of the Dockerfile)
CMD [ "lambda_handler.lambda_handler" ]

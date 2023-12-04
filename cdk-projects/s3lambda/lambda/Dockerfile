FROM ghcr.io/pnnl/exago:exago-develop-x436raytwpfhkihapsbpto3jd3dpfyqg.spack
# Copy requirements.txt
#COPY requirements.txt ${LAMBDA_TASK_ROOT}
 #Define lambda task root
ENV LAMBDA_TASK_ROOT /var/task
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
RUN apt-get -yqq update && apt-get -yqq upgrade \
&& apt-get -yqq install python3-pip \
&& rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip
#RUN pip install -r ${LAMBDA_TASK_ROOT}/requirements.txt
RUN pip install awslambdaric
# Make port 9000 available to the world outside this container
EXPOSE 9000
#ADD aws-lambda-rie /usr/local/bin/aws-lambda-rie
ENTRYPOINT ["python", "-m", "awslambdaric"]
ENV PYTHONPATH "${PYTHONPATH}:/usr/local/lib/python3.10/dist-packages"
 
# Set the CMD to your handler- filename. lambda_handler-function name (could also be done as a parameter override outside of the Dockerfile)
CMD [ "lambda_handler.lambda_handler" ]

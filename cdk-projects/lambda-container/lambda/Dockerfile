FROM public.ecr.aws/lambda/python:3.11
 
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
RUN pip install --upgrade pip
#RUN pip install -r ${LAMBDA_TASK_ROOT}/requirements.txt

 
# Set the CMD to your handler (could also be done as a parameter override outside of the Dockerfile)
CMD [ "lambda_handler.main" ]

#ENTRYPOINT [ "executable" ]
#WORKDIR /the/workdir/path
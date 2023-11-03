#!/bin/bash

cp ../../Dockerfile ./lambda

# TODO:

# 0. Get ExaGO building with `docker_bash.sh`, along with a `git submodule update --init --recursive` after a `git rebase -i develop` after a `git checkout develop && git pull`

# 0.5. Run base built image and verify `python -c "import exago"` works (exago is on python path)

# 1. Add necessary AWS configuration to trigger lambda handler

# Here is what is currently added
# ENV LAMBDA_TASK_ROOT /var/root
# RUN echo ${LAMBDA_TASK_ROOT}
# 
# # Copy function code
# RUN mkdir -p ${LAMBDA_TASK_ROOT}
# COPY lambda_handler.py ${LAMBDA_TASK_ROOT}
# ADD python_wrapper ${LAMBDA_TASK_ROOT}/python_wrapper
#  
# #Setup Proxy
# RUN pip install --upgrade pip && \
#     pip install \
#     --target ${LAMBDA_TASK_ROOT} \
#         awslambdaric
# 
# # Set runtime interface client as default command for the container runtime
# ENTRYPOINT [ "/var/lang/bin/python", "-m", "awslambdaric" ]
#  
# # Set the CMD to your handler (could also be done as a parameter override outside of the Dockerfile)
# WORKDIR ${LAMBDA_TASK_ROOT}
# CMD [ "lambda_handler.main" ]

# 2. I expect `import exago` will fail. That is ok, work with Paul and Cameron on fixing entrypoint/cmd

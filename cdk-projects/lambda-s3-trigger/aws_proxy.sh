#!/bin/bash

# Setting up proxy and git
export HTTPS_PROXY=http://proxy01.pnl.gov:3128 
export https_proxy=http://proxy01.pnl.gov:3128
export AWS_PROFILE=AdministratorAccess-305402452870
export AWS_DEFAULT_REGION=us-west-2

# Configure Git credential helper
git config --global credential.helper '!aws --profile AdministratorAccess-305402452870 codecommit credential-helper $@'

# Start docker engine 
docker --version

#virtual environment setup
python -m venv venv
source venv/bin/activate


#cdk library setup
pip install aws-cdk-lib==2.90.0

# AWS configuration on CLI
aws configure sso

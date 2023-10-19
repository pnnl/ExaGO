#!/bin/bash

# Setting up proxy and git
#export HTTPS_PROXY=http://proxy01.pnl.gov:3128 export https_proxy=http://proxy01.pnl.gov:3128
#export AWS_PROFILE=AdministratorAccess-305402452870

# Configure Git credential helper
git config --global credential.helper '!aws --profile AdministratorAccess-305402452870 codecommit credential-helper $@'

# Start docker engine 
docker --version

# AWS configuration on CLI
aws configure sso

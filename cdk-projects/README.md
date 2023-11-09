# ExaGO Cloud Deployment and Containerization

Deploying ExaGO to cloud runtimes using Spack brings an ease of maintenance and efficient dependency management because of reproducibility and resource management.
The cloud workflow contains the following. 

1. Proxy configuration
2. AWS authentication
3. Bootstrap (Framework building)
4. Deployment
5. Containerization
6. Code Commit CI/CD

## Proxy Configuration 
Steps to enable AWS CLI tools and SSO while outside the PNNL network and within the VPN. The workflow assumes a Linux or OSX environment with a functional Python installation. Please refer to the reference links for instructions for PowerShell.
The instructions describes how to configure the AWS CLI to authenticate users with AWS IAM Identity Center (successor to AWS Single Sign-On) (IAM Identity Center) using the SSO token provider configuration, your AWS SDK or tool can automatically retrieve refreshed authentication tokens.

### Section 1: Configure the AWS CLI tool and Proxy
1.	Install the latest version of the AWS CLI v2 tool: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html
2.	Install the Python library "git-remote-codecommit": pip install git-remote-codecommit (or via your favorite python package manager)
3.	Set up the proxy (in this example, they’re commands to be run in the shell environment; check your system documentation)

```
export HTTPS_PROXY=http://proxy01.pnl.gov:3128 
export https_proxy=http://proxy01.pnl.gov:3128
```

## Directory Structure

Here is the output of `tree -L 2`:

```
.
├── README.md
├── aws_proxy.sh
├── bash_scripts
├── code_build
├── container_codebuild
├── lambda-s3-trigger
├── new_exago
├── spack_cloud_env
├── test_codebuild
└── testcont_codebuild
```

Contents of each directory:

TODO: Have each folder have it's own small README

- README.md
    - This is the document you are reading
- aws_proxy.sh
    - This is for use on your local machine to configure proxy while on VPN
- bash_scripts
    - docker_bash.sh
        - Used to re-generate Dockerfile with `spack.yaml`
    - Dockerfile
        - Output from `docker_bash.sh`
    - spack.yaml
        - Spec of ExaGO used to generate dockerfile
- code_build (TODO: need to delete)
- container_codebuild
    - Contains CodeBuild pipelines to build image and push to ECR
    - TODO:
        - Remove duplicate requirements.txt
        - Document need to keep:
            - codeartifact
            - tests
        - Move python venv config into main `aws_proxy.sh`
- lambda-s3-trigger (TODO: need to delete)
    - Scripts to deploy lambda and associated s3
- new_exago (ExaGO backend lambda + s3 trigger)
    - Represents latest ExaGO builds/development
- spack_cloud_env (TODO: delete in favor of bash_scripts/spack.yaml)
- test_codebuild (TODO: need to delete)
- testcont_codebuild (TODO: need to delete)

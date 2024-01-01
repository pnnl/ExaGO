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

## Proxy and Prerequesties Setup File
The aws_proxy.sh script performs the following tasks

+ Proxy configuration.
+ Install the AWS Cloud Development Kit via the command line.
+ Start the Docker engine.
+ Create a python virtual environment and activate it.
+ Initialize the AWS session configuration.

Execute the aws_proxy.sh script with the commands below

```
 source aws_proxy.sh
```


Note: Start Docker daemon application separately for any containerization tasks.

## Steps to Deploy and Import ExaGO on Cloud

Once the AWS session configuration is done then create a aws bootstrap environment.The cdk bootstrap command is part of the AWS Cloud Development Kit (CDK) and is used to set up the necessary resources in your AWS account for deploying CDK applications.  

```
cdk bootstrap
```

The command ``` cdk synth ``` is used to synthesize or generate AWS CloudFormation templates from the CDK application code. 

```
cdk synth
```

Deploy the CDK application with  ``` cdk deploy ``` command. 

``` 
cdk deploy 
````

## Directory Structure

Here is the output of `tree -L 2`:

```
.
├── README.md
├── aws_proxy.sh
├── bash_scripts
├── container_codebuild
├── s3lambda
```

Contents of each directory:


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
    - build-test.sh
        - This file is to build a docker container with ExaGO base image.
        Note : AWS credentials needs to renew for every run.
        
- container_codebuild
    - Contains CodeBuild pipelines to build image and push to ECR

- s3lambda 
    - This stack is to build a local container with ExaGO and trigger s3lambda. 

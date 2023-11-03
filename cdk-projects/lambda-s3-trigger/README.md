# lambda-s3-trigger

## Goal / Purpose

- Build custom Dockerfile for use in lambda (when using `DockerImageCode.from_image_asset`)
- Eventually, use ECR pre-built image, and have lambda `import exago` and construct OPFLOW object

## Contents

```
.
├── README.md
├── app.py
├── cdk.json
├── lambda
│   ├── Dockerfile
│   ├── lambda_handler.py
│   └── python_wrapper
├── requirements.txt
└── s3trigger
    ├── __init__.py
    ├── __pycache__
    └── s3trigger_stack.py
```

- `app.py`
    - Determines the type of application configured during `cdk-synth`
    - In this case, configured an S3 trigger to run using a Docker image
- `cdk.json`
    - Core configuration for `cdk-synth`, points to `python3 app.py`
- `requirements.txt`
    - Local venv file that ensures the correct version of `aws-cdk-lib` is used
- `s3trigger`
    - `__init__.py`
        - Even though file is empty, it is required
        - Auto-generated, so no need to track this in git?
    - `s3trigger_stack.py`
        - This defines the:
            - Notification that the S3 bucket generates
            - Defines lambda function to be ran
            - If using `DockerImageCode.from_image_asset`, this defines folder where a Dockerfile and configuration live:
                - This configuration is built with a unique name during `cdk-deploy`, and our lambda automatically uses that unique name
            - TODO: If this uses `from_ecr`, then just pulls from existing ecr image
- `lambda`
    - `Dockerfile`
        - Contents to be build if using `DockerImageCode.from_image_asset`
    - `lambda_handler.py`
        - CURRENTLY UNUSED
        - When working, our lambda will call this Python function to handle the trigger
    - `python_wrapper`
        - CURRENTLY UNUSED
        - When working, might include useful utility libraries

## Workflows

1. cdk-synth 
    1. Generates concrete cloud configuration from abstract specification in `.`
1. cdk-deploy
    1. Pushes changes to AWS with updated configuration from `cdk-synth`

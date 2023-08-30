# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0

"""
Re-usable construct to create codebuild job for publishing packages
"""

from constructs import Construct
from aws_cdk import aws_codebuild as codebuild
from aws_cdk import aws_iam as iam


class BuildAndPublishPackage(Construct):

    def __init__(self, scope: Construct, construct_id: str, project_name=None, codeartifact_domain_arn=None, codeartifact_repo_arn=None, artifact_bucket_encryption_key=None,
                 **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)

        self.project = codebuild.PipelineProject(
            self,
            project_name,
            environment=codebuild.BuildEnvironment(
                privileged=False,
                compute_type=codebuild.ComputeType.MEDIUM,
                build_image=codebuild.LinuxBuildImage.STANDARD_5_0
            ),
            encryption_key=artifact_bucket_encryption_key,
            build_spec=codebuild.BuildSpec.from_object({
                "version": "0.2",
                "phases": {
                    "pre_build": {
                        "commands": [
                            "python3 -m venv .venv",
                            ". .venv/bin/activate",
                            "aws codeartifact login --tool pip --repository pip --domain aws-sample-domain",
                            "pip3 install -r requirements.txt",
                            "pip3 install -r requirements-dev.txt",
                            "pip3 install setuptools wheel build twine",
                            "cd ./packages/" + project_name,
                        ],
                    },
                    "build": {
                        "commands": [
                            "python3 -m build .",
                            "export TWINE_USERNAME=aws",
                            "export TWINE_PASSWORD=`aws codeartifact get-authorization-token --domain aws-sample-domain --query authorizationToken --output text`",
                            "export TWINE_REPOSITORY_URL=`aws codeartifact get-repository-endpoint --domain aws-sample-domain --repository pip --format pypi --query repositoryEndpoint --output text`",
                            "twine upload --repository codeartifact dist/*"
                        ],
                    },
                },
            })
        )

        self.project.role.attach_inline_policy(
            iam.Policy(
                self,
                "PublishPolicy",
                statements=[
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=["*"],
                        actions=["sts:GetServiceBearerToken"],
                        conditions={
                            "StringEquals": {
                                "sts:AWSServiceName": "codeartifact.amazonaws.com"
                            },
                        }
                    ),
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=[codeartifact_domain_arn],
                        actions=["codeartifact:GetAuthorizationToken"],
                    ),
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=[codeartifact_repo_arn],
                        actions=[
                            "codeartifact:ReadFromRepository",
                            "codeartifact:GetRepositoryEndpoint",
                            "codeartifact:List*",
                        ],
                    ),
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=["*"],
                        actions=[
                            "codeartifact:PublishPackageVersion"
                        ],
                    )
                ]
            )
        )
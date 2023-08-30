# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0

"""
Sample stack for publishing python packages
"""

from aws_cdk import (
    RemovalPolicy,
    Stack
)
from constructs import Construct
from aws_cdk import aws_iam as iam
from aws_cdk import aws_codebuild as codebuild
from aws_cdk import aws_codecommit as codecommit
from aws_cdk import aws_codeartifact as codeartifact
from aws_cdk import aws_codepipeline as codepipeline
from aws_cdk import aws_codepipeline_actions as codepipeline_actions
from aws_cdk import aws_s3 as s3
from aws_cdk import aws_kms as kms
#from codeartifact.custom_constructs.build_and_publish_package import BuildAndPublishPackage
#from cdk_nag import NagSuppressions
#from cdk_nag import NagPackSuppression


class CodeartifactStack(Stack):

    def __init__(self, scope: Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)
    #edit the repo class with aws documentation 
        repo = codecommit.Repository.from_repository_name(self,"exagocodecommitrepo",
                                                          repository_name="exago-codebuild",
        )

        codeartifact_domain = codeartifact.CfnDomain(
            self,
            "CodeArtifactDomain",
            domain_name="aws-sample-domain",
        )

        pip_private_codeartifact_repository = codeartifact.CfnRepository(
            self,
            "PipPrivateCodeArtifactRepository",
            domain_name=codeartifact_domain.domain_name,
            repository_name="pip",
            description="Private PyPi repo",
            external_connections=["public:pypi"],
        )

        pip_private_codeartifact_repository.add_depends_on(codeartifact_domain)

        codebuild_encryption_key = kms.Key(
            self,
            'codeBuildEncryptionKey',
            enable_key_rotation=True
        )

        access_logs_bucket = s3.Bucket(
            self,
            "AccessLogsBucket",
            bucket_name="sample-cdk-access-logs-" + self.account,
            block_public_access=s3.BlockPublicAccess.BLOCK_ALL,
            encryption=s3.BucketEncryption.KMS,
            encryption_key=codebuild_encryption_key,
            enforce_ssl=True,
            removal_policy=RemovalPolicy.DESTROY,
            auto_delete_objects=True,
        )

        pipeline_artifact_bucket = s3.Bucket(
            self,
            "PipelineArtifactBucket",
            bucket_name="sample-cdk-artifact-" + self.account,
            server_access_logs_bucket=access_logs_bucket,
            block_public_access=s3.BlockPublicAccess.BLOCK_ALL,
            encryption=s3.BucketEncryption.KMS,
            encryption_key=codebuild_encryption_key,
            enforce_ssl=True,
            removal_policy=RemovalPolicy.DESTROY,
            auto_delete_objects=True,
        )

        #NagSuppressions.add_resource_suppressions(
        #    access_logs_bucket,
        #    [NagPackSuppression(id="AwsSolutions-S1", reason="Cannot log to itself")],
        #   True
        #)

        pipeline = codepipeline.Pipeline(
            self,
            "PackagePipeline",
            pipeline_name="python-sample-pipeline",
            restart_execution_on_update=True,
            artifact_bucket=pipeline_artifact_bucket,
        )

        source_output = codepipeline.Artifact()

        source_action = codepipeline_actions.CodeCommitSourceAction(
            action_name="CodeCommit",
            repository=repo,
            output=source_output,
            branch="cloud-dev"
        )

        pipeline.add_stage(
            stage_name="Source",
            actions=[source_action],
        )

        run_build_exago_project = codebuild.PipelineProject(
            self,
            "RunBuildExaGO",
            environment=codebuild.BuildEnvironment(
                privileged=False,
                compute_type=codebuild.ComputeType.MEDIUM,
                build_image=codebuild.LinuxBuildImage.STANDARD_5_0
            ),
    
            encryption_key=pipeline.artifact_bucket.encryption_key,
            build_spec=codebuild.BuildSpec.from_object({
                "version": "0.2",
                "phases": {
                    "pre_build": {
                        "commands": [
                            "python3 -m venv .venv",

                            #"pip3 install -r requirements-dev.txt",
                        ],
                    },
                    "build": {
                        "commands": [
                            ". .venv/bin/activate",
                            "aws codeartifact login --tool pip --repository pip --domain aws-sample-domain",
                            "ls *",
                            #"ls cdk-projects/bash_scripts",
                            #"cd cdk-projects/bash_scripts",
                           "pip3 install -r requirements.txt",
                        ],
                    },
                },
            })
        )

        run_build_exago_project.role.attach_inline_policy(
            iam.Policy(
                self,
                "RunBuildExaGOPolicy",
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
                        resources=[codeartifact_domain.attr_arn],
                        actions=["codeartifact:GetAuthorizationToken"],
                    ),
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=[pip_private_codeartifact_repository.attr_arn],
                        actions=[
                            "codeartifact:ReadFromRepository",
                            "codeartifact:GetRepositoryEndpoint",
                            "codeartifact:List*"
                            #"ecr:GetAuthorizationToken"
                        ],
                    ),
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=["*"],
                        actions=[
                            "ecr:GetAuthorizationToken"
                        ],
                    )
                ]
            )
        )

        pipeline.add_stage(
            stage_name="Build",
            actions=[
                codepipeline_actions.CodeBuildAction(
                    action_name="build-exago-container",
                    project=run_build_exago_project,
                    input=source_output,
                )
            ],
        )

  

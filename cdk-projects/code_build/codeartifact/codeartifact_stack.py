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
from codeartifact1.custom_constructs.build_and_publish_package import BuildAndPublishPackage
from cdk_nag import NagSuppressions
from cdk_nag import NagPackSuppression


class CodeartifactStack(Stack):

    def __init__(self, scope: Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)
    # edit the repo class with aws documentation
        repo = codecommit.Repository.from_repository_name(self, "exagocodecommitrepo",
                                                          repository_name="test_repo",
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

        NagSuppressions.add_resource_suppressions(
            access_logs_bucket,
            [NagPackSuppression(id="AwsSolutions-S1",
                                reason="Cannot log to itself")],
            True
        )

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

        run_unit_tests_project = codebuild.PipelineProject(
            self,
            "RunUnitTests",
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
                            ". .venv/bin/activate",
                            "aws codeartifact login --tool pip --repository pip --domain aws-sample-domain",
                            "ls *",
                            "ls cdk-projects/lambda-s3-trigger",
                            "cd cdk-projects/lambda-s3-trigger",
                           "pip3 install -r requirements.txt",
                            # "pip3 install -r requirements-dev.txt",
                        ],
                    },
                    "build": {
                        "commands": ["pytest"],
                    },
                },
            })
        )

        run_unit_tests_project.role.attach_inline_policy(
            iam.Policy(
                self,
                "RunUnitTestsPolicy",
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
                        resources=[
                            pip_private_codeartifact_repository.attr_arn],
                        actions=[
                            "codeartifact:ReadFromRepository",
                            "codeartifact:GetRepositoryEndpoint",
                            "codeartifact:List*"
                        ],
                    )
                ]
            )
        )

        pipeline.add_stage(
            stage_name="Test",
            actions=[
                codepipeline_actions.CodeBuildAction(
                    action_name="run-unit-tests",
                    project=run_unit_tests_project,
                    input=source_output,
                )
            ],
        )

        self_mutate_project = codebuild.PipelineProject(
            self,
            "SelfMutate",
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
                            "npm install -g aws-cdk",
                            "python3 -m venv .venv",
                            ". .venv/bin/activate",
                            "aws codeartifact login --tool pip --repository pip --domain aws-sample-domain",
                            "ls *",
                            "echo $PWD",
                            "cd test_repo/cdk-projects/new_exago"
                            "pip3 install -r requirements.txt",
                            # "pip3 install -r requirements-dev.txt",
                        ],
                    },
                    "build": {
                        "commands": ["cdk deploy --require-approval=never"],
                    },
                },
            })
        )

        self_mutate_project.role.attach_inline_policy(
            iam.Policy(
                self,
                "SelfMutatePolicy",
                statements=[
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=["*"],
                        actions=["cloudformation:DescribeStacks"],
                    ),
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=["*"],
                        actions=["iam:PassRole"],
                    ),
                    iam.PolicyStatement(
                        effect=iam.Effect.ALLOW,
                        resources=["arn:aws:iam::*:role/cdk-*"],
                        actions=["sts:AssumeRole"],
                    ),
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
                        resources=[
                            pip_private_codeartifact_repository.attr_arn],
                        actions=[
                            "codeartifact:ReadFromRepository",
                            "codeartifact:GetRepositoryEndpoint",
                            "codeartifact:List*",
                        ],
                    ),
                ]
            )
        )

        pipeline.add_stage(
            stage_name="UpdatePipeline",
            actions=[
                codepipeline_actions.CodeBuildAction(
                    action_name="self-mutate",
                    project=self_mutate_project,
                    input=source_output,
                )
            ],
        )

        sample_package_project = BuildAndPublishPackage(
            self,
            "BuildSamplePackage",
            project_name="sample-package",
            artifact_bucket_encryption_key=pipeline.artifact_bucket.encryption_key,
            codeartifact_domain_arn=codeartifact_domain.attr_arn,
            codeartifact_repo_arn=pip_private_codeartifact_repository.attr_arn
        )

        pipeline.add_stage(
            stage_name="BuildAndPublishPackages",
            actions=[
                codepipeline_actions.CodeBuildAction(
                    action_name="sample-package",
                    project=sample_package_project.project,
                    input=source_output,
                )
            ],
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/PackagePipeline/Role/DefaultPolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Action::s3:Abort*",
                    "Action::s3:DeleteObject*",
                    "Action::s3:GetBucket*",
                    "Action::s3:GetObject*",
                    "Action::s3:List*",
                    "Action::kms:GenerateDataKey*",
                    "Action::kms:ReEncrypt*",
                    "Resource::<PipelineArtifactBucketD127CCF6.Arn>/*"
                ])
            ]
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/PackagePipeline/Source/CodeCommit/CodePipelineActionRole/DefaultPolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Action::s3:Abort*",
                    "Action::s3:DeleteObject*",
                    "Action::s3:GetBucket*",
                    "Action::s3:GetObject*",
                    "Action::s3:List*",
                    "Action::kms:GenerateDataKey*",
                    "Action::kms:ReEncrypt*",
                    "Resource::<PipelineArtifactBucketD127CCF6.Arn>/*"
                ])
            ]
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/RunUnitTests/Role/DefaultPolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Resource::arn:<AWS::Partition>:logs:<AWS::Region>:<AWS::AccountId>:log-group:/aws/codebuild/<RunUnitTests2AD5FFEA>:*",
                    "Resource::arn:<AWS::Partition>:codebuild:<AWS::Region>:<AWS::AccountId>:report-group/<RunUnitTests2AD5FFEA>-*",
                    "Action::s3:GetBucket*",
                    "Action::s3:GetObject*",
                    "Action::s3:List*",
                    "Action::kms:GenerateDataKey*",
                    "Action::kms:ReEncrypt*",
                    "Resource::<PipelineArtifactBucketD127CCF6.Arn>/*"
                ])
            ]
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/RunUnitTestsPolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Resource::*",
                    "Action::codeartifact:List*",
                ])
            ]
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/SelfMutate/Role/DefaultPolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Resource::arn:<AWS::Partition>:logs:<AWS::Region>:<AWS::AccountId>:log-group:/aws/codebuild/<SelfMutate95ADA46F>:*",
                    "Resource::arn:<AWS::Partition>:codebuild:<AWS::Region>:<AWS::AccountId>:report-group/<SelfMutate95ADA46F>-*",
                    "Action::s3:GetBucket*",
                    "Action::s3:GetObject*",
                    "Action::s3:List*",
                    "Action::kms:GenerateDataKey*",
                    "Action::kms:ReEncrypt*",
                    "Resource::<PipelineArtifactBucketD127CCF6.Arn>/*"
                ])
            ]
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/SelfMutatePolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Resource::*",
                    "Resource::arn:aws:iam::*:role/cdk-*",
                    "Action::codeartifact:List*"
                ])
            ]
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/BuildSamplePackage/sample-package/Role/DefaultPolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Resource::arn:<AWS::Partition>:logs:<AWS::Region>:<AWS::AccountId>:log-group:/aws/codebuild/<BuildSamplePackagesamplepackageB2962058>:*",
                    "Resource::arn:<AWS::Partition>:codebuild:<AWS::Region>:<AWS::AccountId>:report-group/<BuildSamplePackagesamplepackageB2962058>-*",
                    "Action::s3:GetBucket*",
                    "Action::s3:GetObject*",
                    "Action::s3:List*",
                    "Action::kms:GenerateDataKey*",
                    "Action::kms:ReEncrypt*",
                    "Resource::<PipelineArtifactBucketD127CCF6.Arn>/*"
                ])
            ]
        )

        NagSuppressions.add_resource_suppressions_by_path(
            self,
            "/CodeartifactStack/BuildSamplePackage/PublishPolicy/Resource",
            [
                NagPackSuppression(id="AwsSolutions-IAM5", reason="Defined by a default policy", applies_to=[
                    "Resource::*",
                    "Action::codeartifact:List*"
                ])
            ]
        )

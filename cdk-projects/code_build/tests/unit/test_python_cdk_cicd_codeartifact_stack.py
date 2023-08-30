# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0

import aws_cdk as core
from aws_cdk import assertions

from python_cdk_cicd_codeartifact.python_cdk_cicd_codeartifact_stack import PythonCdkCicdCodeartifactStack


def test_codecommit_created():
    app = core.App()
    stack = PythonCdkCicdCodeartifactStack(app, "PythonCdkCicdCodeartifactStack")
    template = assertions.Template.from_stack(stack)

    template.has_resource_properties("AWS::CodeCommit::Repository", {
        "RepositoryName": "PythonSampleRepository"
    })
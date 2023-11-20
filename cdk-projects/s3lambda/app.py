#!/usr/bin/env python3
import aws_cdk as cdk

from s3lambdatrigger.s3lambdatrigger_stack import S3TriggerStack

app = cdk.App()
S3TriggerStack(app, "S3TriggerStack")

app.synth()

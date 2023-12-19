#!/usr/bin/env python3
import aws_cdk as cdk

from s3lambdatrigger.s3lambdatriggerstack import S3lambdaTriggerStack

app = cdk.App()
S3lambdaTriggerStack(app, "S3lambdaTriggerStack")

app.synth()

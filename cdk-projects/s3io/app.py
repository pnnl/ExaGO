#!/usr/bin/env python3
from aws_cdk import App

from s3io.s3iotrigger_stack import S3ioTriggerStack

app = App()
S3ioTriggerStack(app, "S3ioTriggerStack")

app.synth()

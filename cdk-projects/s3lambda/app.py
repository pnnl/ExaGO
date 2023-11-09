#!/usr/bin/env python3
from aws_cdk import App

from s3lambdatrigger.s3lambdatrigger_stack import S3TriggerStack

app = App()
S3TriggerStack(app, "S3ioTriggerStack")

app.synth()

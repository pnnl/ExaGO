#!/bin/bash
chmod +x /aws-lambda/aws-lambda-rie
exec /aws-lambda/aws-lambda-rie $@

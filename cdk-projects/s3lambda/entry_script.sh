#!/bin/sh
if [ -z "${AWS_LAMBDA_RUNTIME_API}" ]; then
  exec /usr/local/bin/aws-lambda-rie /usr/bin/npx aws-lambda-ric
else
  exec /usr/bin/npx aws-lambda-ric
fi

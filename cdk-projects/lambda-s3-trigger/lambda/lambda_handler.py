import exago
import boto3
import subprocess


def main(event, context):
    # save event to logs
    print(event)
# event source, bucket, key
    return {
        'statusCode': 200,
        'body': event
    }

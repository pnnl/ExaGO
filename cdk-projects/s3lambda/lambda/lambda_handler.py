import boto3
import subprocess
import os

import exago
# GOAL : Copy input file to a different location in S3
print(f"Hello world")
# return {'statusCode': 200, 'body': 'Processing complete.'}
s3_client = boto3.client('s3')


def lambda_handler(event, context):
    # Retrieve input file from S3
    exago.initialize("app")
    print(event)
    print('print event')
    input_bucket = event['Records'][0]['s3']['bucket']['name']
    input_file_key = event['Records'][0]['s3']['object']['key']
    # Download input file -> file path should target s3 bucket file.
    local_input_path = f"/tmp/{input_file_key}"
    s3_client.download_file(input_bucket, input_file_key, local_input_path)
    # Process input file using ExaGO (assuming ExaGO processing logic here)
    exago_output = subprocess.check_output(
        ['cat', local_input_path]).decode('utf-8')
    print(exago_output)
    # Upload output to S3 -> Path to
    output_bucket = os.environ['OUTPUT_BUCKET_NAME']
    # Output file path in the output bucket
    output_file_key = 'processed/' + os.path.basename(local_input_path)
    s3_client.put_object(
        Body=exago_output, Bucket=output_bucket, Key=output_file_key)
#
    return {
        'statusCode': 200,
        'body': 'Processing complete.'
    }

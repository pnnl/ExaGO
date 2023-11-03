import boto3
import subprocess
import os

# GOAL : Copy input file to a different location in S3

s3_client = boto3.client('s3')

def lambda_handler(event, context):
    # Retrieve input file from S3
    #input_bucket = event['Records'][0]['s3']['bucket']['exago']
    #input_file_key = event['Records'][0]['s3']['object']['key']

    # Download input file -> file path should target s3 bucket file.
    local_input_path = '/Users/moha907/projects/github-exago/cdk-projects/s3-io/lambda/python_wrapper'
    s3_client.download_file(local_input_path)

    # Process input file using ExaGO (assuming ExaGO processing logic here)
    exago_output = subprocess.check_output(['exago', 's3']).decode('utf-8')

    # Upload output to S3 -> Path to 
    output_bucket = 'output-bucket'
    output_file_key = 'processed/' + os.path.basename(local_input_path)  # Output file path in the output bucket
    s3_client.put_object(Body=exago_output, Bucket=output_bucket, Key=output_file_key)

   # return {
    #    'statusCode': 200,
        #'body': 'Processing complete.'
    #}
    return {
        'statusCode': 200,
        'body': event
    }

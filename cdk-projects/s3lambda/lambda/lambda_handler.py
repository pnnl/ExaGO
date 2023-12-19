import boto3
import os

# Ensure ExaGO and other necessary libraries are imported and used as needed
import exago
# print exago location
s3_client = boto3.client('s3')


def lambda_handler(event, context):
    print(f"Hello world")
    return {'statusCode': 500, 'body': f"exago location: {os.path.abspath(exago.__file__)}"}

    exago.initialize("app")

    # return {'statusCode': 200, 'body': 'Processing complete.'}
    try:
        # Retrieve input file from S3
        input_bucket = event['Records'][0]['s3']['bucket']['name']
        input_file_key = event['Records'][0]['s3']['object']['key']

        # Download input file
        local_input_path = f"/tmp/{input_file_key}"
        s3_client.download_file(input_bucket, input_file_key, local_input_path)

        # Process input file using your logic
        # (e.g., reading file contents, using ExaGO for processing, etc.)

        # Upload processed output to S3
        output_bucket = os.environ['OUTPUT_BUCKET_NAME']
        output_file_key = 'processed/' + os.path.basename(input_file_key)
        with open(local_input_path, 'rb') as file_data:
            s3_client.put_object(
                Body=file_data, Bucket=output_bucket, Key=output_file_key)

        return {'statusCode': 200, 'body': 'Processing complete.'}

    except Exception as e:
        print(e)
        return {'statusCode': 500, 'body': 'Error processing the S3 event.'}

from aws_cdk import (
    aws_lambda as _lambda,
    aws_s3 as _s3,
    aws_s3_notifications,
    aws_iam as iam,
    Stack
)
from constructs import Construct


class S3TriggerStack(Stack):

    def __init__(self, scope: Construct, id: str, **kwargs) -> None:
        super().__init__(scope, id, **kwargs)

        # ecr_repository = ecr.Repository.from_repository_name(self, "ExagoRepo",
        # repository_name="test_repo")  # Specify your repository name

     # Create an S3 bucket for input
        input_bucket = _s3.Bucket(self,
                                  "InputBucket",
                                  bucket_name="s3in6"
                                  )

        # Create an S3 bucket for output
        output_bucket = _s3.Bucket(self,
                                   "OutputBucket",
                                   bucket_name="s3out6"
                                   )

        # Define the Lambda function using a Docker container from ECR
        function = _lambda.DockerImageFunction(self, "ExagoLambdaFunction",
                                               code=_lambda.DockerImageCode.from_image_asset(
                                                   "lambda"),
                                               environment={
                                                   "OUTPUT_BUCKET_NAME": output_bucket.bucket_name
                                               })

        function.add_to_role_policy(iam.PolicyStatement(
            actions=[
                "s3:GetObject",
                "s3:ListBucket"
            ],
            resources=[
                f"{input_bucket.bucket_arn}",
                f"{input_bucket.bucket_arn}/*"
            ]
        ))
        function.add_to_role_policy(iam.PolicyStatement(
            actions=[
                "s3:GetObject",
                "s3:ListBucket",
                "s3:PutObject"
            ],
            resources=[
                f"{output_bucket.bucket_arn}",
                f"{output_bucket.bucket_arn}/*"
            ]
        ))

    # create s3 notification for lambda function
        notification = aws_s3_notifications.LambdaDestination(function)

    # assign notification for the s3 event type (ex: OBJECT_CREATED)
        input_bucket.add_event_notification(
            _s3.EventType.OBJECT_CREATED, notification)

from aws_cdk import (
    aws_lambda as _lambda,
    aws_s3 as _s3,
    aws_s3_notifications,
    aws_ecr as ecr,
    Stack
)
from constructs import Construct


class S3TriggerStack(Stack):

    def __init__(self, scope: Construct, id: str, **kwargs) -> None:
        super().__init__(scope, id, **kwargs)

        ecr_repository = ecr.Repository.from_repository_name(self, "ExagoRepo",
                                                             repository_name="test_repo")  # Specify your repository name

    # Define the Lambda function using a Docker container from ECR
        function = _lambda.DockerImageFunction(self, "ExagoLambdaFunction",
                                               # environment=dict(), # dictionary of environment variables such as LAMBDA_TASK_ROOT and PYTHONPATH
                                               # to specify the location of the ExaGO python wrappers and lambda handler within the Docker image
                                               code=_lambda.DockerImageCode.from_ecr(repository=ecr_repository))

    # create s3 bucket
        s3 = _s3.Bucket(self, "s3bucket")

    # create s3 notification for lambda function
        notification = aws_s3_notifications.LambdaDestination(function)

    # assign notification for the s3 event type (ex: OBJECT_CREATED)
        s3.add_event_notification(_s3.EventType.OBJECT_CREATED, notification)

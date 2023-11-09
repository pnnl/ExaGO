from aws_cdk import (
    aws_lambda as _lambda,
    aws_s3 as _s3,
    aws_s3_notifications,
    Stack
)
from constructs import Construct


class S3TriggerStackExago(Stack):

    def __init__(self, scope: Construct, id: str, **kwargs) -> None:
        super().__init__(scope, id, **kwargs)

        # create lambda function
        # function = _lambda.Function(self, "exago_lambda_function",
        #                           runtime=_lambda.Runtime.PYTHON_3_10,
        #                          handler="lambda_handler.main",
        #
        #                     code=_lambda.DockerImageCode.from_image_asset("./lambda_exago"))
        function = _lambda.DockerImageFunction(self, "ExagoLambdaFunction",
                                               code=_lambda.DockerImageCode.from_image_asset("lambda_exago"))
        # create s3 bucket
        s3 = _s3.Bucket(self, "s3bucket_exago")

        # create s3 notification for lambda function
        notification = aws_s3_notifications.LambdaDestination(function)

        # assign notification for the s3 event type (ex: OBJECT_CREATED)
        s3.add_event_notification(_s3.EventType.OBJECT_CREATED, notification)

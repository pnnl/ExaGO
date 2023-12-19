#!/bin/bash

# Remove all old docker conatiners/images
docker container prune -f &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker container stop &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker rm -f &&
docker images | tail -n +2 | awk '{print $3}' | xargs docker rmi -f &&
# Build docker container locally
docker build -t testcont:latest . --progress=plain

export AWS_ACCESS_KEY_ID="ASIAUOG3HTODAFZBIELW"
export AWS_SECRET_ACCESS_KEY="+E5E2/n8xlN7RUA8qgc+XCY177FTE7Er6hDTKiFG"
export AWS_SESSION_TOKEN="IQoJb3JpZ2luX2VjEAUaCXVzLXdlc3QtMiJIMEYCIQC1J6MGbyg+5Spbf29oprGKESuE0YcTCOMlNlgiejkK9QIhAJwhTFWHTV7BYYw++YggY/g2QtwcITrIPsC8qsJmRsF9KqADCH4QABoMMzA1NDAyNDUyODcwIgwjluErR0VaJVQT8iUq/QKsv2AqwZEYK0u7IAPqVaEo+/InW64+HXK5yMXW5GAXcaQSw+LJm0CLv38TR8uMWNotOam06WvLnG0UrktqIYhYMRYs8kHBV1YSv+Xpcrbb1H2y2NLqBo75qtHLuVdZf/bJxRFmbF5uCFQlCO5ay4SFTrEYAlnOkXnZQh9Si0XfpF3h1DMFb/DGyMqgkAfnMECsXBWBa+iO0GIs0ZPfT0YE0Ub7+G+dCDT3/PKTXi98UQ7LX1LFTRXuT4aaA1GZIdxswA5kgeSMsOfvIkUH2GAOaDW1d3Jr/4j0fld/uastBHLYWKUwIeJQ2MSVxLGQPMlN6nhJRocX5oWsuigzT3c1K9tRqcIE9RhnIXVKhsT+5PTizaK7Gvv800pWGfkqcApRQrVou1AppkczamgWlWsp4DlAnQpKJ7eGGcYx/Adu6d5fZt7CWC5JrkdLBnq2b4fod7QnfTnO+t9qM8lMBxR/4Myq87biA2fWSiY9y91nib3bUc4wP6DqTE/5Oksw0IHzqwY6pQHhAM0NosHOiRd0SKaK9zFg+Q7lKmdLhADAtqhJFzRq8IQFGtsC0m+F1i+gFCEBGcYKihYWHFA3ZzU5LCnPd0Vc3xJaI1D82u6fimOIdVUApmoZANrmfHnxnN1NpHs32OJ7kYvBe3RxTlqkpj+cqgqOLQlb1EfvWphTIXhywZ0nGgZ75k08E07FoPP/RsL05hVSQMegFas5Uhaapn0XzKRRXon1ge0="

# Run the container and test with awslambdarie
docker run -d \
    -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
    -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
    -e AWS_DEFAULT_REGION=us-west-2 \
    -e AWS_SESSION_TOKEN=$AWS_SESSION_TOKEN \
    -v "$HOME/.aws-lambda-rie:/aws-lambda" \
    -v "$PWD/entrypoint.sh:/entrypoint.sh" -p 9000:8080 \
    --entrypoint /bin/bash testcont:latest /entrypoint.sh python -m awslambdaric lambda_handler.lambda_handler &&

#Invoke S3lambda 
curl -XPOST "http://localhost:9000/2015-03-31/functions/function/invocations" -d '{}'

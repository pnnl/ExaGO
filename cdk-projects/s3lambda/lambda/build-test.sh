#!/bin/bash

# Remove all old docker conatiners/images
docker container prune -f &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker container stop &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker rm -f &&
docker images | tail -n +2 | awk '{print $3}' | xargs docker rmi -f &&
# Build docker container locally
docker build -t testcont:latest . --progress=plain &&

export AWS_ACCESS_KEY_ID="ASIAUOG3HTODMKDAPBUM"
export AWS_SECRET_ACCESS_KEY="5gfOpnocoTJ3CRqvfwrRIgbvlXIlU5g4TgVQzoI8"
export AWS_SESSION_TOKEN="IQoJb3JpZ2luX2VjEKr//////////wEaCXVzLXdlc3QtMiJHMEUCIEWvkDCVII9zZfjsqwkGj/Jy3ELxHkfKuDStRNduZY0MAiEAt30M6khnlNEV8RDF88D6PmQq5e+FWKylOnmfWkpzSAEqoAMIMxABGgwzMDU0MDI0NTI4NzAiDFRUcHmdb6iJB/M2hyr9AlvN7BZ1kHcgQ3t5f5yg0JXR7ttc3QVm8O6XV89x7YxSa6lOOnKZjjIqJXQSXFtZsrrYkQAgN2w6ewx10UXZK3PIkgZ7IMeyYNcCJT5oka04rro0JpIXQWQh6CiluFThGVgccWfE5+GuUuec1qnQ0Ia9sXQOE5ImmksPLDd7zO/XVXInplRe3v2m8ttWCaxjXTzWx856lbV9OZlOieKjj2R3QP7So9v02ngRtDVmBMgAtNplxbrAGtNtRnzTCYh1vsitFwTEkJW6a0hq/pMTiL7wL4BiSCpUkvqA6xMbmwMixgtIPXRw6FqAMBFdKwi1kIdIX6qsmn2Yq4wjpHlKngC3w3Xs2WPdWDBgtasxaVkiZcAXa0u4CXjH6j+Kg4zbP2gxjoTVYU3Xgp8sTBvEXvgTntT+zSbny3uBtCUYE2URYFKttXtt9r2W4B1iByGieammi5cNwgeBvwjTRFBUoyYAP/r1+6RnEwFAxqJ3OKM8S5FdCSD//a065ohNyjCcl5esBjqmAe1DDIc1mwMxVFymahiQy/a333aOHAQ+Ly00ot/ZcdlZk/U4esjJog+shYXhGq2JFcPN8ibEZ/5qR57E4v/TjFXg0LgEL9F6QPfR+0QLYJMDW5Mj3D+aT/t2bundZ4Km+YSXxC99tNkA7WErtZZA3bMn8OPyu/m1Savp3XfIBBgbdYIqCXBecEFOpFkRPUX+SmjyB1Fsv+v0BScY7qk5XU+uzMT/OnM="
    -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
    -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
    -e AWS_DEFAULT_REGION=us-west-2 \
    -e AWS_SESSION_TOKEN=$AWS_SESSION_TOKEN \
    -v "$HOME/.aws-lambda-rie:/aws-lambda" \
    -v "$PWD/entrypoint.sh:/entrypoint.sh" -p 9000:8080 \
    --entrypoint /bin/bash testcont:latest /entrypoint.sh python3.11 -m awslambdaric lambda_handler.lambda_handler &&

#Invoke S3lambda 
curl -XPOST "http://localhost:9000/2015-03-31/functions/function/invocations" -d '{}'

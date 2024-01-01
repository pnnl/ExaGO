#!/bin/bash

# Remove all old docker conatiners/images
docker container prune -f &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker container stop &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker rm -f &&
docker images | tail -n +2 | awk '{print $3}' | xargs docker rmi -f &&
# Build docker container locally
docker build -t testcont:latest . --progress=plain &&

export AWS_ACCESS_KEY_ID="ASIAUOG3HTODJ6PGJJWM"
export AWS_SECRET_ACCESS_KEY="naNPkfeHiMrsB62jLLQUq9G2y0ZJVr5nJwfHltAG"
export AWS_SESSION_TOKEN="IQoJb3JpZ2luX2VjEDwaCXVzLXdlc3QtMiJIMEYCIQCZxg4fibd1B9ombyiZ+Yon7ooGnQ4ZYTw67xeZRrbxAgIhAPMwfV74BGkqdl756+VnthLYGiSeXF9jhUqJq6BCZS++KqkDCMX//////////wEQARoMMzA1NDAyNDUyODcwIgxm7YTxqohkRv2Prg8q/QJuKJnj8EBTuWDuzHFD4KporufAWGlM1jUShpxdFeKXQ8QAdiTBpPm3az8kocJtKku2grqyX1W1NoymPFOlmFUF53MpQKXmHD8dxtrC8EIU7K7AfQqEVrL9r5aC6T16QMw4OUzDc4ATNgzkQnK83RlEy2COvetVXuaMZlH3T5/w8BNo/QkIUTxgYC0THQ0xTWLtSjD1vySWKcuzRlKJAZPWLO+uXlu7BJPOe8G0CKU5/3lyZg5u73Y1cvk3Xj6CRrj+Wgt09OYvWbSQsIHMrM5Cy870xXSm/rJAewn0ztmYr8sdfn0YiQh+te66JPjyc3G/zJFQjaqumMQgO49Y2vORUDLp9aEmPvG1Ik3ZmQQJt9Ar7hbRJKtjfaru0PFG+7oPlyIDdzlAuJSjJGDZ0ZNYwj/HzAegFF6AVoIs2ycLnfj/hEXrtePyCZGqFWQRCzGKhP1jCbDUewyim0t9PuEevyoEBJmJGjvzeUE5DlIwpchRd+fwgclpuNobrz4w06y3rAY6pQEtFMltmiCBKx/6e4OAelTETM4VqWiJEI1QEW+s/TdzLUV0Lm7sMEAqtAPow5Fp+VqsxNsVVJNtMProHConRHBA6Jg38Yu2jXpeKhsdDJ8l9vDuygomorJGjNFdPTMJJksrX7beBTDwDEv+/OeTPYaqEaCmWDF7S6SF2Prlx5mTIQEoQgqDyJ8EfFVMmgd11F1i5JDPRrc/YNw1lKv8reNBKd3hpkU="    -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
    -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
    -e AWS_DEFAULT_REGION=us-west-2 \
    -e AWS_SESSION_TOKEN=$AWS_SESSION_TOKEN \
    -v "$HOME/.aws-lambda-rie:/aws-lambda" \
    -v "$PWD/entrypoint.sh:/entrypoint.sh" -p 9000:8080 \
    --entrypoint /bin/bash testcont:latest /entrypoint.sh python3.11 -m awslambdaric lambda_handler.lambda_handler &&

#Invoke S3lambda 
curl -XPOST "http://localhost:9000/2015-03-31/functions/function/invocations" -d '{}'

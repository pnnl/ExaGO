#!/bin/bash

# Remove all old docker conatiners/images
docker container prune -f &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker container stop &&
docker container ls | tail -n +2 | awk '{print  $1}' | xargs docker rm -f &&
docker images | tail -n +2 | awk '{print $3}' | xargs docker rmi -f &&
# Build docker container locally
docker build -t testcont:latest . --progress=plain &&

export AWS_ACCESS_KEY_ID="ASIAUOG3HTODNB3AHHRI" &&
export AWS_SECRET_ACCESS_KEY="+heBIscSdtHGnlKvQDyACrPpVYSOhwnjbyXnyXMS" &&
export AWS_SESSION_TOKEN="IQoJb3JpZ2luX2VjEHgaCXVzLXdlc3QtMiJGMEQCIGtlJsEOmlm+/WfCuIHYHQS7A+Uw4JHsVUvHpd9/nxctAiB7P89sxdtXPgr00ndpxUsAR31hEXBDQ8rceHVBIPGAnCqpAwjx//////////8BEAEaDDMwNTQwMjQ1Mjg3MCIM3rykMlIcEP8Nbs/GKv0Crj7RjI3RzdNZ0c7J7fouXaDD/wOhbOx6p0Zbzh4jw97SjnB9AnJsDT6xCDXorVsP+ppdcvkBEcCBgBNlxFb8Bp38a4CCFSESScI/ricC4QbuZDeN+gtYOFMwbiaYlKiZ1PaKcMplHugTDBhKYG2iP02WIOVbc9mpZ+zZawdwbwA2Zqfv/wNT2dVOeCbErkcZMMxII3c3sFKmzHp0qwbSr99m+oab5p8NLKwVa9Ly2glJeshM5Pfz3ADPrtIUCKtwdltbCb4qSI9rv9TYK6MvgfnhldzAkDE9VlDuQzLDKIq+fGBeMRLc+uymGp+m1zLJZoiDDyYQh9hVzzHkrnV+5U8ESAbjnv86HQN0nwNaTEcX4LwD/6ONu4Z8N4YUCaXC0HhdGkjC5OUDucU6z9lvk2VoukEhhhsmrx+5jUbNKrebFp5udKuzTh3NGqLM8QByO967movw7uXvkjR8hgEqV9gTFJmK2XVjeP2Am2qYn36pNHSuZscy9ljeiOdUMM6PjKwGOqcB1FRysHdjMN+9Ob6nKHZaXBTYL9DRYlRwVYdwe3vgz7pHRaa9hsvuHWl3TQ9gKMOKdxVuem58JSiVVtmJaNOfbOVXNiNa7Supala5dMxN+ZkVaTyrtcvo2WWFGfl5zfK9hSmpHAF3OAdPG4HqcvhZPwvIOSGyP0FRL0tdGfkUhOuaSgNBUSGXOuUDLq8jQmjLeBRF+L3QNojSfbw7MSnzPwtjyWnB86c=" &&
# Run the container and test with awslambdarie
docker run -d \
    -e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
    -e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
    -e AWS_DEFAULT_REGION=us-west-2 \
    -e AWS_SESSION_TOKEN=$AWS_SESSION_TOKEN \
    -v "$HOME/.aws-lambda-rie:/aws-lambda" \
    -v "$PWD/entrypoint.sh:/entrypoint.sh" -p 9000:8080 \
    --entrypoint /bin/bash testcont:latest /entrypoint.sh python3.11 -m awslambdaric lambda_handler.lambda_handler &&

#Invoke S3lambda 
curl -XPOST "http://localhost:9000/2015-03-31/functions/function/invocations" -d '{}'

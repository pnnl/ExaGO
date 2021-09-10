FROM minio/mc

RUN echo "some text" >> file.txt && \
  . /kaniko/s3env.sh && \
  cat /kaniko/s3env.sh && \
  echo $AWS_ACCESS_KEY_ID && \
  mc alias set minio $S3_ENDPOINT_URL $AWS_ACCESS_KEY_ID $AWS_SECRET_ACCESS_KEY && \
  mc --debug --insecure cp ./file.txt minio/spack

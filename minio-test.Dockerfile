FROM minio/mc

ARG AWS_ACCESS_KEY_ID
ARG AWS_SECRET_ACCESS_KEY

ENV S3_ENDPOINT_URL=http://cache.exasgd.pnl.gov

RUN echo "some text" >> file.txt && \
  mc alias set minio $S3_ENDPOINT_URL $AWS_ACCESS_KEY_ID $AWS_SECRET_ACCESS_KEY && \
  mc --debug --insecure cp ./file.txt minio/spack

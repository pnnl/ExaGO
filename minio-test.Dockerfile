FROM minio/mc

ARG _AWS_ACCESS_KEY_ID
ENV AWS_ACCESS_KEY_ID=$_AWS_ACCESS_KEY_ID

ARG _AWS_SECRET_ACCESS_KEY
ENV AWS_SECRET_ACCESS_KEY=$_AWS_SECRET_ACCESS_KEY

ENV S3_ENDPOINT_URL=http://cache.exasgd.pnl.gov

# mc alias set minio $S3_ENDPOINT_URL $AWS_SECRET_ACCESS_KEY $AWS_ACCESS_KEY_ID && \
RUN echo "some text" >> file.txt && \
  mc alias set minio http://cache.exasgd.pnl.gov pXuBGQLYoa syLgFihiDLd90zO9pAOIZrndRvXpG15cY4fOhfgc && \
  mc --debug --insecure cp ./file.txt minio/spack

#!/usr/bin/bash

cat <<EOD
stages:
  - build
EOD

for filename in *.yaml
do
  export environment="${filename%.*}"
  cp $filename spack.yaml
  export dockerfile="Dockerfile.${environment}"
  cat <<EOD
job_${environment}:
  stage: build
  tags: [k8s, ikp, exasgd, basic]
  image: ubuntu:20.04
  script:
    |
    set -x
    echo "Building dockerfile for ${environment}"
EOD
  rm spack.yaml
done

#!/usr/bin/bash

set -x
for filename in *.yaml
do
  export environment="${filename%.*}"
  cp $filename spack.yaml
  export dockerfile="Dockerfile.${environment}"
  echo "spack containerize > '${dockerfile}'"
  cat <<EOD
job_${environment}:
  stage: spack-container-build
  tags: [k8s, ikp, exasgd, basic]
  image: ubuntu:20.04
  script:
    |
    set -x
    echo "Building dockerfile for ${environment}"
EOD
  rm spack.yaml
  echo "rm $dockerfile"
done

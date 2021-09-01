#!/usr/bin/bash

# git clone --single-branch -b develop https://github.com/LLNL/spack.git
source ./spack/share/spack/setup-env.sh
set -e

cat <<EOD
# Include PNNL GitLab stdlib
include:
  - remote: 'https://raw.githubusercontent.com/pnnl-miscscripts/gitlab-lib/v1/gitlab-lib.yaml'

stages:
  - build
EOD

for filename in environments/*.yaml
do
  # Strip extension
  export tmp_filename=$(basename $filename)
  export environment="${tmp_filename%.*}"

  # Create dockerfile for spack environment
  cp $filename spack.yaml
  export dockerfile="Dockerfile.${environment}"
  spack containerize > $dockerfile

  # Generate gitlab-ci job for container build
  cat <<EOD
${environment}-build:
  stage: build
  tags: [k8s, ikp, exasgd, basic]
  extends:
    - .pnnllib-gitlab-build-container-image
  variables:
    CONTAINER_TAG: $environment
  before_script:
    - cp buildsystem/container/$dockerfile Dockerfile
  needs:
    - pipeline: \$PARENT_PIPELINE_ID
      job: spack-generate-job
EOD

  # Remove temp spack environment, keep dockerfile as an artifact
  rm spack.yaml
done

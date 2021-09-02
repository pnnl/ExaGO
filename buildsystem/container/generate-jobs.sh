#!/usr/bin/bash

cat <<EOD
# Include PNNL GitLab stdlib
include:
  - remote: 'https://raw.githubusercontent.com/pnnl-miscscripts/gitlab-lib/v1/gitlab-lib.yaml'

stages:
  - build
EOD

find environments -type d -depth 1 | while read -r environment
do
  name=$(basename $environment)
  # Generate gitlab-ci job for container build
  cat <<EOD
${name}-build:
  stage: build
  tags: [k8s, ikp, exasgd, basic]
  extends:
    - .pnnllib-gitlab-build-container-image
  before_script:
    - cp buildsystem/container/${environment}/Dockerfile Dockerfile
  needs:
    - pipeline: \$PARENT_PIPELINE_ID
      job: spack-generate-job
EOD
done

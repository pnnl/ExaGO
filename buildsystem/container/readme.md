## Container Build System

This directory contains a script which generates a gitlab ci pipeline.
This file is read by gitlab to create a dynamic pipeline for each file in the
`environments` subdirectory. To add a spack environment to the build system,
just add a `.yaml` file to the `environments` subdirectory and it will be picked
up by the job.

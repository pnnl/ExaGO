# Using ExaGO inside a VSCode Dev Container

## Prerequisites

1. install Docker Desktop and launch the app
1. install the "Remote Development" extension in VSCode
1. open your local clone of exago in VSCode
1. copy the CoinHSL 2019 tarball to the top level of your exago directory

## Generate the container

Spack is used to generate the base image, with some sed/echo commands used to fine-tune Dockerfile for use case within a Dev Container. Run `.devcontainer/create_dockerfile.sh` in order to re-generate the Dockerfile. This should be done semi-frequently whenever a spack version or exago version needs updating.

Note that you should not run this script from within a working devcontainer as it might be quite slow.

## Build Container

The build info for this container is in `.devcontainer/`. There is a Dockerfile and json file associated with the configuration.

1. if connected, disconnect from the PNNL VPN
1. launch the container build  
    * `cmd shift p` to open the command pallette in vscode
    * click `> Dev Container: rebuild and reopen container`
    * this will start building the container, taking about 40 minutes
    * click on the pop up with `(show log)` to view the progress

## Run Notebook
1. once the conatiner has build and launched, open `docs/devcontainer/tutorial.ipynb`
1. select the (newly built) existing jupyter kernel "ExaGO"
1. run all cells!

## Devcontainer Quickstart

A devcontainer is configured through three things here:
- `create_dockerfile.sh` generates the Dockerfile, and needs to be run from ExaGO root with spack submodule cloned with `$ .devcontainer/create_dockerfile.sh`. Note that this currently uses sed commands that only work on Mac shell... so hopefully the Dockerfile is easy enough to fix if slightly broken, but ideally we also have CI for this.
- `Dockerfile` which is generated from the above bash script, and defines a container with ExaGO. Currently this build takes 1.5 hours, but ideally this is cached and successfully pulled from CI builds.
- `devcontainer.json` is unfortunately a json that forbids comments, however should be self-explanatory. You can define extensions to load here which is really helpful for common environments, and you also mount the local folder so you can play with ExaGO clone in container easily.

Note - pushing/pulling from git is not supported in a devcontainer, and should be done independently.

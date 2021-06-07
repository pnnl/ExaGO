## Containerized ExaGO Development


### Building Yourself

To get started, assuming you already have docker installed:

```console
$ cd exago/buildsystem/container
$ docker volume create exago-devel
$ docker build -t exago-devel:latest .

$ # This command will drop you into the docker contianer
$ docker run \
  --mount source=exago-devel,target=/workspace \
  -it exago-devel bash

[root@77929c3c7a65 /]# spack env list
==> 1 environment
    exago-devel
[root@77929c3c7a65 /]# spack env activate exago-devel
```

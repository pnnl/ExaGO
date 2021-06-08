## Containerized ExaGO Development

### Using the Image

```console
$ # Use your PNNL GitLab login
$ docker login gitlab.pnnl.gov:4567

$ # Pull + run the image
$ docker pull gitlab.pnnl.gov:4567/exasgd/frameworks/exago/with-repos:latest
$ docker tag gitlab.pnnl.gov:4567/exasgd/frameworks/exago/with-repos:latest exago
$ docker run -it exago bash
Loading spack environment exago-devel...
Setting environment variables EXAGO_SRC_DIR, EXAGO_BUILD_DIR, and EXAGO_INSTALL_DIR...
Setting environment variables HIOP_SRC_DIR, HIOP_BUILD_DIR, and HIOP_INSTALL_DIR...

# Now you're able to use ExaGO and HiOp installations
[root@b7ffc29b7dcf ~]# cd $EXAGO_INSTALL_DIR
[root@b7ffc29b7dcf install-exago]# ls
bin  include  lib  share  testAcopf  testErrorHandler  tests
[root@b7ffc29b7dcf install-exago]# ./bin/opflow -netfile ./share/exago/datafiles/case9/case9mod.m -opflow_solver HIOP -opflow_model PBPOLRAJAHIOP
[ExaGO INFO]: -options_file not passed.
[ExaGO INFO]: -- Checking /root/install-exago/share/exago/options/opflowoptions exists: yes
[ExaGO INFO]: Creating OPFlow
[ExaGO INFO]: Finalizing opflow application.
```

### Building the Images Yourself

This documentation assumes you already have docker installed.

#### Image with HiOp and ExaGO Pre-Built

```console
$ cd /path/to/exago
$ cd buildsystem/container/with-repos
$ docker build \
  -t gitlab.pnnl.gov:4567/exasgd/frameworks/exago/with-repos:devel .
$ docker push gitlab.pnnl.gov:4567/exasgd/frameworks/exago/with-repos:devel
```

Note: this image depends on the `deps-only` image.

#### Dependency-Only Image

```console
$ cd /path/to/exago
$ cd buildsystem/container/deps-only

$ # This container build requires COINHSL!
$ cp /path/to/coinhsl-archive-2015.06.23.tar.gz ./coinhsl-archive-2015.06.23.tar.gz

$ docker build \
  -t gitlab.pnnl.gov:4567/exasgd/frameworks/exago/deps-only:devel .
$ docker push gitlab.pnnl.gov:4567/exasgd/frameworks/exago/deps-only:devel
```

Note: this image depends on the `exago-base` image.

#### Image with only the OS packages installed

```console
$ cd /path/to/exago
$ cd buildsystem/container/exago-base
$ docker build \
  -t gitlab.pnnl.gov:4567/exasgd/frameworks/exago/exago-base:devel .
$ docker push gitlab.pnnl.gov:4567/exasgd/frameworks/exago/exago-base:devel
```

Note: all other ExaGO images are built on this image. Changing this image
requires rebuilding all other images.

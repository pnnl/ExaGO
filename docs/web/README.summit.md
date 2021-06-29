## Build with everything enabled

The following commands will build ExaGO with most options enabled on Summit:
```console
$ source ./buildsystem/gcc-cuda/summitVariables.sh
$ mkdir build
$ cd build
$ cmake -C ../buildsystem/gcc-cuda/cache.cmake ..
$ make -j 12 install
$ make test
```

The script `./buildsystem/gcc-cuda/summitVariables.sh` will load all the needed modules
to build ExaGO on Summit as of 6/23/2021. The system modules change somewhat
frequently on Summit, so if some modules are not avialable but you need to build
there, please contact Asher Mancinelli <asher.mancinelli@pnnl.gov> or better yet
file an issue here: https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/issues.

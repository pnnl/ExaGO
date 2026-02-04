## Build with everything enabled

The following commands will build ExaGO with most options enabled on Frontier:
```console
$ source ./buildsystem/clang-hip/frontierVariables.sh
$ mkdir build
$ cd build
$ cmake -C ../buildsystem/clang-hip/cache.cmake ..
$ make -j 12 install
$ make test
```

The script `./buildsystem/clang-hip/frontierVariables.sh` will load all the modules
needed to build ExaGO on Frontier. The system modules change somewhat
frequently on Frontier, so if some modules are not avialable but you need to build
there, please file a Github issue.

Alternatively, the `buildsystem/build.sh` can be used as follows to build. 
See the script's documentation for more details.
```console
$ ./buildsystem/build.sh --job=clang-hip --build-only
```

## HiOp

### Download
Download HiOp from its' [github site](https://github.com/LLNL/hiop) and switch to `dev/NewtonMDS` branch.
```
git clone https://github.com/LLNL/hiop.git
cd hiop
git checkout dev/NewtonMDS
```

### Install
Refer to installation instructions to build and compile HiOp

#### Set environment variables (only with `make` build)
SCOPFLOW needs to know the location of the HiOp directory to correctly include its header files and link with the libraries. Set the environment variable HIOP_BUILD_DIR to point to the location of HiOp's build directory. Note that this is not the location of IPOPT's source code directory, it is the location where IPOPT was built.
```
export HIOP_BUILD_DIR=<location-of-HIOP-build-dir>
```
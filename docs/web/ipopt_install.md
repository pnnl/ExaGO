## IPOPT

IPOPT is a popular nonlinear interior point method package for solving general nonlinear optimization problems. Like PIPS-NLP, it is written in C++ and has a C interface available.

### Download
Visit the IPOPT site `https://github.com/coin-or/Ipopt` for download instructions for IPOPT. Note that IPOPT has several dependencies (in particular the HSL libraries) that need to be downloaded and put in the correct IPOPT folders before installing IPOPT. IPOPT tarballs are available at `https://www.coin-or.org/download/source/Ipopt/`

### Install
Refer to the IPOPT installation instructions at `https://coin-or.github.io/Ipopt/INSTALL.html`

#### Set environment variables (only with `make` build)
SCOPFLOW needs to know the location of the IPOPT directory to correctly include its header files and link with the libraries. Set the environment variable IPOPT_BUILD_DIR to point to the location of IPOPT's build directory. Note that this the not the location of IPOPT's source code directory, it is the location where IPOPT was built.
```
export IPOPT_BUILD_DIR=<location-of-IPOPT-build-dir>
```
# Spack Installation

Spack is a package manager for installing packages from source and maintaining many versions of many packages.

Spack may be installed with:
```console
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
export PATH=$PWD/spack/bin:$PATH
spack install zlib
```

The rest of this documentation will assume you have the `spack` executable on your path.

More information on Spack can be found at the [GitHub repo](https://github.com/spack/spack) and on their [readthedocs page](https://spack.readthedocs.io/en/latest/).

## Getting Started

The ExaGO Spack package was first avilable in spack version `0.16.1-2103-c4a83aa22c`.
To use the ExaGO Spack package, download spack and check out a recent version like so:

```console
git clone https://github.com/spack/spack.git
export PATH="$PATH:$PWD/spack/bin"
```

## Building with Spack

The following configuration spack configuration is reccommended for a minimal install with Ipopt solver.
The command `spack info exago` will provide the most up-to-date information about the ExaGO spack package.

***NOTE:*** Please see [the section on installing CoinHSL](#coinhsl-dependency) if you would like to build with Ipopt.
CoinHSL requires a separate license and download process.

```console
source /path-to-spack/share/spack/setup-env.sh
spack compiler find
spack install exago@develop%gcc \
  ^openmpi ^ipopt@3.12.10+coinhsl~mumps ^coinhsl+blas \
  ^petsc@3.13.6+mpi~hypre~superlu-dist~mumps+shared
spack load exago
opflow -help
```

## Installation on OSX

The simplest version of ExaGO running on OSX uses GCC 8-10 and Ipopt with the COINHSL solver library:

```console
uname -a
Darwin WE40281 19.6.0 Darwin Kernel Version 19.6.0: Tue Jan 12 22:13:05 PST 2021; root:xnu-6153.141.16~1/RELEASE_X86_64 x86_64

spack install -j 16 exago%gcc@10.2.0+ipopt ^ipopt+coinhsl~mumps

spack load -r exago@1.0.0
opflow -help
===================================================================
ExaGO Version Info:

ExaGO version 1.0.0 released on 2021-04-07
built with:
	PETSC                                YES
	MPI                                  YES
	Ipopt                                YES
	HiOp                                  NO
	GPU                                   NO
	RAJA                                  NO
============== Help Options for application opflow ==============

 General usage: mpiexec -n <N> ./opflow <options>
 Options:
	 -netfile <netfilename>
	 -opflow_model <POWER_BALANCE_POLAR|...>
   ...
```

The installation prefix may be found with `spack location -i exago`.
See the [Spack documentation](https://spack.readthedocs.io/) for more information.

# FAQ/Common Issues

## Python on OSX

Especially on OSX, we prefer using GCC. This can be an issue when installing ExaGO since ExaGO relies on PETSc, which relies on Python, which does not compile on OSX with GCC.
You may get an error message like this if you attempt to do so:

```shell
==> Error: Conflicts in concretized spec "exago@1.0.0%gcc@10.2.0..."
List of matching conflicts for spec:

    python@3.8.8%gcc@10.2.0...

1. "%gcc platform=darwin" conflicts with "python" [CPython does not compile with GCC on macOS yet, use clang. See: https://github.com/python/cpython/pull/13306]
```

The best way around this is usually to compile python yourself or use Brew/Port to install python and inform spack of where the package is:

```shell
$ brew info python
python@3.9: stable 3.9.2 (bottled)
Interpreted, interactive, object-oriented programming language
https://www.python.org/
/usr/local/Cellar/python@3.9/3.9.2_1 (3,936 files, 66.2MB) *    <--- this is what we're looking for
  Poured from bottle on 2021-03-05 at 18:01:22
From: https://github.com/Homebrew/homebrew-core/blob/HEAD/Formula/python@3.9.rb
License: Python-2.0

$ # The following command may find many other preinstalled packages as well.
$ PATH="$PATH:/usr/local/Cellar/python@3.9/3.9.2_1" spack external find
==> The following specs have been detected on this system and added to /Users/manc568/.spack/packages.yaml
python@3.9.2
```

After finding this external python, spack will no longer complain about the GCC/Python OSX conflict.
This strategy may also be used for other packages that are more difficult to install, such as OpenMPI or PETSc.

## CoinHSL Dependency

Some dependencies in Spack are contingent on the user obtaining a proper licence for the software. 
That is currently the case for `CoinHSL`.
[More information on downloading CoinHSL can be found at this link](https://www.hsl.rl.ac.uk/ipopt/).
In order to activate CoinHSL, add the corresponding variant to HiOp (`+sparse`) and
Ipopt (`+coinhsl`) and run `spack install` with the archive `coinhsl*tar.gz` in
the current working directory.

Sometimes, the CoinHSL tarball does not have the name the spack package expects it to have.
If this is the case, you may run into error messages like `all fetchers failed!`.
In this case, please rename the CoinHSL tarball to `coinhsl-archive-YYYY.MM.DD.tar.gz`.

If you view the CoinHSL spack packge (with the command `spack edit coinhsl`), you'll find the line:
```python
    url = "file://{0}/coinhsl-archive-2014.01.17.tar.gz".format(os.getcwd())
```

This `url` variable determines the pattern spack uses to find the CoinHSL tarball.
As long as your tarball matches this pattern, spack should be able to find the package.

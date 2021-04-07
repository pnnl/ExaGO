# Spack Installation

## Getting Started

The ExaGO Spack package was first avilable in spack version `0.16.1-2103-c4a83aa22c`.
To use the ExaGO Spack package, download spack and check out a recent version like so:

```shell
$ git clone https://github.com/spack/spack.git
$ export PATH="$PATH:$PWD/spack/bin"

$ # The exact version may differ, but the full version will be at least 0.16.1
$ spack --version
0.16.1-2103-7c87ebeb91

$ spack info exago
CMakePackage:   exago

Description:
    ExaGO is a package for solving large-scale power grid optimization
    problems on parallel and distributed architectures, particularly
    targeted for exascale machines.

Homepage: https://gitlab.pnnl.gov/exasgd/frameworks/exago

Tags:
    None

Preferred version:
    1.0.0     [git] https://gitlab.pnnl.gov/exasgd/frameworks/exago.git at tag v1.0.0

Safe versions:
    master    [git] https://gitlab.pnnl.gov/exasgd/frameworks/exago.git on branch master
    1.0.0     [git] https://gitlab.pnnl.gov/exasgd/frameworks/exago.git at tag v1.0.0
    0.99.2    [git] https://gitlab.pnnl.gov/exasgd/frameworks/exago.git at tag v0.99.2
    0.99.1    [git] https://gitlab.pnnl.gov/exasgd/frameworks/exago.git at tag v0.99.1

Variants:
    Name [Default]                 Allowed values          Description
    ===========================    ====================    ==================================

    build_type [RelWithDebInfo]    Debug, Release,         CMake build type
                                   RelWithDebInfo,
                                   MinSizeRel
    cuda [off]                     on, off                 Build with CUDA
    cuda_arch [none]               none, 12, 20, 35,       CUDA architecture
                                   60, 32, 80, 86, 72,
                                   75, 13, 52, 70, 37,
                                   53, 11, 50, 21, 10,
                                   61, 30, 62
    hiop [off]                     on, off                 Enable/Disable HiOp
    ipo [off]                      on, off                 CMake interprocedural optimization
    ipopt [off]                    on, off                 Enable/Disable IPOPT
    mpi [on]                       on, off                 Enable/Disable MPI
    petsc [on]                     on, off                 Enable/Disable PETSc
    raja [off]                     on, off                 Enable/Disable RAJA

Installation Phases:
    cmake    build    install

Build Dependencies:
    blas  camp  cmake  cuda  hiop  ipopt  mpi  petsc  raja  umpire

Link Dependencies:
    blas  camp  cuda  hiop  ipopt  mpi  petsc  raja  umpire

Run Dependencies:
    None

Virtual Packages:
    None
```

## Installation on OSX

The simplest version of ExaGO running on OSX uses GCC 8-10 and Ipopt with the COINHSL solver library:

```shell
$ uname -a
Darwin WE40281 19.6.0 Darwin Kernel Version 19.6.0: Tue Jan 12 22:13:05 PST 2021; root:xnu-6153.141.16~1/RELEASE_X86_64 x86_64

$ spack install -j 16 exago%gcc@10.2.0+ipopt ^ipopt+coinhsl~mumps

$ spack load -r exago@1.0.0
$ opflow -help
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
	 -opflow_solver <IPOPT|...>
	 -opflow_initialization <MIDPOINT|...>
	 -opflow_ignore_lineflow_constraints <0|1>
	 -opflow_include_loadloss_variables <0|1>
	 -opflow_include_powerimbalance_variables <0|1>
	 -opflow_loadloss_penalty <Penalty ($)>
	 -opflow_powerimbalance_penalty <Penalty ($)>
	 -opflow_genbusvoltage <FIXED_WITHIN_QBOUNDS|...>
	 -opflow_has_gensetpoint <0|1>
	 -opflow_objective <MIN_GEN_COST|...>
	 -opflow_use_agc <0|1>
	 -opflow_tolerance <1e-6|...>
	 -hiop_compute_mode <hybrid|...>
	 -hiop_verbosity_level <0-10>
	 -hiop_tolerance <1e-6|...>
	 -print_output <0|1>
	 -save_output <0|1>
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

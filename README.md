[![pipeline status](https://gitlab.pnnl.gov/exasgd/frameworks/exago/badges/master/pipeline.svg)](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/commits/master)

# <b>Exa</b>scale <b>G</b>rid <b>O</b>ptimization toolkit (ExaGO<sup>TM</sup>)
ExaGO<sup>TM</sup> is a package for solving large-scale power grid optimization problems on parallel and distributed architectures, particularly targeted for exascale machines with heteregenous architectures (GPU). It is written in C/C++ using the [PETSc](https://www.mcs.anl.gov/petsc/) library. An overview of the package is given on this page and the different links provided. For detailed information, including the different formulations used, see the [ExaGO manual](docs/manual/manual.pdf). 

ExaGO<sup>TM</sup> includes the following applications for solving different power grid optimization problems:

- [OPFLOW](docs/web/opflow.md) solves an AC optimal power flow either on CPU and GPU
- [TCOPFLOW](docs/web/tcopflow.md) solves a multi-period optimal power flow
- [SCOPFLOW](docs/web/scopflow.md) solves a multi-period security-constrained (contingency-constrained) optimal power
- [SOPFLOW](docs/web/sopflow.md) solves a stochastic security-constrained multi-period optimal power flow

ExaGO<sup>TM</sup> applications can use the following optimization packaages:

- [IPOPT](https://github.com/coin-or/Ipopt) is a popular optimization package for solving nonlinear optimization problems that uses an interior-point algorithm.
- [HiOp](https://github.com/LLNL/hiop) is a HPC package for optimization. ExaGO interfaces with two of its solvers -- a mixed sparse-dense interior-point solver (NewtonMDS) and a sparse interior-point solver (HiOPSparse). NewtonMDS  allows execution of the optimization either on CPU and GPU. The sparse HiOp solver is currently supported on CPU only.

Note that not all applications can utilize all solvers yet. The following table lists the solver-application compatibility.

|  Solver    | OPFLOW  | TCOPFLOW | SCOPLOW | SOPFLOW |
|:------:|:---------:|:-----:|:-------:|:-------:|
| IPOPT      | Y         |  Y     | Y       | Y       | 
| HIOP       | Y          |       |   Y      |  Y       |

Note: The support for solving SCOPFLOW and SOPFLOW in parallel using the HiOP package is currently in development and is expected to be added in the next release.

## Installing

See [INSTALL.md](./INSTALL.md) for information on acquiring, building and installing ExaGO.

## Usage
Instructions for executing the different ExaGO<sup>TM</sup> applications is given below.
- [OPFLOW](docs/web/opflow.md)
- [TCOPFLOW](docs/web/tcopflow.md)
- [SOPFLOW](docs/web/sopflow.md)
- [SCOPFLOW](docs/web/scopflow.md)
- [PFLOW](docs/web/pflow.md)

### Options

Each application has a different set of options that are described in depth in the usage notes. These options can be passed optionally through an options file (`-options_file <option_file>`), or directly on the command line.

Since options may be specified in more than one location (on the command line, and through an options file), it is worth noting that the option specified **last** is used. For example, if `opflowoptions` specified `-netfile case9mod.m`, the following behavior occurs:

```bash
# This uses case118.m
./opflow -options_file opflowoptions -netfile case118.m

# This uses case9mod.m
./opflow -netfile case118.m -options_file opflowoptions
```

Note that all ExaGO applications must run with an options file passed, and so if none is specified on the command line, ExaGO attempts to use the default application options in the `options` directory. 

## Contributing

Please see [the developer guidelines](docs/DeveloperGuidelines.md) before attempting to contribute.
Feel free to raise an issue or contact the team if the guidelines are ambiguous or you have a particular question.

## Authors
- Shrirang Abhyankar
- Slaven Peles
- Asher Mancinelli
- Cameron Rutherford
- Bruce Palmer

## Acknowledgement
This package is developed as a part of [ExaSGD](https://www.exascaleproject.org/research-project/exasgd/) project under the [Exascale computing project](https://www.exascaleproject.org/).

## Copyright
Copyright &copy; 2020, Battelle Memorial Institute.

ExaGO<sup>TM</sup> is a free software distributed under a BSD 2-clause license. You may reuse, modify, and redistribute the software. See the [license](LICENSE) file for details.


## Disclaimer
This material was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

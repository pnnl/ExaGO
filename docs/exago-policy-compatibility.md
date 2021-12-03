# xSDK Community Policy Compatibility for ExaGO

This document summarizes the efforts of current and future xSDK member packages to achieve compatibility with the xSDK community policies. Below only short descriptions of each policy are provided. The full description is available [here](https://github.com/xsdk-project/xsdk-community-policies)
and should be considered when filling out this form.

Please, provide information on your compability status for each mandatory policy, and if possible also for recommended policies.
If you are not compatible, state what is lacking and what are your plans on how to achieve compliance.

**Website:** https://gitlab.pnnl.gov/exasgd/frameworks/exago

### Mandatory Policies

| Policy                 |Support| Notes                   |
|------------------------|-------|-------------------------|
|**M1.** Support portable installation through Spack. |Full|ExaGO Spack package is available in the official Spack repository. The packages are continuously updated.|
|**M2.** Provide a comprehensive test suite for correctness of installation verification. |Full|Comprehensive test suite with 21 unit tests and 94 integration/funcitonality tests (depending on the branch). See the `tests` directory for the full test suite. `tests/unit` contains unit tests, `tests/functionality` contains full end-to-end tests of the ExaGO applications libraries (eg OPFLOW and SOPFLOW). `tests/interfaces` contains tests for the Python bindings, which are a work in progress.|
|**M3.** Employ user-provided MPI communicator (no `MPI_COMM_WORLD`). Don't assume a full MPI 3 implementation without checking. Provide an option to prevent any changes to MPI error-handling if it is changed by default. |Full| ExaGO is built on PETSc, and uses PETSc's APIs for interacting with MPI. All application structures (eg OPFLOW) are constructed with an MPI communicator, and `MPI_COMM_WORLD`. `grep`ing for `MPI_COMM_WORLD` in ExaGO's repository will identify the test drivers, example application drivers, and an initialization utility as using `MPI_COMM_WORLD`. The utility (`src/utils/utils.cpp`) will only use the global communicator if no communicator is given, and the drivers use the global communicator as examples.|
|**M4.** Give best effort at portability to key architectures (standard Linux distributions, GNU, Clang, vendor compilers, and target machines at ALCF, NERSC, OLCF). |Full| Continuous integration tests each branch on multiple platforms, including IBM Power9 at ORNL and PNNL, and x86 at PNNL. CI on an AMD platform is in progress. |
|**M5.** Provide a documented, reliable way to contact the development team. |Full|[Submit issues on GitLab page linked here.](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/issues) The SUPPORT file in the top-level directory points users to this location as well.|
|**M6.** Respect system resources and settings made by other previously called packages (e.g. signal handling). |Full| No signal hanlders are overridden. |
|**M7.** Come with an open source (BSD style) license. |Full|See LISENCE file in root of source directory. ExaGO uses PNNL/Battelle's BSD-style license.|
|**M8.** Provide a runtime API to return the current version number of the software. |Full|Header `exago_config.h` defines version and configuration information, and we expose various runtime APIs for software information, such as `ExaGOVersionGetVersion` and `ExaGOVersionGetReleaseDate`.|
|**M9.** Use a limited and well-defined symbol, macro, library, and include file name space. |Full| All macros are prefixed with `EXAGO_` and headers installed under `exago/` directory. |
|**M10.** Provide an xSDK team accessible repository (not necessarily publicly available). |Full|[Public GitLab repository linked here](https://gitlab.pnnl.gov/exasgd/frameworks/exago/). |
|**M11.** Have no hardwired print or IO statements that cannot be turned off. |Full|Logging may be disabled with `ExaGOLogSetMinLogLevel(EXAGO_LOG_DISABLE)` or by enabling the CMake option `EXAGO_DISABLE_LOGGING` to ensure the logger is fully disabled.|
|**M12.** For external dependencies, allow installing, building, and linking against an outside copy of external software. |Full| Our CMake configuration allows for discovery of external libraries. We use Spack to install and manage our dependencies, which informs our CMake-based build system of external packages.|
|**M13.** Install headers and libraries under \<prefix\>/include and \<prefix\>/lib. |Full| Build targets are installed using CMake under the recommended directories. Search the top-level `CMakeLists.txt` file for the `include` keyword to view the directories under the installation prefix where build targets are installed.|
|**M14.** Be buildable using 64 bit pointers. 32 bit is optional. |Full| We build using only 64 bit pointers. |
|**M15.** All xSDK compatibility changes should be sustainable. |Full| All changes described in this document have been merged into the key development branches, and do not exist solely in a branch. |
|**M16.** Any xSDK-compatible package that compiles code should have a configuration option to build in Debug mode. |Full| Setting the CMake option `CMAKE_BUILD_TYPE` to `Debug` will build our codebase with debugging symbols. Setting the Spack variant `build_type` to `Debug` will also toggle this option. |

M1 details <a id="m1-details"></a>: optional: provide more details about approach to addressing topic M1.

M2 details <a id="m2-details"></a>: optional: provide more details about approach to addressing topic M2.

### Recommended Policies

| Policy                 |Support| Notes                   |
|------------------------|-------|-------------------------|
|**R1.** Have a public repository. |Full| [Public GitLab repository linked here](https://gitlab.pnnl.gov/exasgd/frameworks/exago/). |
|**R2.** Possible to run test suite under valgrind in order to test for memory corruption issues. |Full| It is possible to run any of the application drivers and test drivers under Valgrind. This has only been test with the leakcheck tool, and not any of the other tools from Valgrind. |
|**R3.** Adopt and document consistent system for error conditions/exceptions. |Full| ExaGO makes thorough use of return codes and error checking, particularly the PETSc macros such as `ERRCHKQ`. |
|**R4.** Free all system resources acquired as soon as they are no longer needed. |Full| Memory for the model is allocated at the beginning of the program and freed at the end. ExaGO also allows for using an external solver, in which case the memory is used by a thrid-party library. These libraries (PETSc, Ipopt, and HiOp) also adequately free memory they allocate. |
|**R5.** Provide a mechanism to export ordered list of library dependencies. |Full| ExaGO exposes two arrays, `ExaGODependencyNames` and `ExaGOIsDependencyEnabled`, allowing users to query dependency information. Only key dependencies are tracked in these arrays, such as RAJA and GPU-related dependencies. |
|**R6.** Document versions of packages that it works with or depends upon, preferably in machine-readable form.  |Full| Our Spack packages document much of this information. Documentation in [`INSTALL.md`](INSTALL.md) and [`docs/InstallingWithSpack.md`](docs/InstallingWithSpack.md) contain additional information about dependencies.|
|**R7.** Have README, SUPPORT, LICENSE, and CHANGELOG files in top directory.  |Full| We currently have README.md, CHANGELOG.md, SUPPORT.md, and LICENSE files in root directory. |
|**R8.** Each xSDK member package should have sufficient documentation to support use and further development.  |Full| The directory `docs/manual` contains thorough documentation in LaTeX with a prebuilt user manual PDF [linked here](docs/manual/manual.pdf). The file [`docs/DeveloperGuidelines`](./docs/DeveloperGuidelines.md) contains documentation on software development best practices that contributors are expected to follow. `docs/web` contains markdown documentation on each of the application libraries and further documentation on some dependencies and platforms. `docs/petsc-dependencies` contains further documentation on PETSc usage. |

# xSDK Community Policy Compatibility for ExaGO

This document summarizes the efforts of current and future xSDK member packages to achieve compatibility with the xSDK community policies. Below only short descriptions of each policy are provided. The full description is available [here](https://github.com/xsdk-project/xsdk-community-policies)
and should be considered when filling out this form.

Please, provide information on your compability status for each mandatory policy, and if possible also for recommended policies.
If you are not compatible, state what is lacking and what are your plans on how to achieve compliance.

**Website:** https://gitlab.pnnl.gov/exasgd/frameworks/exago

### Mandatory Policies

| Policy                 |Support| Notes                   |
|------------------------|-------|-------------------------|
|**M1.** Support portable installation through Spack. |Full|Spack packages are not publicly available yet, but are being reviewed by community members in Spack PRs.|
|**M2.** Provide a comprehensive test suite for correctness of installation verification. |Full|Comprehensive test suite with 21 unit tests and 94 integration/funcitonality tests (depending on the branch)|
|**M3.** Employ user-provided MPI communicator (no MPI_COMM_WORLD). Don't assume a full MPI 3 implementation without checking. Provide an option to prevent any changes to MPI error-handling if it is changed by default. |Full| None. |
|**M4.** Give best effort at portability to key architectures (standard Linux distributions, GNU, Clang, vendor compilers, and target machines at ALCF, NERSC, OLCF). |Full| Continuous integration tests each branch on multiple platforms. |
|**M5.** Provide a documented, reliable way to contact the development team. |Full|[Submit issues on GitLab page linked here.](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/issues)|
|**M6.** Respect system resources and settings made by other previously called packages (e.g. signal handling). |Full| None. |
|**M7.** Come with an open source (BSD style) license. |Full|See LISENCE file in root of source directory.|
|**M8.** Provide a runtime API to return the current version number of the software. |Full|Header `exago_config.h` defines version and configuration information, and we expose various runtime APIs for software information, such as `ExaGOVersionGetVersion` and `ExaGOVersionGetReleaseDate`.|
|**M9.** Use a limited and well-defined symbol, macro, library, and include file name space. |Full| All macros are prefixed with `EXAGO_` and headers installed under `exago/` directory. |
|**M10.** Provide an xSDK team accessible repository (not necessarily publicly available). |Full|[GitLab linked here](https://gitlab.pnnl.gov/exasgd/frameworks/exago/). |
|**M11.** Have no hardwired print or IO statements that cannot be turned off. |Partial|Logging may be disabled with `ExaGOLogSetMinLogLevel(EXAGO_LOG_DISABLE)`.|
|**M12.** For external dependencies, allow installing, building, and linking against an outside copy of external software. |Full| Our CMake configuration allows for discovery of external libraries. |
|**M13.** Install headers and libraries under \<prefix\>/include and \<prefix\>/lib. |Full| None. |
|**M14.** Be buildable using 64 bit pointers. 32 bit is optional. |Full| None. |
|**M15.** All xSDK compatibility changes should be sustainable. |Full| None. |
|**M16.** Any xSDK-compatible package that compiles code should have a configuration option to build in Debug mode. |Full| None. |

M1 details <a id="m1-details"></a>: optional: provide more details about approach to addressing topic M1.

M2 details <a id="m2-details"></a>: optional: provide more details about approach to addressing topic M2.

### Recommended Policies

| Policy                 |Support| Notes                   |
|------------------------|-------|-------------------------|
|**R1.** Have a public repository. |Full| None. |
|**R2.** Possible to run test suite under valgrind in order to test for memory corruption issues. |Full| None. |
|**R3.** Adopt and document consistent system for error conditions/exceptions. |Full| ExaGO makes thorough use of return codes and error checking. |
|**R4.** Free all system resources acquired as soon as they are no longer needed. |Full| None. |
|**R5.** Provide a mechanism to export ordered list of library dependencies. |Full| ExaGO exposes two arrays, `ExaGODependencyNames` and `ExaGOIsDependencyEnabled`, allowing users to query dependency information. |
|**R6.** Document versions of packages that it works with or depends upon, preferably in machine-readable form.  |Partial| Our Spack packages document much of this information. |
|**R7.** Have README, SUPPORT, LICENSE, and CHANGELOG files in top directory.  |Partial| We currently have README.md and LICENSE files in top directory. |
|**R8.** Each xSDK member package should have sufficient documentation to support use and further development.  |Full| None. |

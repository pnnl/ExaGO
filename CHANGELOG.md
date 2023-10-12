# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [develop]

### General
- Transition from GitLab to GitHub
- ExaGO now depends on HiOp 1.0.0. Next release will require HiOp 1.0.1+

### Build system
- Add CI on Incline with Spack module pipelines (AMD / HIP GPU machine at PNNL)
- Add Spack module pipelines to Ascent

### PS

### PFLOW

### OPFLOW

### TCOPFLOW

### SCOPFLOW

### SOPFLOW

### Documentation

### Testing

### Miscallenous

## [1.5.1]

### General
- Removed TAO solver support
- Various fixes for matpower input and PSSE/PowerWorld contingency files
- Support minimal builds without HiOp/Ipopt for OPFLOW/PFLOW only builds
- Added preliminary visualization frontend to ExaGO

### Build system
- Moved to classes for Python wrapper development to enable builds without all solvers enabled
- Allowed disabling of spdlog
- Fixed CMake bugs for MacOS and builds without RAJA.
- Disable SCOPFLOW/SOPFLOW when HiOp/Ipopt is disabled as appropriate.
  
### OPFLOW
- Add minimal output format for OPFLOW solutions
### SCOPFLOW
- Added Python interface for SCOPFLOW.  
- Fixed troubles with Python instance creation and destruction. Destruction of an unsolved SCOPFLOW instance no longer causes a SEGV.  

### SOPFLOW
- Added Python interface for SOPFLOW.

### Documentation
- Updated ExaGO user manual with latest development activities

### Testing
- Enabled Ipopt on all testing platforms (Deception)
- Fully support memory space configuration for HiOp in testing
- Added Ascent pipeline status to GitLab for MRs only, using a push mirror instead of pull mirror
- CI was substantially re-factored, with ability to re-build testing dependencies at PNNL through commit messages.
  
### Miscallenous
- Updated issue and MR templates
- Updated main README with viz and graphics, and fix pipeline status badge

## [v1.5.0]

### General
- OPFLOW, SOPFLOW, and SCOPFLOW can run on AMD GPUs (requires HIP)
- C++ standard changed to 14 in preparation for supporting newest RAJA and Umpire
- PFLOW exclusive build is enabled when building without Ipopt and HiOp
- Python wrapper is temporarily disabled when building without Ipopt and HiOp

### Build system
- Added `EXAGO_ENABLE_LOGGING` flag to enable/disable logging
- Add automatic module generation using Spack to GitLab pipelines for PNNL
- Transition to tracking spack configuration for Newell, Deception and Ascent to repository
- ExaGO now depends on Petsc 3.18.0: and HiOp 0.7.0:
- Remove CMake OpenMP linking when compiling without RAJA
- Refactored GitLab CI scripts into separate files and updated syntax
- Allow building without spdlog


### OPFLOW
- Added OPFLOWSetWeight() for setting weight for OPFLOW. This is done for using with SCOPFLOW and SOPFLOW to set the probabilities or weights for scenarios and contingencies.
- OPFLOW can run on AMD GPUs now.
- The memory space for HIOP can be set with option `-hiop_mem_space`.

### SCOPFLOW
- Support for PTI contingency files. Use .con file extension.

### SOPFLOW
- Read-in weights for scenarios set through the scenario files.
- Solve SOPFLOW with scenario weights

### Documentation
- Added `-hiop_memspace` option

### Testing
- Add tests for OPFLOW Objective function for a known network and reference solution
- Convert to using pytest for Python tests

### Miscallenous
- Add template profiling scripts from SRP HPC Summer Project
- Update Authors

## [v1.4.1]

### Build system
- No longer set the C++ compiler to `hipcc` when building on ROCm platforms.
- Add new CMake option `EXAGO_CTEST_LAUNCH_COMMAND` to allow users to toggle the command used to launch CTest tests

### Documentation
- Document new CMake option `EXAGO_CTEST_LAUNCH_COMMAND`

## [v1.4.0]

### General
- Added mpi4py as a dependency to ExaGO when building with Python and MPI.
- Expanded opflow python wrapper with updated documentation.

### Build system
- Added build system support for finding external python modules.
- Improved option handling mechanisms for Python wrapper.

### PS
- Added API PSCopy() to copy PS object data

### OPFLOW
- DC Optimal power flow (DCOPF) implementation. Can be used as an initialization OR as an OPFLOW model.
- Added API OPFLOWSkipOptions() for skipping options (needed for DCOPF initialization)
- Added select monitor of lines (inequality constraints) for OPFLOW. The selection is currently based by KV levels. The user can provide KV levels to monitor via option -opflow_monitor_kvlevels or OPFLOWSetLinesMonitored() to have OPFLOW include these line flows as inequality constraints. Default is to monitor all lines.
- Added getters for various opflow fields that were previously inaccessible.

### SCOPFLOW
- Fixed bugs in configuring SCOPFLOW subproblems. Added API for SCOPFLOWEnablePowerImbalanceVariables, SCOPFLOWEnablePowerImbalanceVariables, SCOPFLOWSetComputeMode and SCOPFLOWSetComputeMode.

### SOPFLOW
- Read-in weights for scenarios set through the scenario files.
- Solve SOPFLOW with scenario weights

### Testing
- Added testing for python with mpi4py, but still require testing without MPI.
- Re-enabled python tests on Ascent.
- Convert existing smoke tests for SCOPLFOW and SOPFLOW to toml based testing.

## [v1.3.0]

### General
- Python bindings re-released using pybind11 instead of the ctypes bindings.
	This is required to build python bindings when shared libraries are disabled.
- Removed -no_optfile option and searching options file in default path.  Use -optionsfile option to explicitly set the options file.
- Updated help messages with clearer default options

### Build system
- Depend on HiOp's exported cmake configuration to configure HiOp interoperation. Many CMake find package scripts have been removed as a result.
- Use of Pkgconfig instead of PETSc cmake modules.
- ExaGO depends on HiOp v0.5.3 and PETSc v3.16 and above.

### PFLOW
- Native Python bindings to the PFLOW ExaGO library have been added.

### OPFLOW
- Replaced quadratic load loss objective with linear.
- Native Python bindings to the PFLOW ExaGO library have been extended.

### SCOPFLOW
- Defaults to corrective mode of operaation

### SOPFLOW
- Defaults to corrective mode o operation
- Added translator for converting native .cont contingency format to .con PSSE format.	

### Documentation
- Update Python documentation.
- Add release checklist.
- Add and extend developer guidelines.

### Testing
- Improved error messages for functionality tests. Reasons for failure are printed in the TOML output.
- First of a new suite of unit tests has been added.
- Added unit test for new internal logging api.
- Misc fixes for more robust continuous integration.

### Miscallenous
- Updated case9mod_gen3_wind.m to (a) have no reactive power for wind generator, (b) tighter bounds on generator bus voltages.
- Added case9mod_loadloss.m as a test network for non-zero load loss.

### OPFLOW
- Added feature of setting selected lines to monitor 
	- OPFLOWSetMonitoredLines()
	- run-time option -opflow_monitor_line_kvlevels

## [v1.2.0] 

### General
- Added clang format utility script and formatted source code
- Use spdlog for logging instead of hand-rolled logger

### SCOPFLOW
- Added new solver `HIOP` for SCOPFLOW
- New API functions
    - `SCOPFLOWSetSubproblemModel` allows setting of subproblem model. Used with HIOP solver only.
    - `SCOPFLOWSetSubproblemSolver` allows setting of subproblem solver. Used with HIOP solver only.

### SOPFLOW
- Added new solver `HIOP` for SOPFLOW
- New API functions
    - `SOPFLOWSetSubproblemModel` allows setting of subproblem model. Used with HIOP solver only.
    - `SOPFLOWSetSubproblemSolver` allows setting of subproblem solver. Used with HIOP solver only.

### Documentation
- Updated changelog file

### Testing
- Updated tests for OPFLOW, SCOPFLOW, and SOPFLOW

### Miscallenous
- Added a new object `scenariolist` to manage scenarios.
- Updated object `contingencylist` to manage contingencies

## [v1.1.2]

### Documentation
- Patch to bring ExaGO into full compliance with xSDK policies and update our compliance document

## [v1.1.1]

### General
- Hotfix to update versions in user manual and CMake

## [v1.1.0]

### General
- Third-party library compatibiity changes
    - HiOp - 0.5.1
    - Umpire - 0.14.0
    - RAJA - 6.0.0
    - MAGMA - 2.6.1
- Added ExaGOCheckError in lieu of CHKERRQ
- Added ExaGOInitialize and ExaGOFinalize to manage ExaGO initialization and finalization.
- Added python wrapper for ExaGO. Not made public.
- Added developer guidelines `docs/DeveloperGuidelines.md`

### Build system
- Added formatting checks for cmake codebase in CI
- CMake targets exported under `ExaGO::` namespace

### PS
- New API functions
    - PSSet/GetGenStatus - Set/get generator status
    - PSSet/GetLineStatus - Set/get line status
    - PSSett/GeGenPowerLimits - Set/get generator real and reactive power limits
    - PSGetGenDispatch - Get generator dispatch


### OPFLOW
- Modified OPFLOW HIOP interface to use HIOP version 0.5.
- Removed applications opflow-proto.c and opflow-proto2.c
- Added new solver HIOPSPARSE for solving with HiOP sparse solver. Only works on CPU.
- OPFLOWInitializationType enum for managing opflow initialization types.
- New API functions
    - OPFLOWSetHiOPComputMode - Set HiOp compute mode
    - OPFLOWSetHiOpVerbosityLevel - Set HiOp verbosity level
    - OPFLOWHasLoadLoss - Set load loss active
    - OPFLOWHasPowerImbalance - Set power imbalance active
    - OPFLOWSetInitializationType - Set initialization type
    - OPFLOWIgnoreLineFlowConstraints - Activate/Deactivate line flow constraints
    - OPFLOWSetBusPowerImbalancePenalty - Set penalty for power imbalance
    - OPFLOWGetPS - Return the PS object
    - OPFLOWSetUpPS - Set up the PS object
    - OPFLOWSolutionToPS - Transfer OPFLOW solution to PS

### SCOPFLOW
- Added HIOP solver interface. Not public yet.
- New API functions
    - SCOPFLOWSetInitializationType - Set initialization for OPFLOW
    - SCOPFLOWSetGenBusVoltageType - Set generator bus voltage control mode
    - SCOPFLOWSetNumContingencies - Set the number of contingencies
    - SCOPFLOWSetMode - Set the control action mode (preventive or corrective)
- Multiperiod API functions
    - SCOPFLOWSetPload/QloadData - Set data files for Pload/Qload
    - SCOPFLOWSetLoadProfiles - Set the load profiles Pload and Qload files
    - SCOPFLOWSetWindGenProfile - Set wind generation profile file
    - SCOPFLOWSetDuration - Set the duration for multiperiod SCOPFLOW
    - SCOPFLOWSetTimeStepandDuration - Set the time-step and duration

### SOPFLOW
- New API functions
    - SCOPFLOWSetInitializationType - Set initialization for OPFLOW
    - SCOPFLOWSetGenBusVoltageType - Set generator bus voltage control mode
    - SCOPFLOWSetNumContingencies - Set the number of contingencies
    - SCOPFLOWSetMode - Set the control action mode (preventive or corrective)
- Multiperiod API functions
    - SCOPFLOWSetPload/QloadData - Set data files for Pload/Qload
    - SCOPFLOWSetLoadProfiles - Set the load profiles Pload and Qload files
    - SCOPFLOWSetWindGenProfile - Set wind generation profile file
    - SCOPFLOWSetDuration - Set the duration for multiperiod SCOPFLOW
    - SCOPFLOWSetTimeStepandDuration - Set the time-step and duration

### Documentation
- Added top-level INSTALLING.md documentation
- Added file `doc/installing_with_spack.md` with updated spack documentation
- Removed outdated spack documentation
- Removed TAO solver from README and the manual
- Updates to the manual and README for changes to the different ExaGO features

### Testing
- Added TOML-based testing system for applications
- Remove CMake-based functionality tests

### Miscallenous
- Added change log file CHANGELOG.md

## [v1.0.0] - 3-31-2021

- First public release of ExaGO

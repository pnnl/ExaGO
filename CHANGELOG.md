# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

* Added option `EXAGO_GEN_CODECOV` to generate code coverage reports.

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


### PFLOW

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


### TCOPFLOW

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
- Added file `doc/InstallingWithSpack.md` with updated spack documentation
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


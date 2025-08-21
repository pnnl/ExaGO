# ExaGO™ - Exascale Grid Optimization Toolkit
**Advanced High-Performance Computing Platform for Large-Scale Power Grid Optimization**

[![Build Status](https://github.com/pnnl/ExaGO/actions/workflows/spack_cpu_build.yaml/badge.svg)](https://github.com/pnnl/ExaGO/actions) [![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE) [![Version](https://img.shields.io/badge/version-1.6.0-green.svg)](https://github.com/pnnl/ExaGO/releases)

## 🚀 Quick Start
**Get running in 3 steps:**
1. **Install dependencies**: PETSc, CMake, solvers (Ipopt/HiOp)
2. **Build with CMake**: `mkdir build && cd build && cmake .. && make`
3. **Run application**: `./applications/opflow_main -netfile datafiles/case9/case9mod.m`

[Detailed steps below ⬇️](#installation)

---

## 📋 Project Overview

**What it does**: ExaGO™ is a cutting-edge C++/C toolkit for solving large-scale power grid optimization problems on parallel and distributed architectures, specifically designed for exascale machines with heterogeneous architectures (CPU/GPU). It solves complex combinations of stochastic, contingency-constrained, and multi-period AC optimal power flow (ACOPF) problems.

**Who it's for**: 
- **Power system researchers** working on grid optimization
- **HPC developers** building scientific computing applications
- **Operations research scientists** solving large-scale optimization problems
- **Utility companies** performing power system analysis
- **Academic researchers** in electrical engineering and applied mathematics

**Key features**:
- ⚡ **Multi-application suite** for different power flow optimization problems
- 🖥️ **CPU/GPU execution** with CUDA and HIP support  
- 🌐 **Massively parallel** with MPI and distributed computing
- 🧮 **Advanced solvers** integration (Ipopt, HiOp)
- 🐍 **Python bindings** for easy integration
- 📊 **Comprehensive benchmarking** and performance analysis tools
- 🗺️ **Interactive visualization** platform

**Tech stack**:
- **Languages**: C++14, C, Python
- **Build System**: CMake 3.18+
- **Parallel Computing**: MPI, CUDA, HIP, OpenMP
- **Core Libraries**: PETSc, RAJA, Umpire
- **Optimization Solvers**: Ipopt, HiOp (sparse/dense)
- **Python Integration**: pybind11
- **Package Management**: Spack

**Architecture style**: High-performance scientific computing library with modular solver architecture and distributed computing capabilities.

---

## 📁 Project Structure

ExaGO follows a scientific/HPC project structure optimized for large-scale computational applications:

```
/
├── 📁 src/                           # Core source implementations
│   ├── 📁 opflow/                    # Optimal power flow solver
│   ├── 📁 scopflow/                  # Security-constrained OPF
│   ├── 📁 sopflow/                   # Stochastic OPF  
│   ├── 📁 tcopflow/                  # Multi-period OPF
│   ├── 📁 pflow/                     # AC power flow
│   ├── 📁 ps/                        # Power system data structures
│   └── 📁 utils/                     # Utility functions
├── 📁 applications/                  # Executable solver programs
│   ├── 📄 opflow_main.cpp           # AC optimal power flow
│   ├── 📄 scopflow_main.cpp         # Security-constrained OPF
│   ├── 📄 sopflow_main.cpp          # Stochastic OPF
│   ├── 📄 tcopflow_main.cpp         # Multi-period OPF
│   └── 📄 pflow_main.cpp            # AC power flow
├── 📁 interfaces/                    # Language bindings
│   └── 📁 python/                    # Python API bindings
├── 📁 include/                       # Public header files
│   ├── 📄 opflow.h                   # OPFLOW API
│   ├── 📄 scopflow.h                 # SCOPFLOW API  
│   ├── 📄 sopflow.h                  # SOPFLOW API
│   ├── 📄 tcopflow.h                 # TCOPFLOW API
│   └── 📁 private/                   # Internal headers
├── 📁 datafiles/                     # Test cases and validation data
│   ├── 📁 case9/                     # 9-bus test system
│   ├── 📁 unit/                      # Unit test cases
│   └── 📄 case*.m                    # MATPOWER format data files
├── 📁 tests/                         # Comprehensive test suite
│   ├── 📁 unit/                      # Unit tests
│   ├── 📁 functionality/             # Integration tests
│   └── 📁 interfaces/                # API binding tests
├── 📁 performance_analysis/          # Benchmarking and profiling
│   ├── 📄 perf_pipeline.py          # Performance analysis pipeline
│   └── 📄 *.toml                     # Test suite configurations
├── 📁 buildsystem/                   # Multi-platform build scripts
│   ├── 📄 build.sh                   # Main build script
│   ├── 📁 spack/                     # Spack environment configs
│   ├── 📁 cmake/                     # CMake modules
│   └── 📁 [system]/                  # Platform-specific configs
├── 📁 docs/                          # Comprehensive documentation
│   ├── 📁 manual/                    # LaTeX user manual
│   ├── 📁 web/                       # Application documentation
│   └── 📁 developer_guidelines.md   # Development standards
├── 📁 options/                       # Solver configuration files
│   ├── 📄 ipopt.opt                  # Ipopt solver options
│   ├── 📄 hiop.options               # HiOp solver options
│   └── 📄 *options                   # Application-specific configs
├── 📁 viz/                           # Interactive visualization platform
│   ├── 📄 app.js                     # Web application
│   ├── 📄 index.html                 # Main interface
│   └── 📄 geninputfile.py           # Data preprocessing
├── 📁 tutorials/                     # Jupyter notebook tutorials
│   ├── 📄 demo1.ipynb               # Basic usage examples
│   └── 📄 demo2.ipynb               # Advanced features
├── 📄 CMakeLists.txt                # Main CMake configuration
├── 📄 INSTALL.md                    # Installation instructions
└── 📄 DEV-README.md                 # This developer guide
```

### **Folder Explanations:**

**Core Implementation (`src/`)**:
- **`opflow/`**: AC optimal power flow solver with multiple models (polar, cartesian, current balance)
- **`scopflow/`**: Security-constrained OPF handling contingency scenarios  
- **`sopflow/`**: Stochastic OPF for uncertainty quantification
- **`tcopflow/`**: Multi-period OPF for time-series optimization
- **`pflow/`**: Basic AC power flow solver for initialization
- **`ps/`**: Power system data structures and network representation
- **`utils/`**: Common utilities, logging, and helper functions

**Applications (`applications/`)**:
- Stand-alone executable programs for each solver type
- Command-line interfaces with comprehensive option parsing
- MPI-enabled for parallel execution on HPC systems

**Language Bindings (`interfaces/`)**:
- **`python/`**: Complete Python API using pybind11 for all applications
- Provides Pythonic interface to C++/C core functionality

**Build System (`buildsystem/`)**:
- **`build.sh`**: Universal build script with HPC system detection
- **`spack/`**: Spack package manager integration
- **`cmake/`**: Custom CMake modules and find scripts
- **Platform directories**: System-specific configurations (GPU, MPI, etc.)

---

## ⚡ Installation

### Prerequisites

**For C++/CMake/HPC Projects:**
- **CMake** 3.18+ (required for modern C++ features and GPU support)
- **C++14 compatible compiler**: GCC 9+, Clang 10+, Intel 2019+, MSVC 2019+
- **PETSc** ==3.16 (core dependency, must match version exactly)
- **MPI** 3.0+ (OpenMPI, MPICH, Intel MPI for parallel execution)
- **BLAS/LAPACK** (linear algebra backend)

**Optimization Solvers (at least one required)**:
- **Ipopt** ≥3.12 (interior-point nonlinear optimization)
- **HiOp** v0.5.3 (HPC-optimized mixed sparse-dense solver)

**Optional HPC/GPU Dependencies**:
- **CUDA** 11.0+ (for GPU acceleration with HiOp)
- **HIP/ROCm** (for AMD GPU support)  
- **RAJA** + **Umpire** (for portable GPU programming)
- **Spack** (for dependency management)

### 1. Clone the Repository
```bash
git clone https://github.com/pnnl/ExaGO.git
cd ExaGO

# Initialize submodules for third-party libraries
git submodule update --init --recursive
```

### 2. Environment Setup

**For HPC Systems with Spack:**
```bash
# Load spack environment (automatically detects system)
source buildsystem/spack/load_spack.sh

# Install all dependencies
spack install

# Activate ExaGO environment
spack env activate exago-env
```

**For Manual Dependencies:**
```bash
# Set environment variables for dependencies
export PETSC_DIR=/path/to/petsc
export IPOPT_DIR=/path/to/ipopt  
export HIOP_DIR=/path/to/hiop
export MPI_HOME=/path/to/mpi

# Load required modules on HPC systems
module load gcc/9.3.0
module load cmake/3.20.0
module load openmpi/4.1.0
```

**For CUDA/GPU Systems:**
```bash
# Additional GPU environment setup
module load cuda/11.8
export CUDA_HOME=/usr/local/cuda

# For HiOp GPU support
export HIOP_WITH_CUDA=ON
export HIOP_WITH_GPU=ON
```

### 3. Build Process

**Using ExaGO Build Script (Recommended for HPC):**
```bash
# Automatic build with system detection and all options
cd buildsystem
./build.sh --build-only

# For specific systems (override auto-detection)
MY_CLUSTER=crusher ./build.sh --build-only
MY_CLUSTER=ascent ./build.sh --build-only
```

**Manual CMake Build:**
```bash
# Create build directory
mkdir build && cd build

# Configure with CMake (basic configuration)
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_IPOPT=ON \
  -DEXAGO_ENABLE_HIOP=ON \
  -DEXAGO_ENABLE_PYTHON=ON

# Build with parallel compilation
make -j$(nproc)

# Install to system or custom location
make install
```

**Advanced CMake Configuration:**
```bash
# Full-featured build with GPU support
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_IPOPT=ON \
  -DEXAGO_ENABLE_HIOP=ON \
  -DEXAGO_ENABLE_GPU=ON \
  -DEXAGO_ENABLE_CUDA=ON \
  -DEXAGO_ENABLE_PYTHON=ON \
  -DEXAGO_ENABLE_LOGGING=ON \
  -DPETSC_DIR=${PETSC_DIR} \
  -DIPOPT_DIR=${IPOPT_DIR} \
  -DHIOP_DIR=${HIOP_DIR} \
  -DCMAKE_INSTALL_PREFIX=/path/to/install
```

### 4. Verification
```bash
# Run comprehensive test suite
cd build
ctest --output-on-failure

# Test specific components
ctest -R opflow                # OPFLOW tests only
ctest -R unit                  # Unit tests only  
ctest -R functionality         # Integration tests

# Quick verification with example
./applications/opflow_main -netfile ../datafiles/case9/case9mod.m -print_output
```

**Expected Output:**
```
[ExaGO] Creating OPFlow
[ExaGO] Reading network data from ../datafiles/case9/case9mod.m
[ExaGO] OPFLOW: Using IPOPT solver
[ExaGO] OPFLOW: Setup completed
...
[ExaGO] OPFLOW: Optimal solution found
[ExaGO] Objective function value: 5296.689
```

---

## 🔧 Configuration

### Project-Specific Configuration

**CMake Build Configuration (`CMakeLists.txt`)**:
- **Main build settings**: Compiler flags, optimization levels, feature toggles
- **Dependency discovery**: Automatic detection of PETSc, solvers, MPI
- **Component selection**: Enable/disable applications and bindings
- **Installation paths**: Configure install directories and data file locations

**Solver Option Files (`options/`)**:
- **`ipopt.opt`**: Ipopt solver parameters (convergence, iterations, verbosity)
- **`hiop.options`**: HiOp solver configuration (compute mode, memory space)
- **`*options`**: Application-specific runtime configuration files
- **Custom solver files**: User-defined parameter sets for specific problems

**Build System Configurations (`buildsystem/`)**:
- **`build.sh`**: Universal build script with system auto-detection
- **Platform directories**: HPC system-specific environment setups
- **Spack environments**: Dependency management and environment activation
- **Container configurations**: Docker and Singularity build files

### Configuration Files Explained

**Core Configuration Files:**
- **`CMakeLists.txt`**: Primary build configuration defining dependencies, features, and build targets
- **`buildsystem/build.sh`**: Intelligent build script that auto-detects HPC systems and loads appropriate environments
- **`options/ipopt.opt`**: Ipopt nonlinear solver parameters controlling convergence criteria and algorithm behavior
- **`options/hiop.options`**: HiOp solver configuration for CPU/GPU execution modes and memory management

**System-Specific Configurations:**
- **`buildsystem/crusher/`**: ORNL Crusher (AMD GPU) specific environment and build settings
- **`buildsystem/ascent/`**: ORNL Ascent (NVIDIA GPU) system configuration  
- **`buildsystem/spack/`**: Spack package manager environments for different system architectures
- **Platform detection**: Automatic hostname-based system detection and environment loading

**Runtime Configuration Options:**
- **Application options**: Command-line parameters for solver selection, model types, and problem setup
- **Output formats**: MATPOWER, JSON, CSV output configuration
- **Logging levels**: Comprehensive logging system with configurable verbosity
- **Performance profiling**: PETSc logging and profiling configuration options

---

## 📄 Key Files Explained

### Core Application Files

**Executable Applications (`applications/`)**:
- **`opflow_main.cpp`**: AC optimal power flow solver with model/solver selection and comprehensive output options
- **`scopflow_main.cpp`**: Security-constrained OPF handling contingency scenarios with parallel processing capabilities  
- **`sopflow_main.cpp`**: Stochastic OPF for uncertainty quantification in renewable generation and load
- **`tcopflow_main.cpp`**: Multi-period OPF for time-series optimization with temporal constraints
- **`pflow_main.cpp`**: Basic AC power flow solver used for initialization and feasibility analysis

**Core Library Headers (`include/`)**:
- **`opflow.h`**: Complete OPFLOW API with model definitions, solver options, and function declarations
- **`scopflow.h`**: SCOPFLOW interface for contingency-constrained optimization problems
- **`sopflow.h`**: SOPFLOW API for stochastic optimization with scenario-based uncertainty modeling
- **`tcopflow.h`**: TCOPFLOW interface for multi-period optimization problems
- **`ps.h`**: Power system data structure definitions and network topology management

**Implementation Core (`src/`)**:
- **`src/opflow/interface/opflow.cpp`**: Main OPFLOW implementation with solver integration and problem setup
- **`src/scopflow/interface/scopflow.cpp`**: SCOPFLOW core logic for contingency scenario management
- **`src/ps/`**: Power system data structures, network parsing, and constraint formulation
- **`src/utils/`**: Logging, error handling, mathematical utilities, and helper functions

### Build & Development Files

**Build Configuration:**
- **`CMakeLists.txt`**: Master build configuration with dependency discovery, feature flags, and target definitions
- **`buildsystem/build.sh`**: Intelligent build script with HPC system auto-detection and parallel build management
- **`buildsystem/cmake/`**: Custom CMake modules for finding dependencies and configuring complex builds
- **`INSTALL.md`**: Comprehensive installation guide with dependency requirements and build instructions

**Development & Testing:**
- **`tests/CMakeLists.txt`**: Test suite configuration with unit, functionality, and integration test definitions
- **`tests/functionality/`**: End-to-end application testing with real power system test cases
- **`tests/unit/`**: Component-level testing for individual functions and classes
- **`docs/developer_guidelines.md`**: Coding standards, contribution guidelines, and development best practices

**Performance & Analysis:**
- **`performance_analysis/perf_pipeline.py`**: Automated performance benchmarking pipeline with TOML configuration
- **`performance_analysis/*.toml`**: Test suite configurations for scalability analysis and solver comparison
- **Profiling integration**: PETSc logging stages and performance event tracking

---

## 💻 Development Workflow

### Daily Development (HPC/Scientific Computing Focus)

**For ExaGO C++/HPC Development:**
1. **Pull latest changes**: `git pull origin develop`
2. **Update submodules**: `git submodule update --init --recursive`  
3. **Load environment**: `source buildsystem/spack/load_spack.sh` or `module load [system-modules]`
4. **Rebuild if needed**: `cd build && make -j$(nproc)` or `cd buildsystem && ./build.sh`
5. **Run tests**: `cd build && ctest` or `ctest -R [component]`

**For Python Interface Development:**
1. **Activate environment**: `spack env activate exago-env`
2. **Build Python bindings**: `make -j$(nproc) && make install`
3. **Test Python interface**: `python -c "import exago; print('Success')"`
4. **Run Python examples**: `python interfaces/python/example_opflow.py`

**For Performance/HPC Development:**
1. **Load HPC modules**: `module load [compiler] [mpi] [cuda]`
2. **Build with profiling**: `cmake .. -DEXAGO_ENABLE_LOGGING=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo`  
3. **Run performance tests**: `cd performance_analysis && python perf_pipeline.py`
4. **Analyze scaling**: Review generated performance reports and scaling plots

### Making Changes

1. **Create feature branch**: `git checkout -b feature/your-feature`
2. **Make your changes** following ExaGO coding standards:
   - Use clang-format for C++ code formatting
   - Follow PETSc error handling conventions  
   - Add comprehensive documentation for new functions
   - Include unit tests for new functionality
3. **Build and test**: 
   ```bash
   cd build
   make -j$(nproc)
   ctest --output-on-failure
   ```
4. **Test on target platforms**: Verify changes work on relevant HPC systems
5. **Run performance benchmarks**: Ensure no performance regressions
6. **Commit and push**: `git add . && git commit -m "Descriptive message" && git push`
7. **Create Pull Request** with detailed description and test results

### Testing Framework

**Comprehensive Test Suite:**
```bash
# Full test suite (all components and platforms)
cd build && ctest

# Component-specific testing
ctest -R opflow                # OPFLOW solver tests
ctest -R scopflow              # SCOPFLOW tests  
ctest -R python                # Python binding tests
ctest -R unit                  # Unit tests only
ctest -R functionality         # Integration tests

# Performance testing
cd performance_analysis
python perf_pipeline.py opflow_testsuite.toml

# Parallel testing (MPI)
mpiexec -n 4 ctest -R parallel

# GPU testing (if available)
ctest -R gpu
```

**Test Categories:**
- **Unit tests**: Individual function and class testing with isolated test cases
- **Functionality tests**: End-to-end application testing with real power system problems  
- **Interface tests**: Python binding validation and API consistency checks
- **Performance tests**: Scaling analysis and solver comparison benchmarks
- **Regression tests**: Ensuring solution accuracy and convergence behavior

### Code Quality Tools

**C++ Code Quality:**
- **clang-format**: Automatic code formatting following ExaGO style guidelines
- **clang-tidy**: Static analysis for code quality and potential issues
- **PETSc conventions**: Error handling and memory management following PETSc patterns
- **Documentation**: Doxygen-style comments for all public interfaces

**Python Code Quality:**
- **Black**: Python code formatting for binding implementations
- **Type hints**: Complete type annotations for Python API
- **Docstrings**: NumPy-style documentation for all Python functions
- **Testing**: pytest-based testing framework for Python components

**Build and CI:**
- **GitHub Actions**: Automated testing on multiple platforms and compilers
- **Spack CI**: Package manager integration testing
- **Pre-commit hooks**: Automatic formatting and linting before commits
- **Performance monitoring**: Automated performance regression detection

---

## 🔌 API/Interface Reference

### Command Line Interface

**Basic Application Execution:**
```bash
# AC Optimal Power Flow
./applications/opflow_main -netfile datafiles/case9/case9mod.m -print_output

# Security-Constrained OPF  
mpiexec -n 4 ./applications/scopflow_main \
  -netfile datafiles/case9/case9mod.m \
  -ctgcfile datafiles/case9/case9.cont \
  -scopflow_Nc 4 -print_output

# Stochastic OPF
./applications/sopflow_main \
  -netfile datafiles/case9/case9mod.m \
  -windgen datafiles/case9/10_scenarios_9bus.csv \
  -sopflow_Ns 4

# Multi-period OPF
./applications/tcopflow_main \
  -netfile datafiles/case9/case9mod.m \
  -tcopflow_dT 1.0 -tcopflow_duration 24.0
```

**Advanced Options and Solver Configuration:**
```bash
# Solver and model selection
./applications/opflow_main \
  -netfile input.m \
  -opflow_solver HIOP \
  -opflow_model POWER_BALANCE_HIOP \
  -hiop_compute_mode gpu

# Output format options
./applications/opflow_main \
  -netfile input.m \
  -save_output solution \
  -opflow_output_format JSON \
  -print_output

# Optimization parameters
./applications/opflow_main \
  -netfile input.m \
  -opflow_tolerance 1e-8 \
  -opflow_initialization ACPF \
  -opflow_ignore_lineflow_constraints
```

### C++ Library Interface

**Basic OPFLOW Usage:**
```cpp
#include <opflow.h>

int main(int argc, char **argv) {
    // Initialize ExaGO and MPI
    ExaGOInitialize(MPI_COMM_WORLD, &argc, &argv, "app", help);
    
    // Create OPFLOW object
    OPFLOW opflow;
    OPFLOWCreate(MPI_COMM_WORLD, &opflow);
    
    // Read network data
    OPFLOWReadMatPowerData(opflow, "datafiles/case9/case9mod.m");
    
    // Configure solver and model
    OPFLOWSetSolver(opflow, "IPOPT");
    OPFLOWSetModel(opflow, "POWER_BALANCE_POLAR");
    
    // Setup and solve
    OPFLOWSetUp(opflow);
    OPFLOWSolve(opflow);
    
    // Get results
    PetscScalar objective;
    OPFLOWGetObjective(opflow, &objective);
    
    // Cleanup
    OPFLOWDestroy(&opflow);
    ExaGOFinalize();
    return 0;
}
```

**Advanced SCOPFLOW Configuration:**
```cpp
#include <scopflow.h>

// Create security-constrained OPF
SCOPFLOW scopflow;
SCOPFLOWCreate(MPI_COMM_WORLD, &scopflow);

// Set input files
SCOPFLOWSetNetworkData(scopflow, "network.m");
SCOPFLOWSetContingencyData(scopflow, "contingencies.cont", NATIVE);

// Configure parallel solver
SCOPFLOWSetSolver(scopflow, "HIOP");
SCOPFLOWSetSubproblemSolver(scopflow, "IPOPT");
SCOPFLOWSetNumContingencies(scopflow, 10);

// Solve and get results
SCOPFLOWSetUp(scopflow);
SCOPFLOWSolve(scopflow);

PetscScalar total_obj, base_obj;
SCOPFLOWGetTotalObjective(scopflow, &total_obj);
SCOPFLOWGetBaseObjective(scopflow, &base_obj);
```

### Python API Interface

**Basic Python Usage:**
```python
import exago
import os

# Initialize ExaGO
exago.initialize("opflow_app")

# Create optimal power flow object
opf = exago.OPFLOW()

# Load network data
path = exago.prefix()
datafile = os.path.join(path, 'share', 'exago', 'datafiles', 'case9', 'case9mod.m')
opf.read_mat_power_data(datafile)

# Configure and solve
opf.set_solver("IPOPT")
opf.set_model("POWER_BALANCE_POLAR")
opf.solve()

# Get results
objective = opf.get_objective()
converged = opf.get_convergence_status()
iterations = opf.get_num_iterations()

print(f"Objective: {objective}, Converged: {converged}, Iterations: {iterations}")

# Save solution
opf.save_solution(exago.OutputFormat.MATPOWER, "solution.m")

# Cleanup
del opf
exago.finalize()
```

**Advanced Security-Constrained Python Example:**
```python
import exago
from mpi4py import MPI

# Initialize with MPI
comm = MPI.COMM_WORLD
exago.initialize("scopflow_app", comm)

# Create SCOPFLOW object
scopf = exago.SCOPFLOW()

# Configure problem
scopf.set_network_data("network.m")
scopf.set_contingency_data("contingencies.cont", exago.ContingencyFileInputFormat.NATIVE)
scopf.set_solver("HIOP")
scopf.set_subproblem_solver("IPOPT")
scopf.set_num_contingencies(4)

# Advanced options
scopf.set_subproblem_model("POWER_BALANCE_POLAR")
scopf.set_mode(1)  # Corrective mode
scopf.set_tolerance(1e-6)

# Solve
scopf.solve()

# Analysis results
total_obj = scopf.get_total_objective()
base_obj = scopf.get_base_objective()
converged = scopf.get_convergence_status()

print(f"Total Objective: {total_obj}")
print(f"Base Objective: {base_obj}")
print(f"Converged: {converged}")

# Save all solutions
scopf.save_solution_all_default("results/")
```

### Available Applications and Models

**OPFLOW (AC Optimal Power Flow):**
- **Models**: `POWER_BALANCE_POLAR`, `POWER_BALANCE_CARTESIAN`, `POWER_BALANCE_HIOP`
- **Solvers**: `IPOPT`, `HIOP`, `HIOPSPARSE`
- **Usage**: Single-period AC OPF with various formulations and solvers

**SCOPFLOW (Security-Constrained OPF):**
- **Models**: `GENRAMP`, `GENRAMPT` (multi-period)
- **Solvers**: `IPOPT`, `HIOP`, `EMPAR` (embarrassingly parallel)
- **Usage**: Contingency-constrained optimization for grid security analysis

**SOPFLOW (Stochastic OPF):**
- **Models**: `GENRAMP`, scenario-based uncertainty modeling
- **Solvers**: `IPOPT`, `HIOP` with stochastic programming
- **Usage**: Renewable integration and uncertainty quantification

**TCOPFLOW (Multi-period OPF):**
- **Models**: Time-coupled optimization with ramping constraints
- **Solvers**: `IPOPT` for temporal optimization problems
- **Usage**: Day-ahead scheduling and multi-period planning

### Example Usage Patterns

**Performance Benchmarking:**
```bash
# Run automated performance analysis
cd performance_analysis
python perf_pipeline.py opflow_testsuite.toml

# Custom solver comparison
./applications/opflow_main -netfile large_case.m -opflow_solver IPOPT -log_view
./applications/opflow_main -netfile large_case.m -opflow_solver HIOP -log_view
```

**GPU Acceleration:**
```bash
# HiOp GPU execution
./applications/opflow_main \
  -netfile large_case.m \
  -opflow_solver HIOP \
  -opflow_model PBPOLRAJAHIOP \
  -hiop_compute_mode gpu \
  -hiop_mem_space device
```

**Large-Scale Parallel Execution:**
```bash
# Massive parallel SCOPFLOW
mpiexec -n 128 ./applications/scopflow_main \
  -netfile huge_case.m \
  -ctgcfile contingencies.cont \
  -scopflow_solver HIOP \
  -scopflow_Nc 100
```

---

## 🌍 Deployment

### HPC Production Deployment

**For Large-Scale HPC Systems:**
```bash
# Build optimized version for production
cd buildsystem
./build.sh --build-only

# Or with specific HPC system optimization
MY_CLUSTER=crusher ./build.sh --build-only
MY_CLUSTER=ascent ./build.sh --build-only
MY_CLUSTER=frontier ./build.sh --build-only

# Install to shared filesystem
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/shared/exago \
  -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_GPU=ON

make -j$(nproc) install
```

**Optimized Production Build:**
```bash
# Maximum optimization for performance
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-O3 -march=native -mtune=native" \
  -DEXAGO_ENABLE_MPI=ON \
  -DEXAGO_ENABLE_GPU=ON \
  -DEXAGO_ENABLE_HIOP=ON \
  -DEXAGO_ENABLE_CUDA=ON \
  -DHIOP_WITH_GPU=ON

make -j$(nproc)
make install
```

### Spack-Based Deployment

**Production Spack Installation:**
```bash
# Install ExaGO with all variants
spack install exago +mpi +cuda +hiop +python

# Create environment for users
spack env create exago-prod
spack env activate exago-prod
spack add exago +mpi +cuda +hiop +python
spack install

# Module generation for user access
spack module tcl refresh exago
```

**Custom Spack Configuration:**
```yaml
# spack.yaml for ExaGO environment
spack:
  specs:
  - exago +mpi +cuda +hiop +python ^petsc@3.16 +mpi +cuda
  - hiop +cuda +mpi ^cuda@11.8
  - ipopt +mumps
  concretizer:
    unify: true
```

### Container Deployment

**Docker Container Build:**
```bash
# Build ExaGO container
docker build -t exago:latest -f buildsystem/container/Dockerfile .

# Run container with GPU support
docker run --gpus all -it exago:latest

# Run with data volume mounted
docker run -v /data:/data -it exago:latest \
  ./applications/opflow_main -netfile /data/network.m
```

**Singularity for HPC:**
```bash
# Build Singularity container
singularity build exago.sif buildsystem/container/Singularity

# Run on HPC cluster with MPI
mpiexec -n 8 singularity exec exago.sif \
  ./applications/scopflow_main -netfile network.m
```

### Job Scheduler Integration

**SLURM Job Script:**
```bash
#!/bin/bash
#SBATCH --job-name=exago-opf
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --gpu=4
#SBATCH --time=01:00:00
#SBATCH --partition=gpu

# Load modules
module load gcc/9.3.0
module load openmpi/4.1.0
module load cuda/11.8
module load exago/1.6.0

# Set environment
export OMP_NUM_THREADS=2
export CUDA_VISIBLE_DEVICES=0,1,2,3

# Run ExaGO application
mpiexec -n 32 ./applications/scopflow_main \
  -netfile $CASE_FILE \
  -ctgcfile $CTGC_FILE \
  -scopflow_solver HIOP \
  -hiop_compute_mode gpu \
  -scopflow_Nc 100
```

**PBS/Torque Script:**
```bash
#!/bin/bash
#PBS -N exago-job
#PBS -l nodes=8:ppn=16:gpus=2
#PBS -l walltime=02:00:00
#PBS -q gpu

cd $PBS_O_WORKDIR

# Load ExaGO environment
source /shared/exago/env.sh

# Run large-scale stochastic OPF
mpiexec -n 128 ./applications/sopflow_main \
  -netfile large_network.m \
  -windgen scenarios.csv \
  -sopflow_Ns 1000 \
  -sopflow_solver HIOP
```

**Module File Creation:**
```tcl
#%Module1.0
##
## ExaGO 1.6.0
##
proc ModulesHelp { } {
    puts stderr "ExaGO - Exascale Grid Optimization Toolkit"
}

module-whatis "ExaGO power grid optimization toolkit"

set basedir /shared/exago/1.6.0

prepend-path PATH $basedir/bin
prepend-path LD_LIBRARY_PATH $basedir/lib
prepend-path PYTHONPATH $basedir/lib/python3.8/site-packages

setenv EXAGO_DIR $basedir
setenv EXAGO_DATA_DIR $basedir/share/exago/datafiles
```

---

## 🚨 Troubleshooting

### Common Issues (HPC/Scientific Computing)

**For C++/CMake Build Issues:**
```bash
# CMake configuration errors
cd build && rm -rf * && cmake ..

# Missing dependencies - check system-specific instructions
ls buildsystem/
cat buildsystem/README.md

# Check PETSc installation
echo $PETSC_DIR
ls $PETSC_DIR/lib

# Verbose compilation to see detailed errors
make VERBOSE=1

# Library linking issues (check dependencies)
ldd ./applications/opflow_main
```

**For MPI/Parallel Execution Issues:**
```bash
# MPI environment problems
mpirun --version
which mpiexec
echo $MPI_HOME

# Check MPI compilation
mpicc --version
mpicxx --version

# Test basic MPI functionality  
mpiexec -n 2 hostname

# Debug MPI issues with ExaGO
mpiexec -n 2 ./applications/opflow_main -netfile datafiles/case9/case9mod.m -print_output
```

**For GPU/CUDA Issues:**
```bash
# Check GPU availability
nvidia-smi
nvcc --version

# CUDA environment verification
echo $CUDA_HOME
ls $CUDA_HOME/lib64

# Test GPU memory and compute capability
deviceQuery  # CUDA samples utility

# HiOp GPU solver debugging
./applications/opflow_main \
  -netfile small_case.m \
  -opflow_solver HIOP \
  -hiop_compute_mode gpu \
  -hiop_verbosity_level 5
```

**For Solver-Related Issues:**
```bash
# Ipopt solver problems
./applications/opflow_main -netfile case.m -opflow_solver IPOPT -opflow_tolerance 1e-4

# HiOp solver issues
./applications/opflow_main -netfile case.m -opflow_solver HIOP -hiop_verbosity_level 3

# Solver convergence problems
./applications/opflow_main \
  -netfile case.m \
  -opflow_initialization ACPF \
  -opflow_tolerance 1e-6 \
  -print_output

# Check solver installation
ls $IPOPT_DIR/lib
ls $HIOP_DIR/lib
```

### Build/Compilation Issues

**Dependency Resolution:**
```bash
# Verify all required dependencies
cd buildsystem
./build.sh --check-deps

# PETSc version mismatch
pkg-config --modversion petsc  # Should be 3.16.x

# Spack dependency issues
spack find exago
spack spec exago +mpi +cuda +hiop

# Clean rebuild process
rm -rf build install
mkdir build && cd build
cmake .. -DCMAKE_VERBOSE_MAKEFILE=ON
```

**Performance Issues:**
```bash
# Memory usage profiling
valgrind --tool=memcheck ./applications/opflow_main -netfile large_case.m

# Performance profiling with PETSc
./applications/opflow_main -netfile case.m -log_view -log_summary

# MPI scaling analysis
for np in 1 2 4 8; do
  echo "Testing with $np processes"
  mpiexec -n $np ./applications/scopflow_main -netfile case.m -log_view
done

# GPU utilization monitoring
nvidia-smi -l 1  # Run while executing GPU jobs
```

**Python Interface Issues:**
```bash
# Python import errors
python -c "import sys; print('\n'.join(sys.path))"
python -c "import exago; print('Success')"

# PYTHONPATH configuration
export PYTHONPATH=$EXAGO_INSTALL_DIR/lib/python3.8/site-packages:$PYTHONPATH

# Rebuild Python bindings
cd build
make -j$(nproc)
make install

# Test Python examples
cd interfaces/python
python example_opflow.py
```

### Solver-Specific Debugging

**Ipopt Debugging:**
```bash
# Enable Ipopt detailed output
./applications/opflow_main \
  -netfile problematic_case.m \
  -opflow_solver IPOPT \
  -options_file ipopt_debug.opt

# ipopt_debug.opt contents:
# print_level 5
# output_file ipopt_debug.out
# derivative_test first-order
```

**HiOp Debugging:**
```bash
# HiOp verbose debugging
./applications/opflow_main \
  -netfile case.m \
  -opflow_solver HIOP \
  -hiop_verbosity_level 10 \
  -hiop_compute_mode cpu

# GPU memory issues
export CUDA_VISIBLE_DEVICES=0
./applications/opflow_main \
  -netfile case.m \
  -opflow_solver HIOP \
  -hiop_compute_mode gpu \
  -hiop_mem_space device
```

**Large-Scale Problem Debugging:**
```bash
# Memory requirements estimation
./applications/opflow_main -netfile huge_case.m -print_problem_size

# Parallel load balancing
mpiexec -n 16 ./applications/scopflow_main \
  -netfile large_case.m \
  -scopflow_solver HIOP \
  -log_view -log_summary

# I/O and data loading issues
./applications/opflow_main -netfile case.m -opflow_output_format MINIMAL
```

### Debug Mode and Logging

**Enable Debug Builds:**
```bash
# Debug build with symbols
cmake .. -DCMAKE_BUILD_TYPE=Debug -DEXAGO_ENABLE_LOGGING=ON
make -j$(nproc)

# Run with GDB
gdb --args ./applications/opflow_main -netfile case.m
(gdb) run
(gdb) bt  # If crash occurs
```

**Advanced Logging:**
```bash
# Set log file output
export EXAGO_LOG_FILE=exago_debug.log

# Enable detailed PETSc logging  
./applications/opflow_main \
  -netfile case.m \
  -log_view -log_summary \
  -info -options_view

# Performance event logging
./applications/scopflow_main \
  -netfile case.m \
  -log_trace performance.trace \
  -log_view
```

---

## 🤝 Contributing

### Development Setup
1. **Fork the repository** on GitHub
2. **Clone your fork**: `git clone https://github.com/yourusername/ExaGO.git`
3. **Set up development environment**: 
   ```bash
   cd ExaGO
   source buildsystem/spack/load_spack.sh  # Or load manual dependencies
   ```
4. **Create feature branch**: `git checkout -b feature/your-feature-name`
5. **Make your changes** following coding standards below
6. **Add comprehensive tests** for new functionality
7. **Ensure all tests pass**: `cd build && ctest`
8. **Submit pull request** with detailed description

### Code Standards

**C++ Coding Standards:**
- **Follow PETSc conventions**: Use PETSc error handling (`CHKERRQ`, `PetscFunctionBegin/Return`)
- **Use clang-format**: Automatic formatting with project `.clang-format` configuration
- **Memory management**: Follow PETSc object lifecycle patterns, no raw pointers
- **Documentation**: Doxygen-style comments for all public functions and classes
- **Naming conventions**: PetscCase for functions, lowercase_underscore for variables

**Python Coding Standards:**
- **Follow PEP 8**: Use Black formatter for automatic compliance
- **Type hints**: Complete type annotations for all Python bindings
- **Docstrings**: NumPy-style documentation with examples
- **Error handling**: Proper exception handling and informative error messages

**General Standards:**
- **Commit messages**: Use conventional commit format (`feat:`, `fix:`, `docs:`, etc.)
- **Documentation**: Update relevant documentation for any user-facing changes
- **Performance**: Benchmark new features to ensure no regressions
- **Testing**: Comprehensive unit and integration tests required

### Pull Request Process

1. **Update documentation** if your changes affect user-facing functionality
2. **Add tests** for new features and ensure existing tests pass
3. **Run performance benchmarks** to verify no significant regressions
4. **Ensure CI passes** on all supported platforms and compilers
5. **Request review** from relevant maintainers based on changed components:
   - Core solvers: @core-team
   - Python bindings: @python-team  
   - HPC/performance: @hpc-team
   - Documentation: @docs-team

**Review Checklist:**
- [ ] Code follows ExaGO style guidelines
- [ ] All tests pass on relevant platforms
- [ ] Documentation updated for user-facing changes
- [ ] Performance impact assessed and acceptable
- [ ] Security considerations addressed for new features
- [ ] Backward compatibility maintained unless explicitly breaking

### Testing Requirements

**Required Testing:**
- **Unit tests**: For all new functions and classes
- **Integration tests**: End-to-end testing with real power system cases  
- **Python binding tests**: API consistency and error handling
- **Performance tests**: Scaling and regression analysis for performance-critical changes
- **Platform testing**: Verification on major HPC systems when applicable

**Test Execution:**
```bash
# Run all tests
cd build && ctest

# Component-specific tests
ctest -R component_name

# Performance regression testing
cd performance_analysis
python perf_pipeline.py --compare baseline_results/
```

---

## 📖 Additional Resources

### Documentation

**Core Documentation:**
- **[ExaGO Manual](docs/manual/manual.pdf)**: Comprehensive mathematical formulations and theoretical background
- **[OPFLOW Documentation](docs/web/opflow.md)**: Complete guide to AC optimal power flow application
- **[SCOPFLOW Documentation](docs/web/scopflow.md)**: Security-constrained OPF user guide and examples
- **[Installation Guide](INSTALL.md)**: Detailed dependency and build instructions
- **[Developer Guidelines](docs/developer_guidelines.md)**: Coding standards and contribution process

**API References:**
- **[Python API Examples](interfaces/python/)**: Complete Python binding usage examples
- **[C++ API Headers](include/)**: Full C++/C interface documentation
- **[Jupyter Tutorials](tutorials/)**: Interactive notebooks with step-by-step examples
- **[Application Options](docs/web/)**: Comprehensive command-line reference

### External Dependencies Documentation

**Core Libraries:**
- **[PETSc Documentation](https://petsc.org/release/)**: Distributed linear algebra and solvers
- **[MPI Documentation](https://www.mpi-forum.org/)**: Message Passing Interface standard
- **[CMake Documentation](https://cmake.org/documentation/)**: Build system configuration

**Optimization Solvers:**
- **[Ipopt Documentation](https://coin-or.github.io/Ipopt/)**: Interior-point nonlinear optimization
- **[HiOp Documentation](https://github.com/LLNL/hiop)**: HPC-optimized optimization solver
- **[RAJA Documentation](https://raja.readthedocs.io/)**: Portable GPU programming model
- **[Spack Documentation](https://spack.readthedocs.io/)**: Package manager for HPC

**HPC and GPU Programming:**
- **[CUDA Toolkit Documentation](https://docs.nvidia.com/cuda/)**: NVIDIA GPU programming
- **[HIP Documentation](https://rocmdocs.amd.com/en/latest/Programming_Guides/HIP-GUIDE.html)**: AMD GPU programming
- **[OpenMPI Documentation](https://www.open-mpi.org/doc/)**: MPI implementation

### Community and Support

**Getting Help:**
- **[GitHub Issues](https://github.com/pnnl/ExaGO/issues)**: Bug reports and feature requests
- **[GitHub Discussions](https://github.com/pnnl/ExaGO/discussions)**: Community Q&A and general discussions
- **[Developer Mailing List](mailto:exago-dev@pnnl.gov)**: Development coordination and technical discussions

**Contributing:**
- **[Contributing Guide](docs/developer_guidelines.md)**: How to contribute code and documentation
- **[Code of Conduct](CODE_OF_CONDUCT.md)**: Community standards and expectations
- **[Security Policy](SECURITY.md)**: How to report security vulnerabilities

**Project Information:**
- **[ExaSGD Project](https://www.exascaleproject.org/research-project/exasgd/)**: Broader exascale computing initiative
- **[PNNL Research](https://www.pnnl.gov/projects/exascale-grid-optimization-exago)**: Institutional project page
- **[Publications](docs/publications.md)**: Academic papers and conference presentations

### Related Tools and Ecosystems

**Power System Analysis:**
- **[MATPOWER](https://matpower.org/)**: MATLAB-based power system simulation
- **[PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl)**: Julia-based power network optimization
- **[PYPOWER](https://github.com/rwl/PYPOWER)**: Python port of MATPOWER

**HPC Scientific Computing:**
- **[Trilinos](https://trilinos.github.io/)**: Collection of HPC solver libraries
- **[SUNDIALS](https://computing.llnl.gov/projects/sundials)**: Suite of differential equation solvers
- **[TAO](https://www.mcs.anl.gov/petsc/petsc-tao/)**: Toolkit for Advanced Optimization

**Visualization and Analysis:**
- **[ParaView](https://www.paraview.org/)**: Scientific data visualization
- **[VisIt](https://wci.llnl.gov/simulation/computer-codes/visit)**: Interactive visualization tool
- **[Jupyter](https://jupyter.org/)**: Interactive computing environment

---

## 🏗️ Architecture and Design

### System Architecture

ExaGO follows a **modular, hierarchical architecture** designed for HPC environments:

```
┌─────────────────────────────────────────────────────────────┐
│                    APPLICATION LAYER                        │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌────────┐ │
│  │   OPFLOW    │ │  SCOPFLOW   │ │   SOPFLOW   │ │TCOPFLOW│ │
│  └─────────────┘ └─────────────┘ └─────────────┘ └────────┘ │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│                      SOLVER LAYER                           │
│  ┌─────────────────┐           ┌─────────────────────────┐   │
│  │      IPOPT      │           │         HIOP            │   │
│  │  Interior-Point │           │ ┌─────────┐ ┌─────────┐ │   │
│  │   Nonlinear     │           │ │ Dense   │ │ Sparse  │ │   │
│  │  Optimization   │           │ │CPU/GPU  │ │   CPU   │ │   │
│  └─────────────────┘           │ └─────────┘ └─────────┘ │   │
└─────────────────────────────────┘───────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│                     CORE LIBRARY                            │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────────────────┐ │
│  │ Power System│ │   Utilities │ │    Data Structures      │ │
│  │   (PS)      │ │   Logging   │ │   Network Topology      │ │
│  │ Data Model  │ │   Math      │ │   Constraint Matrices   │ │
│  └─────────────┘ └─────────────┘ └─────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│                  FOUNDATION LAYER                           │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────────────────┐ │
│  │    PETSc    │ │     MPI     │ │       GPU Runtime       │ │
│  │Linear Algebra│ │Parallel Comm│ │    CUDA / HIP / OpenMP  │ │
│  │   Solvers   │ │             │ │      RAJA / Umpire      │ │
│  └─────────────┘ └─────────────┘ └─────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

### Key Design Principles

1. **Modularity**: Each application (OPFLOW, SCOPFLOW, etc.) is independent but shares common infrastructure
2. **Scalability**: MPI-based parallelization with support for thousands of cores
3. **Portability**: Runs on diverse HPC architectures (CPU, GPU, heterogeneous systems)
4. **Extensibility**: Plugin architecture for new solvers and models
5. **Performance**: Zero-copy interfaces and memory-efficient data structures

### Data Flow and Processing

**Typical ExaGO Workflow:**
1. **Data Input**: Parse MATPOWER files or custom network formats
2. **Problem Setup**: Create optimization variables and constraint matrices
3. **Solver Interface**: Interface with Ipopt/HiOp through standardized API
4. **Parallel Execution**: Distribute computations across MPI processes/GPU threads
5. **Solution Processing**: Extract results and generate output in various formats
6. **Visualization**: Optional web-based visualization of results

This architecture enables ExaGO to handle power networks with millions of buses and thousands of contingencies while maintaining numerical accuracy and computational efficiency.

---

*This DEV-README.md provides comprehensive developer documentation for the ExaGO project. For user documentation, see the [main README.md](README.md) and [ExaGO manual](docs/manual/manual.pdf).*

**ExaGO Development Team**  
Pacific Northwest National Laboratory  
Battelle Memorial Institute

---

**Copyright © 2020-2025, Battelle Memorial Institute**  
**ExaGO™ is developed as part of the ExaSGD project under the Exascale Computing Project.**

#[[#############################################################################

Preamble

#]]

cmake_minimum_required(VERSION 3.18)

project(ExaGO VERSION "1.2.0")

string(TIMESTAMP EXAGO_RELEASE_DATE "%Y-%m-%d")

set(CMAKE_CXX_STANDARD 11)
# Enable setting Mac OS X rpath
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # CMake 3.12
endif()

if(POLICY CMP0104)
  cmake_policy(SET CMP0104 NEW) # CMake 3.18
endif()

# Make CMake debugging available globally
include(CMakePrintHelpers)

#[[#############################################################################

Configuration and library discovery.

#]]

# Specify where CMake modules are to be found
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/buildsystem/cmake)
include(ExaGOCheckGitSubmodules)

# Set the path to options file in EXAGO installation
set(EXAGO_OPTIONS_DIR "${CMAKE_INSTALL_PREFIX}/share/exago/options")
set(EXAGO_DATAFILES_DIR "${CMAKE_INSTALL_PREFIX}/share/exago/datafiles")

# Specify whether to build shared libraries, static libraries, or both
option(EXAGO_BUILD_SHARED "Build shared libraries" OFF)
option(EXAGO_BUILD_STATIC "Build static libraries" ON)
include(ExaGOAddLibrary)

# Enable PETSc by default
option(EXAGO_ENABLE_PETSC "Enable PETSc support" ON)

# MPI support is enabled by default
option(EXAGO_ENABLE_MPI "Enable MPI support" ON)

# Ipopt support is disabled by default
option(EXAGO_ENABLE_IPOPT "Enable Ipopt support" OFF)

# HiOP support is disabled by default
option(EXAGO_ENABLE_HIOP "Enable HiOP support" OFF)

# MPI support is enabled by default
option(EXAGO_ENABLE_MPI "Enable MPI support" ON)

# GPU support is disabled by default
option(EXAGO_ENABLE_GPU "Enable GPU support" OFF)

include(CMakeDependentOption)

# CUDA is disabled by default, but is default when EXAGO_ENABLE_GPU=ON
cmake_dependent_option(
  EXAGO_ENABLE_CUDA "Build with support for CUDA" OFF "EXAGO_ENABLE_GPU" OFF
)

# HIP is disabled by default
cmake_dependent_option(
  EXAGO_ENABLE_HIP "Build with support for HIP" OFF "EXAGO_ENABLE_GPU" OFF
)

# Enable RAJA support
option(EXAGO_ENABLE_RAJA "Enables RAJA and Umpire libraries" OFF)

# Enable Python interface
option(EXAGO_ENABLE_PYTHON "Enables Python bindings and tests" ON)

# Enable OpenMP support. This should not be toggled directly by the user, which
# is why it's marked as advanced.
option(EXAGO_ENABLE_OMP "Enables OpenMP" ON)
mark_as_advanced(FORCE EXAGO_ENABLE_OMP)

# Enable tests
option(EXAGO_RUN_TESTS "Enable tests" ON)

# If tests are enabled, decides if tests are installed
option(EXAGO_INSTALL_TESTS "Install testing binaries" ON)

# For testing with BSUB/jsrun on certain clusters
option(EXAGO_TEST_WITH_BSUB
       "Use `jsrun` instead of `mpirun` commands when running tests" OFF
)

# For testing with Slurm/srun on certain clusters
option(EXAGO_TEST_WITH_SLURM
       "Use `srun` instead of `mpirun` commands when running tests" OFF
)

option(EXAGO_DISABLE_LOGGING "Disable all internal logging" OFF)

# When building with HIP support, OpenMP is not supported
if(NOT EXAGO_ENABLE_RAJA OR (EXAGO_ENABLE_GPU AND EXAGO_ENABLE_HIP))
  set(EXAGO_ENABLE_OMP OFF)
endif()

# Estimate number of flops per trig op
set(EXAGO_FLOPS_SINOP
    15
    CACHE STRING "Estimated number of flops per sin op"
)
set(EXAGO_FLOPS_COSOP
    15
    CACHE STRING "Estimated number of flops per cos op"
)
set(EXAGO_FLOPS_TANOP
    18
    CACHE STRING "Estimated number of flops per tan op"
)
set(EXAGO_FLOPS_SQRTOP
    6
    CACHE STRING "Estimated number of flops per sqrt op"
)

# Set default include path for the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})

# Create alias for -lm which is needed when linking math on UNIX
if(UNIX)
  set(EXAGO_MATH_LIB m)
endif()

# If MPI support is enabled, find MPI includes and libraries
if(EXAGO_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C CXX)
  if(NOT DEFINED MPI_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  endif(NOT DEFINED MPI_CXX_COMPILER)
  include_directories(
    ${MPI_CXX_HEADER_DIR} ${MPI_CXX_ADDITIONAL_INCLUDE_DIRS}
    ${MPI_CXX_COMPILER_INCLUDE_DIRS}
  )
  set(EXAGO_EXTRA_MPI_FLAGS
      ""
      CACHE STRING "Extra args to mpiexec when running tests"
  )
endif(EXAGO_ENABLE_MPI)

# If GPU support is enabled, find CUDA
if(EXAGO_ENABLE_GPU)

  if(NOT EXAGO_ENABLE_CUDA AND NOT EXAGO_ENABLE_HIP)
    message(
      WARNING
        "\nGPU is enabled, but both CUDA and HIP are disabled. "
        "Enabling CUDA by default. If you would like to enable a particular GPU "
        "platform, set one of the following variables to ON: "
        "\n* EXAGO_ENABLE_CUDA"
        "\n* EXAGO_ENABLE_HIP"
    )
    set(EXAGO_ENABLE_CUDA ON)
  endif()

  if(EXAGO_ENABLE_CUDA)
    include(CheckLanguage)
    enable_language(CUDA)
    check_language(CUDA)

    find_package(CUDA REQUIRED)

    if(NOT DEFINED CMAKE_CUDA_STANDARD)
      set(CMAKE_CUDA_STANDARD 11)
      set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES 35)
    endif()

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")

    # HiOp requires this since we require static
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

  elseif(EXAGO_ENABLE_HIP)
    set(EXAGO_HAVE_HIP 1)
    add_definitions(-DHAVE_HIP) # For hipmagma
    find_package(HIP REQUIRED)
    set(CMAKE_CXX_COMPILER hipcc)
  else()
    message(
      FATAL_ERROR "EXAGO_ENABLE_GPU enabled, but neither EXAGO_ENABLE_CUDA "
                  "nor EXAGO_ENABLE_HIP is enabled!"
    )
  endif()
endif(EXAGO_ENABLE_GPU)

# Find PETSc and configure related variables
if(EXAGO_ENABLE_PETSC)
  include(FindPETSC)
  include_directories(${PETSC_INCLUDE_CONF})
endif()

# Set install rpath to the locations where EXAGO and PETSc libraries reside.
# TODO: Automatically skip if PETSc and EXAGO are already on system lib path.
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib ${PETSC_DIR}/lib)

if(EXAGO_ENABLE_IPOPT)
  set(IPOPT_ROOT_DIR CACHE PATH "Path to Ipopt installation root.")
  include(FindIpopt)
  set(CMAKE_INSTALL_RPATH ${IPOPT_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
endif(EXAGO_ENABLE_IPOPT)

if(EXAGO_ENABLE_HIOP)
  set(HIOP_ROOT_DIR CACHE PATH "Path to HiOp installation root.")
  include(FindHiop)
  set(CMAKE_INSTALL_RPATH ${HIOP_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
endif(EXAGO_ENABLE_HIOP)

if(EXAGO_ENABLE_RAJA)

  if(EXAGO_ENABLE_OMP)
    find_package(OpenMP REQUIRED)
  endif()

  # Look for CMake configuration file in RAJA installation
  find_package(
    RAJA CONFIG PATHS ${RAJA_DIR} ${RAJA_DIR}/share/raja/cmake REQUIRED
  )
  # Look for CMake configuration file in Umpire installation
  find_package(
    umpire CONFIG PATHS ${umpire_DIR} ${umpire_DIR}/share/umpire/cmake REQUIRED
  )
  target_include_directories(umpire INTERFACE ${camp_INSTALL_PREFIX}/include)
endif(EXAGO_ENABLE_RAJA)

if(NOT DEFINED LAPACK_LIBRARIES)
  # in case the toolchain defines them
  find_package(LAPACK REQUIRED)
  message(STATUS "Found LAPACK libraries: ${LAPACK_LIBRARIES}")
endif(NOT DEFINED LAPACK_LIBRARIES)

# Set up configuration header file
configure_file(include/exago_config.h.in exago_config.h)

#[[#############################################################################

Build libraries and testing framework.

#]]

if(EXAGO_RUN_TESTS)
  enable_testing()
  include(ExaGOTestingUtilities)
endif()

include_directories(${PROJECT_SOURCE_DIR}/tpl/spdlog/include)

# Build libraries (specified in src/CMakeLists.txt)
add_subdirectory(src)

# Build applications (specified in applications/CMakeLists.txt)
add_subdirectory(applications)

#[[#############################################################################

Installation of libraries, headers, and datafiles used for testing.

#]]

# EXAGO options files to install
set(EXAGO_OPTIONS_FILES
    options/pflowoptions
    options/opflowoptions
    options/scopflowoptions
    options/sopflowoptions
    options/tcopflowoptions
    options/hiop.options
    options/ipopt.opt
)

# Install options files
install(FILES ${EXAGO_OPTIONS_FILES} DESTINATION ${EXAGO_OPTIONS_DIR})
# Install options file to binary dir for testing
file(INSTALL ${EXAGO_OPTIONS_FILES} DESTINATION options)

# EXAGO datafiles to install
set(EXAGO_DATA_FILES
    ${PROJECT_SOURCE_DIR}/datafiles/case118.m
    ${PROJECT_SOURCE_DIR}/datafiles/case118.cont
    ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
)

# Install data files
install(FILES ${EXAGO_DATA_FILES} DESTINATION ${EXAGO_DATAFILES_DIR})
install(DIRECTORY datafiles/case9 DESTINATION ${EXAGO_DATAFILES_DIR})
install(DIRECTORY datafiles/TAMU200_scenarios
        DESTINATION ${EXAGO_DATAFILES_DIR}
)
install(DIRECTORY datafiles/test_validation DESTINATION ${EXAGO_DATAFILES_DIR})

# Install data files in build directory for testing
file(INSTALL ${EXAGO_DATA_FILES} DESTINATION datafiles)
file(INSTALL datafiles/case9 DESTINATION datafiles)
file(INSTALL datafiles/TAMU200_scenarios DESTINATION datafiles)
# Custom testing data files
file(INSTALL datafiles/test_validation DESTINATION datafiles)

# EXAGO headers to install
set(EXAGO_COMMON_INCLUDE include/constants.h include/version.hpp
                         ${PROJECT_BINARY_DIR}/exago_config.h
)

# Install header files
install(FILES ${EXAGO_COMMON_INCLUDE} DESTINATION include)

#[[#############################################################################

Test definitions.

@note The tests rely on utilities in @file cmake/ExaGOTestingUtilities.cmake.

@note If running with BSUB, we should use jsrun commands
  Use in CMake as:
    - ${MPICMD} "-n" "<num_proc>" ${EXAGO_EXTRA_MPI_FLAGS}
      for running with more than one process.
        OR
    - ${RUNCMD} when running with only one.

@note EXAGO_EXTRA_MPI_FLAGS should be passed as a ';' delimited string, e.g.
  export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS
  -DEXAGO_EXTRA_MPI_FLAGS=-mca;oob;tcp" will result in "-mca" "oob" "tcp" being
  passed to mpiexec/jsrun

#]]

if(EXAGO_TEST_WITH_BSUB)
  set(MPICMD "jsrun")
  set(RUNCMD "jsrun" "-n" "1")
  if(EXAGO_ENABLE_GPU)
    set(MPICMD ${MPICMD} "-g" "1")
    set(RUNCMD ${RUNCMD} "-g" "1" ${EXAGO_EXTRA_MPI_FLAGS})
  endif()
  # Otherwise, set mpi command + extra flags
elseif(EXAGO_TEST_WITH_SLURM)
  set(MPICMD "srun")
  set(RUNCMD "srun" "-n" "1")
  if(EXAGO_ENABLE_GPU)
    set(MPICMD ${MPICMD} "--gres=gpu:1")
    set(RUNCMD ${RUNCMD} "--gres=gpu:1" ${EXAGO_EXTRA_MPI_FLAGS})
  endif()
else()
  if(EXAGO_ENABLE_MPI)
    set(MPICMD "mpiexec")
    set(RUNCMD ${MPICMD} "-n" "1" ${EXAGO_EXTRA_MPI_FLAGS})
  else()
    set(MPICMD "")
    set(RUNCMD "")
  endif()
endif()

include(ExaGOCheckPython)
if(EXAGO_ENABLE_PYTHON)
  # Build pybind11 target manually
  add_subdirectory(tpl/pybind11)
  add_subdirectory(interfaces/python)
endif()

# Build tests
add_subdirectory(tests)

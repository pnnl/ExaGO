#[[#############################################################################

Preamble

#]]

cmake_minimum_required(VERSION 3.10)

project(ExaGO VERSION "1.0.0")

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

#[[#############################################################################

Configuration and library discovery.

#]]

# Specify where CMake modules are to be found
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/buildsystem/cmake)

# Set the path to options file in EXAGO installation
set(EXAGO_OPTIONS_DIR "${CMAKE_INSTALL_PREFIX}/share/exago/options")
set(EXAGO_DATAFILES_DIR "${CMAKE_INSTALL_PREFIX}/share/exago/datafiles")

# Specify whether to build shared libraries, static libraries, or both
option(EXAGO_BUILD_SHARED "Build shared libraries" ON)
option(EXAGO_BUILD_STATIC "Build static libraries" OFF)
include(ExaGOAddLibrary)

# Enable PETSc by default
option(EXAGO_ENABLE_PETSC "Enable PETSc support" ON)

# Enable MPI
option(EXAGO_ENABLE_MPI "Enable MPI support" ON)

# Ipopt support is disabled by default
option(EXAGO_ENABLE_IPOPT "Enable Ipopt support" OFF)

# HiOP support is disabled by default
option(EXAGO_ENABLE_HIOP "Enable HiOP support" OFF)

option(EXAGO_ENABLE_HIOP_SPARSE "Enable HiOp's HIOP_SPARSE option" OFF)

# MPI support is enabled by default
option(EXAGO_ENABLE_MPI "Enable MPI support" ON)

# GPU support is disabled by default
option(EXAGO_ENABLE_GPU "Enable GPU support" OFF)

# Enable RAJA support
option(EXAGO_ENABLE_RAJA "Enables RAJA and Umpire libraries" OFF)

# Enable tests
option(EXAGO_RUN_TESTS "Enable tests" ON)

# If tests are enabled, decides if tests are installed
option(EXAGO_INSTALL_TESTS "Install testing binaries" ON)

# For testing with BSUB on certain clusters
option(EXAGO_TEST_WITH_BSUB
       "Use `jsrun` instead of `mpirun` commands when running tests" OFF
)

option(EXAGO_DISABLE_LOGGING "Disable all internal logging" OFF)

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
  include(CheckLanguage)
  enable_language(CUDA)
  check_language(CUDA)

  find_package(CUDA REQUIRED)

  if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  endif()

  if(HIOP_NVCC_ARCH)
    set(CMAKE_CUDA_FLAGS "-arch=${HIOP_NVCC_ARCH}")
    message(STATUS "Using CUDA arch ${HIOP_NVCC_ARCH}")
  else()
    set(CMAKE_CUDA_FLAGS "-arch=sm_35")
    message(STATUS "Using CUDA arch sm_35")
  endif()
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
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
  set(HIOP_ROOT_DIR CACHE PATH "Path to Ipopt installation root.")
  include(FindHiop)
  set(CMAKE_INSTALL_RPATH ${HIOP_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
endif(EXAGO_ENABLE_HIOP)

if(EXAGO_ENABLE_RAJA)
  # Find OpenMP -- this is a temporary fix
  find_package(OpenMP)
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
set(EXAGO_COMMON_INCLUDE include/constants.h include/version.h
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
else()
  set(MPICMD "mpiexec")
  set(RUNCMD "mpiexec" "-n" "1")
  set(MPICMD ${MPICMD})
  set(RUNCMD ${RUNCMD} ${EXAGO_EXTRA_MPI_FLAGS})
endif()

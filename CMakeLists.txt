#[[#############################################################################

Preamble

#]]

cmake_minimum_required(VERSION 3.18)

project(ExaGO VERSION "1.5.1")

set(CMAKE_CXX_STANDARD 14)
# Enable setting Mac OS X rpath
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # CMake 3.12
endif()

if(POLICY CMP0104)
  cmake_policy(SET CMP0104 NEW) # CMake 3.18
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING ""
  )
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
option(EXAGO_ENABLE_IPOPT "Enable Ipopt support" ON)

# HiOP support is disabled by default
option(EXAGO_ENABLE_HIOP "Enable HiOP support" OFF)

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

set(EXAGO_CTEST_LAUNCH_COMMAND
    ""
    CACHE STRING "Command to use when launching tests"
)

option(EXAGO_ENABLE_LOGGING "Enable internal logging" OFF)

# When building with HIP support, OpenMP is not supported.
if(EXAGO_ENABLE_RAJA AND EXAGO_ENABLE_HIP)
  set(EXAGO_ENABLE_OMP OFF)
endif()
# When building without RAJA, OpenMP is not required.
if(NOT EXAGO_ENABLE_RAJA)
  set(EXAGO_ENABLE_OMP OFF)
endif()

# When building with OpenMP using Apple-Clang, users should be warned about
# numerical precision issues
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  message(
    WARNING
      " When building with apple-clang built OpenMP, you will get numerical precision issues"
  )
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
    set(EXAGO_ENABLE_CUDA
        ON
        CACHE BOOL ""
    )
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
    find_package(hip REQUIRED)
  else()
    message(
      FATAL_ERROR "EXAGO_ENABLE_GPU enabled, but neither EXAGO_ENABLE_CUDA "
                  "nor EXAGO_ENABLE_HIP is enabled!"
    )
  endif()
endif(EXAGO_ENABLE_GPU)

# Find PETSc and configure related variables
if(EXAGO_ENABLE_PETSC)
  include(FindPkgConfig)
  # Include petsc package path in pkg_config_path
  set(ENV{PKG_CONFIG_PATH}
      ${PETSC_DIR}/lib/pkgconfig:${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig
  )
  pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
  set(EXAGO_HAVE_PETSC 1)
endif()

# Set install rpath to the locations where EXAGO and PETSc libraries reside.
# TODO: Automatically skip if PETSc and EXAGO are already on system lib path.
# set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib ${PETSC_DIR}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

if(EXAGO_ENABLE_IPOPT)
  set(IPOPT_ROOT_DIR CACHE PATH "Path to Ipopt installation root.")
  include(FindIpopt)
  set(CMAKE_INSTALL_RPATH ${IPOPT_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
endif(EXAGO_ENABLE_IPOPT)

if(EXAGO_ENABLE_OMP)
  find_package(OpenMP REQUIRED COMPONENTS C CXX)
endif()

if(EXAGO_ENABLE_HIOP)
  include(ExaGOFindHiOp)
  if(NOT TARGET RAJA)
    # In HiOp's exported CMake configuration, there is a call to:
    # find_package(RAJA PATHS @RAJA_DIR@) This means the RAJA target should be
    # defined in HiOp 0.7.0+ if hiop+raja RAJA part of ExaGO will not work if
    # HiOp is not built with RAJA
    message(
      STATUS
        "Warning: HiOp was built without RAJA - ExaGO RAJA module calling HiOp may not work!"
    )
  endif()
endif(EXAGO_ENABLE_HIOP)

if(EXAGO_ENABLE_RAJA)

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
configure_file(include/exago_build_options.h.in exago_build_options.h)

#[[#############################################################################

Build libraries and testing framework.

#]]

if(EXAGO_RUN_TESTS)
  enable_testing()
  include(ExaGOTestingUtilities)
endif()

if(EXAGO_ENABLE_LOGGING)
  add_subdirectory(${PROJECT_SOURCE_DIR}/tpl/spdlog/)
  # -fPIC required for spdlog for use with python wrapper
  set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Build libraries (specified in src/CMakeLists.txt)
add_subdirectory(src)

# Build applications (specified in applications/CMakeLists.txt)
add_subdirectory(applications)

#[[#############################################################################

Installation of libraries, headers, and datafiles used for testing.

#]]

# EXAGO options files to install
set(EXAGO_OPTIONS_FILES options/pflowoptions)

if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)
  set(EXAGO_OPTIONS_FILES ${EXAGO_OPTIONS_FILES} options/opflowoptions
                          options/hiop.options
  )
endif()

if(EXAGO_ENABLE_IPOPT)
  set(EXAGO_OPTIONS_FILES
      ${EXAGO_OPTIONS_FILES} options/scopflowoptions options/sopflowoptions
      options/tcopflowoptions options/ipopt.opt
  )
endif()

# Install options files
install(FILES ${EXAGO_OPTIONS_FILES} DESTINATION ${EXAGO_OPTIONS_DIR})
# Install options file to binary dir for testing
file(INSTALL ${EXAGO_OPTIONS_FILES} DESTINATION options)

# EXAGO datafiles to install
set(EXAGO_DATA_FILES
    ${PROJECT_SOURCE_DIR}/datafiles/case118.m
    ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
    ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg10k.m
)

if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)
  set(EXAGO_DATA_FILES
      ${EXAGO_DATA_FILES}
      ${PROJECT_SOURCE_DIR}/datafiles/case118.cont
      ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
      ${PROJECT_SOURCE_DIR}/datafiles/10_scenarios_ACTIVSg200.csv
      ${PROJECT_SOURCE_DIR}/datafiles/10_scenarios_ACTIVSg2000.csv
      ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg10k.cont
  )
endif()

# Install data files
install(FILES ${EXAGO_DATA_FILES} DESTINATION ${EXAGO_DATAFILES_DIR})

if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)
  install(DIRECTORY datafiles/unit DESTINATION ${EXAGO_DATAFILES_DIR})
  install(DIRECTORY datafiles/TAMU200_scenarios
          DESTINATION ${EXAGO_DATAFILES_DIR}
  )
  install(DIRECTORY datafiles/test_validation
          DESTINATION ${EXAGO_DATAFILES_DIR}
  )
endif()

install(DIRECTORY datafiles/case9 DESTINATION ${EXAGO_DATAFILES_DIR})

# Install data files in build directory for testing
file(INSTALL ${EXAGO_DATA_FILES} DESTINATION datafiles)
file(INSTALL datafiles/case9 DESTINATION datafiles)

if(EXAGO_ENABLE_IPOPT OR EXAGO_ENABLE_HIOP)
  file(INSTALL datafiles/TAMU200_scenarios DESTINATION datafiles)
  file(INSTALL datafiles/unit DESTINATION datafiles)
  # Custom testing data files
  file(INSTALL datafiles/test_validation DESTINATION datafiles)
endif()

# EXAGO headers to install
set(EXAGO_COMMON_INCLUDE
    include/constants.h include/version.h ${PROJECT_BINARY_DIR}/exago_config.h
    ${PROJECT_BINARY_DIR}/exago_build_options.h
)

# Install header files
install(FILES ${EXAGO_COMMON_INCLUDE} DESTINATION include)

#[[#############################################################################

Test definitions.

@note The tests rely on utilities in @file cmake/ExaGOTestingUtilities.cmake.

@note If running with a job scheduler, pass the required commands to
  EXAGO_CTEST_LAUNCH_COMMAND. For example, to launch a slurm job for each test:
  `cmake $SRC_DIR -DEXAGO_CTEST_LAUNCH_COMMAND="srun -N 1 --gres=gpu:1 -t 10:00"`

#]]
if(NOT (EXAGO_CTEST_LAUNCH_COMMAND STREQUAL ""))
  message(
    STATUS
      "Manually overriding launch command for CTest with \"${EXAGO_CTEST_LAUNCH_COMMAND}\""
  )
  string(REPLACE " " ";" RUNCMD_LIST "${EXAGO_CTEST_LAUNCH_COMMAND}")
  set(RUNCMD ${RUNCMD_LIST})
  set(MPICMD ${RUNCMD_LIST})
else()
  if(EXAGO_ENABLE_MPI)
    set(MPICMD "${MPIEXEC_EXECUTABLE}")
    set(RUNCMD ${MPICMD} "${MPIEXEC_NUMPROC_FLAG}" "1" ${EXAGO_EXTRA_MPI_FLAGS})
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

if((NOT EXAGO_ENABLE_IPOPT) AND (NOT EXAGO_ENABLE_HIOP))
  message(
    WARNING
      "ExaGO not configured with any optimization solver. Only PFLOW application will be available. Please enable Ipopt and/or HiOp for other applications."
  )
endif()

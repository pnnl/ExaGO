
cmake_minimum_required(VERSION 3.10)

project(Scopflow)

## Preamble

# Enable setting Mac OS X rpath
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # CMake 3.12
endif ()


## Configure

# Specify where CMake modules are to be found
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set the path to options file in Scopflow installation
set(SCOPFLOW_OPTIONS_DIR "${CMAKE_INSTALL_PREFIX}/options")

# Enable PETSc by default
option(SCOPFLOW_ENABLE_PETSC "Enable PETSc support" ON)

# Ipopt support is disabled by default
option(SCOPFLOW_ENABLE_IPOPT "Enable Ipopt support" OFF)

# HiOP support is disabled by default
option(SCOPFLOW_ENABLE_HIOP "Enable HiOP support" OFF)

# Enable tests
option(SCOPFLOW_RUN_TESTS "Enable tests" OFF)

# Set default include path for the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})

# Create alias for -lm which is needed when linking math on UNIX
if(UNIX)
  set(SCOPFLOW_MATH_LIB m)
endif()

find_package(MPI REQUIRED COMPONENTS C CXX)

# Find PETSc and configure related variables
if(SCOPFLOW_ENABLE_PETSC)
  include(FindPETSC)
  include_directories(${PETSC_INCLUDE_CONF})
endif()

# Set install rpath to the locations where Scopflow and PETSc libraries reside.
# TODO: Automatically skip if PETSc and Scopflow are already on system lib path.
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib ${PETSC_DIR}/lib)

if(SCOPFLOW_ENABLE_IPOPT)
  set(IPOPT_ROOT_DIR CACHE PATH "Path to Ipopt installation root.")
  include(FindIpopt)
  set(CMAKE_INSTALL_RPATH  ${IPOPT_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
  set(SCOPFLOW_HAVE_IPOPT 1)
endif(SCOPFLOW_ENABLE_IPOPT)

if(SCOPFLOW_ENABLE_HIOP)
  set(HIOP_ROOT_DIR CACHE PATH "Path to Ipopt installation root.")
  include(FindHiop)
  set(CMAKE_INSTALL_RPATH  ${HIOP_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
  set(SCOPFLOW_HAVE_HIOP 1)
endif(SCOPFLOW_ENABLE_HIOP)

if(NOT DEFINED LAPACK_LIBRARIES)
  # in case the toolchain defines them
  find_package(LAPACK REQUIRED)
  message(STATUS "Found LAPACK libraries: ${LAPACK_LIBRARIES}")
endif(NOT DEFINED LAPACK_LIBRARIES)

# Set up configuration header file
configure_file(include/scopflow_config.h.in scopflow_config.h)

## Build

# Build libraries (specified in src/CMakeLists.txt)
add_subdirectory(src)

# Build applications (specified in applications/CMakeLists.txt)
add_subdirectory(applications)

# Build unit tests
add_subdirectory(tests)

## Install

# Scopflow options files to install
set(SCOPFLOW_OPTIONS_FILES
  options/pflowoptions
  options/opflowoptions
  options/scopflowoptions
  options/tcopflowoptions
)

# Install options files
install(FILES ${SCOPFLOW_OPTIONS_FILES} DESTINATION options)

# Scopflow datafiles to install
set(SCOPFLOW_DATA_FILES
  datafiles/case9mod.m
  datafiles/case118.m
  datafiles/case118.cont
  datafiles/case_ACTIVSg200.m
  datafiles/case9.cont
  datafiles/case9.con
)

# Install data files
install(FILES ${SCOPFLOW_DATA_FILES} DESTINATION bin/datafiles)
install(FILES ${SCOPFLOW_DATA_FILES} DESTINATION datafiles)
install(FILES ${SCOPFLOW_DATA_FILES} DESTINATION tests/datafiles)
install(DIRECTORY datafiles/TAMU200_scenarios DESTINATION bin/datafiles)
install(DIRECTORY datafiles/TAMU200_scenarios DESTINATION datafiles)

# Scopflow headers to install
set(SCOPFLOW_COMMON_INCLUDE
  include/constants.h
  "${PROJECT_BINARY_DIR}/scopflow_config.h"
)
 
# Install header files
install(FILES ${SCOPFLOW_COMMON_INCLUDE} DESTINATION include)

if(SCOPFLOW_RUN_TESTS)
  # Testing
  enable_testing()

  file(GLOB network_files
    "${PROJECT_BINARY_DIR}/datafiles/*.m"
  )
  message(STATUS "Starting tests")

  # Tests with IPOPT
  if(SCOPFLOW_HAVE_IPOPT)
    # SCOPFLOW tests - Preventive mode
    foreach(nc RANGE 1 9)
      add_test(NAME SCOPFLOW_9bus_${nc}cont_preventive COMMAND mpiexec -n 1 $<TARGET_FILE:app_scopflow> -netfile datafiles/case9mod.m -ctgcfile datafiles/case9.cont -scopflow_solver IPOPT -scopflow_mode 0 -opflow_include_loadloss_variables -scopflow_Nc ${nc})
    endforeach(nc)

    # SCOPFLOW tests - Corrective mode
    foreach(nc RANGE 1 9)
      add_test(NAME SCOPFLOW_9bus_${nc}cont_corrective COMMAND mpiexec -n 1 $<TARGET_FILE:app_scopflow> -netfile datafiles/case9mod.m -ctgcfile datafiles/case9.cont -scopflow_solver IPOPT -scopflow_mode 1 -opflow_include_loadloss_variables -scopflow_Nc ${nc})
    endforeach(nc)

    message(STATUS "Running OPFLOW formulation tests")

    # OPFLOW tests

    foreach(test_net ${network_files})
      foreach(opflow_form POWER_BALANCE_POLAR)
        foreach(opflow_init FROMFILE FLATSTART ACPF MIDPOINT)
          string(REPLACE "${PROJECT_BINARY_DIR}/datafiles/" "" input_file ${test_net})
          add_test(NAME OPFLOW_${opflow_form}_${opflow_init}_${input_file} COMMAND mpiexec -n 1 $<TARGET_FILE:app_opflow> -netfile ${test_net} -opflow_formulation ${opflow_form} -opflow_initialization ${opflow_init} -opflow_solver IPOPT)
        endforeach()
      endforeach()
    endforeach()
  endif(SCOPFLOW_HAVE_IPOPT)

  # Tests with HIOP
  if(SCOPFLOW_HAVE_HIOP)
    foreach(test_net ${network_files})
      foreach(opflow_init FROMFILE FLATSTART ACPF MIDPOINT)
        string(REPLACE "${PROJECT_BINARY_DIR}/datafiles/" "" input_file ${test_net})
        add_test(NAME OPFLOW_HIOP_${opflow_init}_${input_file} COMMAND mpiexec -n 1 $<TARGET_FILE:app_opflow> -netfile ${test_net} -opflow_initialization ${opflow_init} -opflow_solver HIOP)
      endforeach()
    endforeach()
  endif(SCOPFLOW_HAVE_HIOP)

  # Tests with TAO solver
  foreach(np RANGE 1 1)
    # Tests for > 1 proc are failing because DMNetwork in PS needs to be updated
    foreach(test_net ${network_files})
      string(REPLACE "${PROJECT_BINARY_DIR}/" "" input_file ${test_net})
      add_test(NAME OPFLOW_TAO_${input_file} COMMAND mpiexec -n 1 $<TARGET_FILE:app_opflow> -netfile ${test_net} -opflow_solver TAO)
    endforeach(test_net)
  endforeach(np)

  message(STATUS "Running PFLOW tests")
  # PFLOW tests
  foreach(np RANGE 1 2)
    foreach(test_net ${network_files})
      string(REPLACE "${PROJECT_BINARY_DIR}/" "" input_file ${test_net})
      add_test(NAME PFLOW_${input_file}_${np}proc COMMAND mpiexec -n 1 $<TARGET_FILE:app_pflow> -netfile ${test_net})
    endforeach(test_net)
  endforeach(np)

endif(SCOPFLOW_RUN_TESTS)

cmake_minimum_required(VERSION 3.10)

project(Scopflow)

## Preamble

# Enable setting Mac OS X rpath
set(CMAKE_MACOSX_RPATH 1)

## Configure

# Specify where CMake modules are to be found
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set the path to options file in Scopflow installation
set(SCOPFLOW_OPTIONS_DIR "${CMAKE_INSTALL_PREFIX}/options")

# Set up configuration header file
configure_file(include/scopflow_config.h.in scopflow_config.h)

# Enable PETSc by default
option(SCOPFLOW_ENABLE_PETSC "Enable PETSc support" ON)

# Ipopt support is disabled by default
option(SCOPFLOW_ENABLE_IPOPT "Enable Ipopt support" OFF)

# Set default include path for the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})

# Create alias for -lm which is needed when linking math on UNIX
if(UNIX)
  set(SCOPFLOW_MATH_LIB m)
endif()

# Variable to store power flow libraries
set(SCOPFLOW_PFLOW_LIBS
  #
)

find_package(MPI REQUIRED COMPONENTS C CXX)

# Find PETSc and configure related variables
if(SCOPFLOW_ENABLE_PETSC)
  include(FindPETSC)
  include_directories(${PETSC_INCLUDE_CONF})
endif()

if(SCOPFLOW_ENABLE_IPOPT)
  include(FindIpopt)
endif()

# Set install rpath to the locations where Scopflow and PETSc libraries reside.
# TODO: Automatically skip if PETSc and Scopflow are already on system lib path.
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib ${PETSC_DIR}/lib)


## Build

# Build libraries (specified in src/CMakeLists.txt)
add_subdirectory(src)

# Build applications (specified in applications/CMakeLists.txt)
add_subdirectory(applications)


## Install

# Scopflow options files to install (note difference between strings and vars)
set(SCOPFLOW_OPTIONS_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/options/pflowoptions"
  options/opflowoptions
  options/scopflowoptions
)

# Install options files
install(FILES ${SCOPFLOW_OPTIONS_FILES} DESTINATION options)


# Scopflow headers to install
set(SCOPFLOW_COMMON_INCLUDE
  include/common.h
  include/constants.h
  include/ps.h
  "${PROJECT_BINARY_DIR}/scopflow_config.h"
)
 
# Install header files
install(FILES ${SCOPFLOW_COMMON_INCLUDE} DESTINATION include)

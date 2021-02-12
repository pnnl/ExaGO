
cmake_minimum_required(VERSION 3.10)

project(ExaGO VERSION "0.0.1")

string(TIMESTAMP EXAGO_RELEASE_DATE "%Y-%m-%d")

## Preamble

set(CMAKE_CXX_STANDARD 11)
# Enable setting Mac OS X rpath
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # CMake 3.12
endif ()

if (POLICY CMP0104)
  cmake_policy(SET CMP0104 NEW) # CMake 3.18
endif ()

## Configure

#set(MPI_ASSUME_NO_BUILTIN_MPI On)

# Specify where CMake modules are to be found
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set the path to options file in EXAGO installation
set(EXAGO_OPTIONS_DIR "${CMAKE_INSTALL_PREFIX}/options")

# Enable PETSc by default
option(EXAGO_HAVE_PETSC "Enable PETSc support" ON)

# Enable MPI
option(EXAGO_HAVE_MPI "Enable MPI support" ON)

# Ipopt support is disabled by default
option(EXAGO_HAVE_IPOPT "Enable Ipopt support" OFF)

# HiOP support is disabled by default
option(EXAGO_HAVE_HIOP "Enable HiOP support" OFF)

option(EXAGO_HAVE_HIOP_SPARSE "Enable HiOp's HIOP_SPARSE option" OFF)

# MPI support is enabled by default
option(EXAGO_HAVE_MPI "Enable MPI support" ON)

# GPU support is disabled by default
option(EXAGO_HAVE_GPU "Enable GPU support" OFF)

# Enable RAJA support
option(EXAGO_HAVE_RAJA "Enables RAJA and Umpire libraries" OFF)

# Enable tests
option(EXAGO_RUN_TESTS "Enable tests" OFF)

# For testing with BSUB on certain clusters
option(EXAGO_TEST_WITH_BSUB "Use `jsrun` instead of `mpirun` commands when running tests" OFF)

option(EXAGO_DISABLE_LOGGING "Disable all internal logging" OFF)

# Estimate number of flops per trig op
set(EXAGO_FLOPS_SINOP 15 CACHE STRING "Estimated number of flops per sin op")
set(EXAGO_FLOPS_COSOP 15 CACHE STRING "Estimated number of flops per cos op")

# Set default include path for the project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR})

# Create alias for -lm which is needed when linking math on UNIX
if(UNIX)
  set(EXAGO_MATH_LIB m)
endif()

# If MPI support is enabled, find MPI includes and libraries
if(EXAGO_HAVE_MPI)
  find_package(MPI REQUIRED COMPONENTS C CXX)
  if(NOT DEFINED MPI_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  endif(NOT DEFINED MPI_CXX_COMPILER)
  include_directories(${MPI_CXX_ADDITIONAL_INCLUDE_DIRS} ${MPI_CXX_COMPILER_INCLUDE_DIRS})
  set(EXAGO_HAVE_MPI 1)
  set(EXAGO_EXTRA_MPI_FLAGS "" CACHE STRING "Extra args to mpiexec when running tests")
endif(EXAGO_HAVE_MPI)

# If GPU support is enabled, find CUDA
if(EXAGO_HAVE_GPU)
  include(CheckLanguage)
  enable_language(CUDA)
  check_language(CUDA)
  
  find_package(CUDA REQUIRED)

  set(EXAGO_HAVE_GPU 1)

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
endif(EXAGO_HAVE_GPU)


# Find PETSc and configure related variables
if(EXAGO_HAVE_PETSC)
  include(FindPETSC)
  include_directories(${PETSC_INCLUDE_CONF})
  set(EXAGO_HAVE_PETSC 1)
endif()

# Set install rpath to the locations where EXAGO and PETSc libraries reside.
# TODO: Automatically skip if PETSc and EXAGO are already on system lib path.
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib ${PETSC_DIR}/lib)

if(EXAGO_HAVE_IPOPT)
  set(IPOPT_ROOT_DIR CACHE PATH "Path to Ipopt installation root.")
  include(FindIpopt)
  set(CMAKE_INSTALL_RPATH  ${IPOPT_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
  set(EXAGO_HAVE_IPOPT 1)
endif(EXAGO_HAVE_IPOPT)

if(EXAGO_HAVE_HIOP)
  set(HIOP_ROOT_DIR CACHE PATH "Path to Ipopt installation root.")
  include(FindHiop)
  set(CMAKE_INSTALL_RPATH  ${HIOP_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
  set(EXAGO_HAVE_HIOP 1)
endif(EXAGO_HAVE_HIOP)

if(EXAGO_HAVE_RAJA)
  # Find OpenMP -- this is a temporary fix
  find_package(OpenMP)
  # Look for CMake configuration file in RAJA installation
  find_package(RAJA CONFIG
    PATHS ${RAJA_DIR} ${RAJA_DIR}/share/raja/cmake
    REQUIRED)
  # Look for CMake configuration file in Umpire installation
  find_package(umpire CONFIG
    PATHS ${umpire_DIR} ${umpire_DIR}/share/umpire/cmake
    REQUIRED)
  set(EXAGO_HAVE_RAJA 1)
endif(EXAGO_HAVE_RAJA)

if(NOT DEFINED LAPACK_LIBRARIES)
  # in case the toolchain defines them
  find_package(LAPACK REQUIRED)
  message(STATUS "Found LAPACK libraries: ${LAPACK_LIBRARIES}")
endif(NOT DEFINED LAPACK_LIBRARIES)

# Set up configuration header file
configure_file(include/exago_config.h.in exago_config.h)

## Build

if(EXAGO_RUN_TESTS)
  enable_testing()
endif()

# Build libraries (specified in src/CMakeLists.txt)
add_subdirectory(src)

# Build applications (specified in applications/CMakeLists.txt)
add_subdirectory(applications)

# Build tests
add_subdirectory(tests)

## Install

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
install(FILES ${EXAGO_OPTIONS_FILES} DESTINATION options)

# EXAGO datafiles to install
set(EXAGO_DATA_FILES
  ${PROJECT_SOURCE_DIR}/datafiles/case118.m
  ${PROJECT_SOURCE_DIR}/datafiles/case118.cont
  ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg200.m
  ${PROJECT_SOURCE_DIR}/datafiles/case_ACTIVSg200.cont
)

# Install data files
install(FILES ${EXAGO_DATA_FILES} DESTINATION bin/datafiles)
install(FILES ${EXAGO_DATA_FILES} DESTINATION datafiles)
install(FILES ${EXAGO_DATA_FILES} DESTINATION tests/datafiles)
install(DIRECTORY datafiles/case9 DESTINATION bin/datafiles)
install(DIRECTORY datafiles/TAMU200_scenarios DESTINATION bin/datafiles)
install(DIRECTORY datafiles/TAMU200_scenarios DESTINATION datafiles)

set(EXAGO_SRC_DATAFILES_DIR ${PROJECT_SOURCE_DIR}/datafiles)
set(EXAGO_SRC_OPTIONS_DIR   ${PROJECT_SOURCE_DIR}/options)

# Custom testing data files
install(DIRECTORY datafiles/test_validation DESTINATION tests/datafiles)

# EXAGO headers to install
set(EXAGO_COMMON_INCLUDE
  include/constants.h
  include/version.h
  ${PROJECT_BINARY_DIR}/exago_config.h
)
 
# Install header files
install(FILES ${EXAGO_COMMON_INCLUDE} DESTINATION include)

# If running with BSUB, we should use jsrun commands
# Use in CMake as:
#   - ${MPICMD} "-n" "<num_proc>" ${EXAGO_EXTRA_MPI_FLAGS}
#     for running with more than one process.
#       OR
#   - ${RUNCMD} when running with only one.
#
# Note : ${EXAGO_EXTRA_MPI_FLAGS} should be passed as a ';' delimited string
# e.g. export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_EXTRA_MPI_FLAGS=-mca;oob;tcp"
#      will result in "-mca" "oob" "tcp" being passed to mpiexec/jsrun
if(EXAGO_TEST_WITH_BSUB)
  set(MPICMD "jsrun")
  set(RUNCMD "jsrun" "-n" "1")
  if(EXAGO_HAVE_GPU)
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

if(EXAGO_RUN_TESTS)

  # Test cases
  file(GLOB network_files
    datafiles/case9/case9mod.m
    datafiles/case118.m
    datafiles/case_ACTIVSg200.m
    )
  message(STATUS "Starting tests")

  # Unit tests
  if(EXAGO_HAVE_HIOP)
    foreach(test_net ${network_files})
      get_filename_component(input_file ${test_net} NAME)
      add_test(NAME UNIT_TESTS_OPFLOW_MODEL_${input_file}_TEST COMMAND
        ${RUNCMD} $<TARGET_FILE:testAcopf>
        -netfile ${test_net})
    endforeach()
  else()
    #[[ Some of our automated tests only run unit tests with ctest arguments
    `-R UNIT_TESTS`. We have this placeholder such that our ctest invokation
    will fail if only unit tests are ran and no unit tests were available. #]]
    add_test(NAME UNIT_TESTS_PLACEHOLDER_NOHIOP
      COMMAND bash -c "
      echo 'No unit tests exist without HIOP enabled.'
      exit 1")
  endif()

  # Tests with IPOPT
  if(EXAGO_HAVE_IPOPT)

    #  SCOPFLOW single-period tests

       foreach(nc RANGE 1 9)
         add_test(NAME SCOPFLOW_9bus_${nc}cont COMMAND 
                  ${RUNCMD} $<TARGET_FILE:app_scopflow> 
                  -netfile  ${EXAGO_SRC_DATAFILES_DIR}/case9/case9mod.m
                  -ctgcfile ${EXAGO_SRC_DATAFILES_DIR}/case9/case9.cont
                  -scopflow_solver IPOPTNEW
                  -scopflow_Nc ${nc})
         endforeach(nc)

    #  SCOPFLOW multi-period tests

       foreach(nc RANGE 1 9)
         add_test(NAME SCOPFLOW_9bus_multiperiod_${nc}cont COMMAND 
                  ${RUNCMD} $<TARGET_FILE:app_scopflow> 
                  -netfile  ${EXAGO_SRC_DATAFILES_DIR}/case9/case9mod_gen3_wind.m
                  -ctgcfile ${EXAGO_SRC_DATAFILES_DIR}/case9/case9.cont
		  -scopflow_enable_multiperiod 1
		  -scopflow_ploadprofile  ${EXAGO_SRC_DATAFILES_DIR}/case9/load_P.csv
		  -scopflow_qloadprofile  ${EXAGO_SRC_DATAFILES_DIR}/case9/load_Q.csv
		  -scopflow_windgenprofile  ${EXAGO_SRC_DATAFILES_DIR}/case9/scenarios.csv
                  -scopflow_solver IPOPTNEW
                  -scopflow_Nc ${nc})
         endforeach(nc)


   #   SOPFLOW tests

       foreach(ns RANGE 1 3)
         add_test(NAME SOPFLOW_200bus_IPOPT_${ns}scen COMMAND 
                  ${RUNCMD} $<TARGET_FILE:app_sopflow>
                 -netfile  ${EXAGO_SRC_DATAFILES_DIR}/case_ACTIVSg200.m
                 -scenfile ${EXAGO_SRC_DATAFILES_DIR}/TAMU200_scenarios/scenarios.csv
                 -sopflow_solver IPOPT
		 -sopflow_enable_multicontingency 0
                 -sopflow_Ns ${ns})
       endforeach(ns)

       foreach(ns RANGE 1 3)
         add_test(NAME SOPFLOW_200bus_EMPAR_${ns}scen COMMAND 
                 ${RUNCMD} $<TARGET_FILE:app_sopflow>
                 -netfile  ${EXAGO_SRC_DATAFILES_DIR}/case_ACTIVSg200.m
                 -scenfile ${EXAGO_SRC_DATAFILES_DIR}/TAMU200_scenarios/scenarios.csv
                 -sopflow_solver EMPAR
		 -sopflow_enable_multicontingency 0
                 -sopflow_Ns ${ns})
       endforeach(ns)

       #Multicontingency

       add_test(NAME SOPFLOW_9bus_3scen_9cont COMMAND 
                 ${RUNCMD} $<TARGET_FILE:app_sopflow>
                 -netfile  ${EXAGO_SRC_DATAFILES_DIR}/case9/case9mod_gen3_wind.m
                 -scenfile ${EXAGO_SRC_DATAFILES_DIR}/case9/scenarios.csv
		 -ctgcfile ${EXAGO_SRC_DATAFILES_DIR}/case9/case9.cont
                 -sopflow_solver IPOPTNEW
		 -sopflow_enable_multicontingency 1
                 -sopflow_Ns -1
		 -scopflow_Nc -1)


       #Multicontingency-multiperiod

       add_test(NAME SOPFLOW_9bus_3scen_9cont_3steps COMMAND 
                 ${RUNCMD} $<TARGET_FILE:app_sopflow>
                 -netfile  ${EXAGO_SRC_DATAFILES_DIR}/case9/case9mod_gen3_wind.m
                 -scenfile ${EXAGO_SRC_DATAFILES_DIR}/case9/scenarios.csv
		 -ctgcfile ${EXAGO_SRC_DATAFILES_DIR}/case9/case9.cont
		 -scopflow_ploadprofile ${EXAGO_SRC_DATAFILES_DIR}/case9/load_P.csv
		 -scopflow_qloadprofile ${EXAGO_SRC_DATAFILES_DIR}/case9/load_Q.csv
                 -sopflow_solver IPOPTNEW
		 -sopflow_enable_multicontingency 1
		 -scopflow_enable_multiperiod 1
		 -scopflow_dT 5.0
		 -scopflow_duration 0.16666667
                 -sopflow_Ns -1
		 -scopflow_Nc -1)

    message(STATUS "Running OPFLOW model tests")

    # OPFLOW tests

    foreach(test_net ${network_files})
      foreach(opflow_model POWER_BALANCE_POLAR)
        foreach(opflow_init FROMFILE FLATSTART ACPF MIDPOINT)
          foreach(opflow_obj MIN_GEN_COST MIN_GENSETPOINT_DEVIATION)
            get_filename_component(input_file ${test_net} NAME)
            add_test(NAME OPFLOW_${opflow_model}_${opflow_init}_${input_file}_${opflow_obj} COMMAND
                   ${RUNCMD} $<TARGET_FILE:app_opflow>
                   -netfile ${test_net}
                   -opflow_model ${opflow_model}
                   -opflow_initialization ${opflow_init}
                   -opflow_solver IPOPT
		   -opflow_objective ${opflow_obj})
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endif(EXAGO_HAVE_IPOPT)

  # Tests with HIOP
  if(EXAGO_HAVE_HIOP)
    foreach(test_net ${network_files})
      foreach(opflow_init ACPF MIDPOINT)
        get_filename_component(input_file ${test_net} NAME)
        add_test(NAME OPFLOW_HIOP_CPU_${opflow_init}_${input_file} COMMAND
                 ${RUNCMD} $<TARGET_FILE:app_opflow>
                 -netfile ${test_net}
                 -opflow_initialization ${opflow_init}
                 -hiop_compute_mode CPU
                 -opflow_solver HIOP
                 -opflow_ignore_lineflow_constraints 1
                 -opflow_model POWER_BALANCE_HIOP)
      endforeach()
    endforeach()

    foreach(test_net ${network_files})
      foreach(opflow_init ACPF MIDPOINT)
        get_filename_component(input_file ${test_net} NAME)
        set(hiop_compute_mode HYBRID)
        if(EXAGO_HAVE_GPU)
          set(hiop_compute_mode GPU)
        endif()
        add_test(NAME OPFLOW_HIOP_GPU_${opflow_init}_${input_file} COMMAND
                 ${RUNCMD} $<TARGET_FILE:app_opflow>
                 -netfile ${test_net}
                 -opflow_initialization ${opflow_init}
                 -hiop_compute_mode ${hiop_compute_mode}
                 -opflow_solver HIOP
                 -opflow_ignore_lineflow_constraints 1
                 -opflow_model PBPOLRAJAHIOP)
      endforeach()
    endforeach()
  endif(EXAGO_HAVE_HIOP)

  # Tests with TAO solver
  foreach(np RANGE 1 1)
    # Tests for > 1 proc are failing because DMNetwork in PS needs to be updated
    foreach(test_net ${network_files})
      get_filename_component(input_file ${test_net} NAME)
      add_test(NAME OPFLOW_TAO_${input_file} COMMAND
               ${RUNCMD} $<TARGET_FILE:app_opflow>
               -netfile ${test_net}
               -opflow_solver TAO)
    endforeach(test_net)
  endforeach(np)

  message(STATUS "Configuring OPFLOW functionality tests.")
  foreach(test_net ${network_files})
    get_filename_component(test_net_name ${test_net} NAME)
    set(tolerance "1e-3")

    set(solver HIOP)
    foreach(model POWER_BALANCE_HIOP PBPOLRAJAHIOP)
      if("OPFLOW_${solver}_${model}_${test_net_name}_FUNCTIONALITY_TEST"
          STREQUAL
          "OPFLOW_HIOP_POWER_BALANCE_HIOP_case_ACTIVSg200.m_FUNCTIONALITY_TEST")
        message(STATUS "Skipping test `OPFLOW_HIOP_POWER_BALANCE_HIOP_case_ACTI"
          "VSg200.m_FUNCTIONALITY_TEST`.")
        message(STATUS "See issue https://gitlab.pnnl.gov/exasgd/frameworks/exa"
          "go/-/issues/157 for more information.")
        continue()
      endif()
      add_test(NAME OPFLOW_${solver}_${model}_${test_net_name}_FUNCTIONALITY_TEST
        COMMAND ${RUNCMD} $<TARGET_FILE:testOpflowFunctionality>
          -opflow_solver ${solver}
          -opflow_model ${model}
          -netfile ${test_net}
          -opflow_tolerance ${tolerance}
        )
    endforeach()

    set(solver TAO)
    foreach(model POWER_BALANCE_HIOP POWER_BALANCE_CARTESIAN)
      add_test(NAME OPFLOW_${solver}_${model}_${test_net_name}_FUNCTIONALITY_TEST
        COMMAND ${RUNCMD} $<TARGET_FILE:testOpflowFunctionality>
          -opflow_solver ${solver}
          -opflow_model ${model}
          -netfile ${test_net}
          -opflow_tolerance ${tolerance}
        )
    endforeach()

    set(solver IPOPT)
    foreach(model POWER_BALANCE_POLAR POWER_BALANCE_CARTESIAN CURRENT_BALANCE_CARTESIAN CURRENT_BALANCE_CARTESIAN2)
      add_test(NAME OPFLOW_${solver}_${model}_${test_net_name}_FUNCTIONALITY_TEST
        COMMAND ${RUNCMD} $<TARGET_FILE:testOpflowFunctionality>
          -opflow_solver ${solver}
          -opflow_model ${model}
          -netfile ${test_net}
          -opflow_tolerance ${tolerance}
        )
    endforeach()
  endforeach()

  message(STATUS "Configuring PFLOW functionality tests")
  foreach(np RANGE 1 2)
    foreach(test_net ${network_files})
      get_filename_component(input_file ${test_net} NAME)
      add_test(NAME PFLOW_${input_file}_${np}proc_FUNCTIONALITY_TEST COMMAND
              ${MPICMD} "-n" ${np} ${EXAGO_EXTRA_MPI_FLAGS}
              $<TARGET_FILE:testPflowFunctionality>
               -netfile ${test_net})
    endforeach(test_net)
  endforeach(np)

endif(EXAGO_RUN_TESTS)

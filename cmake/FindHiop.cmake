#[[

Finds HiOp include directory and libraries and exports target `HiOp`

User may set:
- HIOP_ROOT_DIR

]]

find_library(HIOP_LIBRARY
  NAMES
  hiop
  PATHS
  ${HIOP_DIR} $ENV{HIOP_DIR} ${HIOP_ROOT_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

if(HIOP_LIBRARY)
  get_filename_component(HIOP_LIBRARY_DIR ${HIOP_LIBRARY} DIRECTORY)
endif()

find_path(HIOP_INCLUDE_DIR
  NAMES
  hiopInterface.hpp
  PATHS
  ${HIOP_DIR} ${HIOP_ROOT_DIR} $ENV{HIOP_DIR} ${HIOP_LIBRARY_DIR}/..
  PATH_SUFFIXES
  include)

if(HIOP_LIBRARY AND HIOP_INCLUDE_DIR)
  message(STATUS "Found HiOp library: " ${HIOP_LIBRARY})
  message(STATUS "Found HiOp include directory: " ${HIOP_INCLUDE_DIR})
  add_library(HiOp INTERFACE)
  target_link_libraries(HiOp INTERFACE ${HIOP_LIBRARY})
  target_include_directories(HiOp INTERFACE ${HIOP_INCLUDE_DIR})
else()
  if(NOT HIOP_LIB)
    message(STATUS "HiOp library not found! Please provide correct filepath.")
  endif()
  if(NOT HIOP_INCLUDE_DIR)
    message(STATUS "HiOp include directory  not found! Please provide correct path.")
  endif()
endif()

set(HIOP_INCLUDE_DIR CACHE PATH "Path to HiOp header files")
set(HIOP_LIBRARY CACHE FILEPATH "Path to HiOp library")

# # If HIOP_DIR is not set then try to find it at default locations
# if(NOT HIOP_DIR)
#    find_path(HIOP_DIR 
#      NAMES 
#      _dist-default-build
#      lib/libhiop.a 
#      HINTS
#      /usr
#      /usr/local
#      /opt/local
#      ~/local/ipopt
#      # PNNL specific paths
#      /qfs/projects/exasgd/newell/hiop_shared_gpu
#      /qfs/projects/exasgd/newell/hiop_shared_cpu
#    )
# endif()

# # Exit if HIOP_DIR is not set and cannot be found in default locations
# if(NOT HIOP_DIR)
#    message(FATAL_ERROR "HiOP directory could not be found. Please specify the directory using the flag -DHIOP_DIR=<location_of_HiOP_install>")
# endif(NOT HIOP_DIR)

# # HiOp library location
# find_library(HIOP_LIBRARY NAME libhiop.so libhiop.dylib libhiop.a
#   HINTS 
#   ${HIOP_DIR}/lib
#   ${HIOP_DIR}/lib64
# )
# message(STATUS "HIOP library directory=${HIOP_LIBRARY}")

# # HiOP include directories
# # Find HiOP header path and ensure all needed files are there
# find_path(HIOP_INCLUDE_DIR NAME hiopNlpFormulation.hpp 
#   HINTS 
#   ${HIOP_DIR}/include
# )
# message(STATUS "HIOP Include directory = ${HIOP_INCLUDE_DIR}")

# include(CheckSymbolExists)

# if(HIOP_INCLUDE_DIR AND HIOP_LIBRARY)
#   include_directories(${HIOP_INCLUDE_DIR})
#   set(SCOPFLOW_HAVE_HIOP 1)
#   check_symbol_exists(HIOP_USE_GPU "${HIOP_INCLUDE_DIR}/hiop_defs.hpp" HIOP_USE_GPU)
#   if(HIOP_USE_GPU)
#     include(FindMagma)
#   endif()
# else()
#   message(FATAL_ERROR "HIOP not found!")
# endif()

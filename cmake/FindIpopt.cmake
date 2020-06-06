#[[

Finds Ipopt include directory and libraries and exports target `Ipopt`

User may set:
- IPOPT_ROOT_DIR

]]

find_library(IPOPT_LIBRARY
  NAMES
  ipopt
  PATHS
  ${IPOPT_DIR} $ENV{IPOPT_DIR} ${IPOPT_ROOT_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

if(IPOPT_LIBRARY)
  get_filename_component(IPOPT_LIBRARY_DIR ${IPOPT_LIBRARY} DIRECTORY)
endif()

find_path(IPOPT_INCLUDE_DIR
  NAMES
  IpTNLP.hpp
  PATHS
  ${IPOPT_DIR} ${IPOPT_ROOT_DIR} $ENV{IPOPT_DIR} ${IPOPT_LIBRARY_DIR}/..
  PATH_SUFFIXES
  include
  include/coin
  include/coin-or)

if(IPOPT_LIBRARY AND IPOPT_INCLUDE_DIR)
  message(STATUS "Found Ipopt library: " ${IPOPT_LIBRARY})
  message(STATUS "Found Ipopt include directory: " ${IPOPT_INCLUDE_DIR})
  add_library(Ipopt INTERFACE)
  target_link_libraries(Ipopt INTERFACE ${IPOPT_LIBRARY})
  target_include_directories(Ipopt INTERFACE ${IPOPT_INCLUDE_DIR})
else()
  if(NOT IPOPT_LIB)
    message(STATUS "Ipopt library not found! Please provide correct filepath.")
  endif()
  if(NOT IPOPT_INCLUDE_DIR)
    message(STATUS "Ipopt include directory  not found! Please provide correct path.")
  endif()
endif()

set(IPOPT_INCLUDE_DIR CACHE PATH "Path to Ipopt header files")
set(IPOPT_LIBRARY CACHE FILEPATH "Path to Ipopt library")
